#!/bin/bash

set -e
export LC_NUMERIC=C

BPP="./BPP/"
BPPdemo="${BPP}/demo"
threads=6    #$(fgrep -c processor /proc/cpuinfo) # Use number of threads equal to CPU cores

genome=$1 # Genome file - open
gtf=$2 # Genome annotation gzipped
vcf=$3 # sample vcf file bgzipped
sample=${vcf%.vcf.gz}

# === Downloading annotation files ===========
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.primary_assembly.genome.fa.gz
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz

# Downloading VCF
# wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
# wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

# Downloading ClinVar
# wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
# wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
# --------------------------------------------


# === Service tasks for development ==========

# BGZIP genome
# zcat GRCh38.primary_assembly.genome.fa.gz | bgzip > GRCh38.fa.gz

# Extract one chromosome from file
# samtools faidx GRCh38.primary_assembly.genome.fa.gz chr1 > GRCh38.chr1.fa

# Subsetting gtf annotation
# zcat gencode.v26.annotation.gtf.gz | grep -P '^chr1\t' > GRCh38.chr1.gtf

# Make a 1-st chromosome subset from vcf file
# zcat HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | grep -w '^#\|^chr1' | bgzip > HG001_GRCh38.chr1.vcf.gz

# Make a small subset from vcf file
# zcat homo_sapiens-chr1.vcf.gz | head -n 100000 > test.vcf

# ============================================

# --- Preparing reference data - need to do only once ---

# Filtering gtf to have only transcripts and exons
echo "Preparing annotations"
filtered_gtf=${gtf%.gtf.gz}.filtered.gtf.gz
[[ ! -f ${filtered_gtf} ]] && zcat ${gtf} | grep -v '#' | awk '$3=="transcript" || $3=="exon"' | bgzip > ${filtered_gtf}

# Sorting preserving order of transcripts and exons
sorted_gtf=${gtf%.gtf.gz}.sorted.gtf.gz
[[ ! -f ${sorted_gtf} ]] && zcat ${filtered_gtf} | sed -e 's/gene/1gene/' -e 's/transcript/2transcript/' -e 's/exon/3exon/' \
    | sort -k1,1 -k4,4n -k3,3 | sed -e 's/1gene/gene/' -e 's/2transcript/transcript/' -e 's/3exon/exon/' | bgzip > ${sorted_gtf}

# Obtaining introns from exons
intrones_gtf=${gtf%.gtf.gz}.intrones.gtf.gz
[[ ! -f ${intrones_gtf} ]] && zcat ${sorted_gtf} | ./get_introns.py | bgzip > ${intrones_gtf}

# Making index for reference genome if not exists
echo "Making reference intrones"
[[ ! -f ${genome}.fai ]] && samtools faidx ${genome}

# Obtaining reference intrones sequences
intrones_fasta_ref=${genome%.fa}.intrones.ref.fa
[[ ! -f ${intrones_fasta_ref} ]] && bedtools getfasta -fi ${genome} -bed ${intrones_gtf} -name+ -s > ${intrones_fasta_ref}


# --- Applying variations - need to do for every new vcf ---


# Converting VCF to branchpointer table and splitting to chunks 
echo "Splitting varшations to chunks"
bp_prefix=${sample}.branchpointer.part_
[[ ! -f ${bp_prefix}aaaa.tsv ]] && echo zcat ${vcf} | grep -v '#' | awk -v FS="\t" -v OFS="\t" '{ if ($3 == ".") ($3 = sprintf("%s%s","varId",NR)); print $3,$1,$2,".",$4,$5}' \
  | split -a 4 -l ${chunk_size} --additional-suffix=".tsv" - ${bp_prefix}


# "Saving run commands to runlist"
echo "Creating a run list for parallel run"
run_str_inject="./bpp_inject_var.py ${intrones_fasta_ref}" # Stop code 255 to exit xargs
run_str_bpp="python2 BPP/BP_PPT.py -b ${BPPdemo}/pwmBP_human.txt -p ${BPPdemo}/scPPT_human.txt"

echo "#!/bin/bash" > runlist.sh
for file in `ls ${bp_prefix}????_??.tsv`; do
    file_bpp=${file%.tsv}.bpp.tsv
    if [ ! -f ${file_bpp} ]; then
        file_fa=${file%.tsv}.intrones.fa
        echo "${run_str_inject} ${file} > ${file_fa} 2> ${file_fa}.log && ${run_str_bpp} -i ${file_fa} > ${file_bpp} &&  rm ${file_fa}" >> runlist.sh
    fi
done

# Running tasks for runlist in parallel
echo "Running tasks for runlist in parallel"
nice xargs --arg-file=runlist.sh --max-procs=$threads --replace /bin/sh -c "{}"

# Combining data
echo "Combining data"
table_combined=${sample}.${genome}.bpp.combined_per_variation.tsv.gz
[[ ! -f ${table_combined} ]] && cat ${bp_prefix}*.bpp.tsv | ./bpp_collect_var.py | gzip > ${table_combined}

echo "Vacuuming final data"
table_combined_vacuumed=${sample}.${genome}.bpp.combined_per_variation.vacuumed.tsv
[[ ! -f ${table_combined_vacuumed} ]] && zcat ${table_combined} | cut -f 1-5,8,14,16,22 | awk -v FS="\t" -v OFS="\t" '($6!=$8) || ($7!=$9)' > ${table_combined_vacuumed}

