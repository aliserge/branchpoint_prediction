#!/bin/bash

set -e
export LC_NUMERIC=C

vcf=${1}

gft_file_name="gencode.v26.annotation.gtf.gz"
transcript_names="gencode.v26.annotation.transcriptID"
threads=4 #$(fgrep -c processor /proc/cpuinfo) # Use number of threads equal to CPU cores
chunk_size=10000


# 
# A script to prepare data to Branchpointer format and run branchpointer

# ===== Reference genome loading from file - Not used ================
### 1. Loading reference genome
# [[ ! -f "gencode.v26.annotation.gtf.gz" ]] && wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz

### 3. дописываем все поля в новом файлике, чтобы программка не ругалась:
### cat GRCh38_latest_genomic.gtf | sed -r 's/gene_biotype ("[A-z_]+")/gene_biotype \1; transcript_biotype \1/g' > GRCh38_latest_genomic_add-biotype.gtf
# ======================================================================

# Creating a file with transcript names
[[ ! -f ${transcript_names}.txt ]] && zcat ${gft_file_name} | cut -f 9 | grep -Po 'transcript_id "[A-z_0-9\.\-]+"' | cut -d ' ' -f 2 | sed 's/"//g' | uniq > ${transcript_names}.txt

# Converting VCF to branchpointer table and splitting to chunks 
bp_prefix=${vcf%.vcf.gz}.branchpointer.part_
[[ ! -f ${bp_prefix}aaaa.tsv ]] && zcat ${vcf} | grep -v '#' | awk -v FS="\t" -v OFS="\t" '{ if ($3 == ".") ($3 = sprintf("%s%s","varId",NR)); print $3,$1,$2,".",$4,$5}' \
  | split -a 4 -l ${chunk_size} --additional-suffix=".tsv" - ${bp_prefix}


# "Saving run commands to runlist"
run_str="./branchpointer.R ${transcript_names}.txt" # Stop code 255 to exit xargs
echo "#!/bin/bash" > runlist.sh
for file in `ls ${bp_prefix}????.tsv`; do
	if [ ! -f "hg38.branchpointer.${file}.summary_original.csv" ]; then
		echo "${run_str} ${file}" >> runlist.sh
	fi
done

# Running tasks for runlist in parallel
nice xargs --arg-file=runlist.sh --max-procs=$threads --replace /bin/sh -c "{}"


# Collecting summary data
(head -1 hg38.branchpointer.${bp_prefix}aaaa.tsv.summary_original.csv ; tail -n +3 -q hg38.branchpointer.${bp_prefix}????.tsv.summary_original.csv ) > hg38.branchpointer.${vcf%.vcf.gz}.summary.combined.csv

