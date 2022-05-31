# Branchpoint effect on variations

The software below calculates the effect if intronic variations
on the position and strength of branchpoint.

Scripts use two existing applications: branchpointer and BPP:
* https://github.com/zhqingit/BPP
* https://www.bioconductor.org/packages/release/bioc/html/branchpointer.html

## Licence
This software is licensed under GNU GPLv3 licence: https://www.gnu.org/licenses/gpl.html

# Installation for Ubuntu 20.04 and 22.04
All scripts here uses Python language version 3 and bash shell. 
The BioPython package also should be installed for python3.
The BPP uses the version 2 of python, which also should be installed.

## BPP scripts dependencies

    sudo apt install python2 python3-pip samtools tabix bedtools 
    pip3 install biopython

## Branchpointer dependencies
    apt install build-essential libssl-dev libxml2-dev libhts-dev 

For ubuntu 20.04 also install development files for libcurl:
`sudo apt install r-base libcurl4-openssl-dev`

## Branchpointer R packages
Run `R` and execute the following commands:

    install.packages("BiocManager")
    BiocManager::install("BSgenome")
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
    BiocManager::install("branchpointer")

# Data dependencies
The software was developed and tested using the following data files:
* Reference genome GRCh38: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.primary_assembly.genome.fa.gz
* Transcripts annotation gencode v2.6: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz

# Usage

## Using bpp
`./bpp_var.sh <reference.fa> <annotation.gtf.gz> <vcf.gz>`
* reference.fa - reference genome, fasta format
* annotation.gtf.gz - transcripts for the reference genome, gtf format, gzipped
* vcf.gz - variations to analyze in VCF format, gzipped

### Example
`./bpp_var.sh GRCh38.fa GRCh38.gtf.gz clinvar.pathogenic.vcf`

## Using branchpointer
`branchpointer.sh <vcf.gz>`
* vcf.gz - variations to analyze in VCF format, gzipped

The branchpointer script has hardcoded dependency for anootation file
`gencode.v26.annotation.gtf`. It should be located in the run folder of 
the `branchpointer.sh` script.

### Example
`./branchpointer.sh clinvar.pathogenic.vcf.gz`

# Detailed description of data processing

## BPP
The original BPP software takes the sequences of introns in plain fasta format 
and predicts the branchpoint location and score in the intron. 
The main script `bpp_var.sh` and set of additional python scripts are used to 
prepare data for BPP run and collecting results. 

Briefly, the `bpp_var.sh` do the following

* Prepare annotation data, sort, filter
* Create introns form exons using `get_introns.py` script
* Create reference intronic sequences in FASTA format
* Convert variations to table format and split it to chunks (useful for very big 
VCF files like GNOMAD). The table format is the same for BPP wrapper scripts and branchpointer, so the same 
split files can be used for both utilities.
* Execute in parallel for each vcf chunk the following actions
    * Create alternative introns with injected variations using `bpp_inject_var.py` script
    * Run BPP on alternative variations
* Collect data from BPP results to per-variation table using `bpp_collect_var.py` script
* Filter data and create summary table

## Branchpointer
The branchpointer library has the build-in ability to process variations file, so 
the wrapper script for branchpointer is much easier 

The main script `branchpointer.sh` does the follows:

* Extract all transcript names form annotation. One can edit this file to
* Convert variations to table format and split it to chunks (useful for very big 
VCF files like GNOMAD). The table format is the same for BPP wrapper scripts and branchpointer, so the same 
split files can be used for both utilities.
* Execute in parallel for each vcf chunk the `branchpointer.R` script.
The code of `branchpointer.R` just loads data, runs functions from 
branchpointer library according to manual and saves final data
* Collect data from branchpointer results to summary table


## References
[1] Leman R, Tubeuf H, Raad S, Tournier I, Derambure C, Lanos R, Gaildrat P, Castelain G, Hauchard J, Killian A, Baert-Desurmont S, Legros A, Goardon N, Quesnelle C, Ricou A, Castera L, Vaur D, Le Gac G, Ka C, Fichou Y, Bonnet-Dorion F, Sevenet N, Guillaud-Bataille M, Boutry-Kryza N, Schultz I, Caux-Moncoutier V, Rossing M, Walker LC, Spurdle AB, Houdayer C, Martins A, Krieger S. Assessment of branch point prediction tools to predict physiological branch points and their alteration by variants. BMC Genomics. 2020 Jan 28;21(1):86. doi: 10.1186/s12864-020-6484-5. PMID: 31992191; PMCID: PMC6988378.

[2] Zhang Q, Fan X, Wang Y, Sun MA, Shao J, Guo D. BPP: a sequence-based algorithm for branch point prediction. Bioinformatics. 2017 Oct 15;33(20):3166-3172. doi: 10.1093/bioinformatics/btx401. PMID: 28633445.

[3] Signal B, Gloss B, Dinger M, Mercer T (2016). “Machine-learning annotation of human splicing branchpoints.” bioRxiv. doi: 10.1101/094003, http://biorxiv.org/content/early/2016/12/14/094003.


