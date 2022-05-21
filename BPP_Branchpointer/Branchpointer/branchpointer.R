#!/usr/bin/Rscript
# The script to run BRANCHPOINTER

# installation

# apt install libssl-dev libxml2-dev libhts-dev #libcurl4-openssl-dev

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("BSgenome")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("branchpointer")

#if (length(args)!=2) {
#  stop("Please provide vcf file for input", call.=FALSE)  #TODO: fix
#}
args = commandArgs(trailingOnly=TRUE)

gft_file_name <- "gencode.v26.annotation.gtf"
transcript_names <- args[1]
querySNPFile <- args[2]


message("Loading libraries, this can take a while...")
suppressPackageStartupMessages(library(branchpointer))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))

message("Loading reference genome")
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

message("Loading transcripts")
transcripts <- scan(transcript_names, character())

message("Loading annotations")
exons <- gtfToExons(gft_file_name)

message("Creating intrones from GTF")
queryIntronFromGTF <- makeBranchpointWindowForExons(transcripts, idType = "transcript_id", exons = exons)

#message("Running branchpointer on reference genome")
#branchpointPredictionsIntron <- predictBranchpoints(queryIntronFromGTF, queryType = "region", BSgenome = g)
#write.csv(branchpointPredictionsIntron,"hg38.branchpointer.ref.csv", row.names = FALSE)

message("Loading variations from: ", querySNPFile)
querySNP <- readQueryFile(querySNPFile, queryType = "SNP", exons = exons, filter = TRUE)

message("Running BrabchPointer on variations")
branchpointPredictionsSNP <- predictBranchpoints(querySNP, queryType = "SNP", BSgenome = g)
#outname <- paste0("hg38.branchpointer.", querySNPFile, ".csv")
#write.csv(branchpointPredictionsSNP, outname, row.names = FALSE)

message("Calculating summary")
querySNPSummary <- predictionsToSummary(querySNP, branchpointPredictionsSNP)
outname2 <- paste0("hg38.branchpointer.", querySNPFile, ".summary_original.csv")
write.csv(querySNPSummary, outname2, row.names = FALSE)
