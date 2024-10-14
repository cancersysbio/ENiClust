#!/usr/bin/env Rscript

#######################################
#   Run ENiClust script for Isabl     #
#   Authors: Lise and Aziz            #
#######################################

library("dplyr")
library("stringr")
suppressMessages(library("argparse", quietly = TRUE))

rm(list=ls())
gc()

parser <- ArgumentParser(description="Prepare inputs for ENiClust")
parser$add_argument(
    "--outdir",
    required=TRUE,
    help="Output Directory")

parser$add_argument(
    "--sample",
    required=TRUE,
    help="Sample name or ID")

parser$add_argument(
    "--segment",
    required=TRUE,
    help="FACETS segments hisens.cncf.txt file")

parser$add_argument(
    "--gene_level",
    required=TRUE,
    help="FACETS gene level calls *.gene_level.txt")

parser$add_argument(
    "--qc",
    required=TRUE,
    help="QC text file")

parser$add_argument(
    "--hg_build",
    required=TRUE,
    default="hg38",
    help="Reference genome (hg38 or hg19)")

parser$add_argument(
    "--ER",
    required=FALSE,
    default=1,
    help="ER status (0 if -ive or 1 for +)")

parser$add_argument(
    "--HER2",
    required=FALSE,
    default=0,
    help="HER2 Status (0 if if -ive or 1 for if +)")

opt <- parser$parse_args()

setwd(opt$outdir) # run in outdir
set.seed(1234) # for reproducibility

if (opt$hg_build == FALSE){
  hg_build = "hg38"
}

if (opt$outdir == ""){
  opt$outdir = "."
}

## get the arguments
hg_build = opt$hg_build
output_dir = paste0(opt$outdir,"/")

print(paste("Runing the script using",hg_build, "annotations. Outputs will be stored in",output_dir))

# First file contains gene level 
gene_cn_data <- read.table(opt$gene_level, sep="\t", header=T)
    
gene_cn_data <- gene_cn_data %>% dplyr::select("gene", "tcn.em")
rownames(gene_cn_data) = gene_cn_data$gene
colnames(gene_cn_data) <- c("gene", opt$sample)

rownames(gene_cn_data) <- gene_cn_data$gene
gene_cn_data$gene <- NULL
gene_cn_data <- as.data.frame(t(gene_cn_data))
gene_cn_data$Sample <- rownames(gene_cn_data)

dim(gene_cn_data)

gene_cn_data["Sample"] = rownames(gene_cn_data)
gene_cn_data = gene_cn_data[,c(ncol(gene_cn_data), 1:(ncol(gene_cn_data)-1))]

# Second file contains segment information 
segment_data <- read.table(opt$segment, header=T, sep="\t")

# Third file contains qc information
qc_data <- read.table(opt$qc, header=T, sep="\t")
qc_data <- qc_data[which(qc_data$run_type == "hisens"),]

colnames(segment_data)[which(colnames(segment_data) == "ID")] = "Sample"
segment_data$Sample = opt$sample
head(segment_data)

colnames(qc_data)[which(colnames(qc_data) == "sample")] = "Sample"
head(qc_data)

write.table(gene_cn_data, paste(output_dir, "01_gene_level.txt", sep=""), sep="\t", row.names=F, quote=F)
write.table(segment_data, paste(output_dir, "02_segments.txt", sep=""), sep="\t", row.names=F, quote=F)
write.table(qc_data, paste(output_dir, "03_qc.txt", sep=""), sep="\t", row.names=F, quote=F)

##-----------------------------------------------------------##
# Match clinical metadata by submitted sample ids
##-----------------------------------------------------------##

receptor_data = data.frame(Sample = opt$sample,
                           ER = opt$ER, 
                           HER2 = opt$HER2, stringsAsFactors = F)

head(receptor_data)
write.table(receptor_data, paste(output_dir, "05_receptors.txt", sep=""), sep="\t", row.names=F, quote=F)
