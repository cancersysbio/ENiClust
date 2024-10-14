


outdir = "/oak/stanford/groups/ccurtis2/users/lisem/git_repo/ENiClust/eniclust/data/input/"

## 01
gene_cn_data_tcga = read.table("/oak/stanford/groups/ccurtis2/users/srinivap/Work/Projects/TCGA/TCGA_WES/pivot_gene_level_best_fit.txt", sep = '\t', header = T)
gene_cn_data_tcga$sample = sapply(strsplit(gene_cn_data_tcga$sample,"_"),"[[",1)

gene_cn_data_hmf = read.table("/oak/stanford/groups/ccurtis2/users/srinivap/Work/Projects/Hartwig/pivot_gene_level_best_fit.txt", sep = '\t', header = T)
gene_cn_data_hmf = gene_cn_data_hmf[which(grepl("_1800_1000_nd20",gene_cn_data_hmf$sample)),]
gene_cn_data_hmf$sample = sapply(strsplit(gene_cn_data_hmf$sample,"_"),"[[",1)
any(duplicated(gene_cn_data_hmf$sample))

dim(gene_cn_data_tcga)
dim(gene_cn_data_hmf)
length(intersect(colnames(gene_cn_data_tcga), colnames(gene_cn_data_hmf)))
gene_cn_data = dplyr::bind_rows(gene_cn_data_tcga, gene_cn_data_hmf)
dim(gene_cn_data)

colnames(gene_cn_data)[1] = "Sample"
head(gene_cn_data[,'Sample'])
tail(gene_cn_data[,'Sample'])

all(mysamples %in% receptors_file$Sample)
gene_cn_data$Sample[which(!gene_cn_data$Sample %in% receptors_file$Sample)]
unique(segments_file$Sample)[which(!unique(segments_file$Sample) %in% receptors_file$Sample)]
qc_file$Sample[which(!qc_file$Sample %in% receptors_file$Sample)]

## 02
segments_file_tcga = read.table("/oak/stanford/groups/ccurtis2/users/srinivap/Work/Projects/TCGA/TCGA_WES/merged_cncf_best_fit.txt", header = T, sep ="\t")
segments_file_tcga$ID = sapply(strsplit(segments_file_tcga$ID,"_"),"[[",1)
segments_file_hmf = read.table("/oak/stanford/groups/ccurtis2/users/srinivap/Work/Projects/Hartwig/merged_cncf_best_fit.txt", header = T, sep ="\t")
segments_file_hmf = segments_file_hmf[which(grepl("_1800_1000_nd20",segments_file_hmf$ID)),]
segments_file_hmf$ID = sapply(strsplit(segments_file_hmf$ID,"_"),"[[",1)

all(colnames(segments_file_tcga) == colnames(segments_file_hmf))
segments_file = rbind(segments_file_tcga, segments_file_hmf)
segments_file$Sample = segments_file$ID
head(segments_file[,'Sample'])
tail(segments_file[,'Sample'])

## 03
qc_file_tcga = read.table("/oak/stanford/groups/ccurtis2/users/lisem/eniclust/wes_wgs/data/TCGA_WES/Facets_qc_summary_best_fit_hrd_corrected.txt", header = T, sep ="\t")
qc_file_hmf = read.table("/oak/stanford/groups/ccurtis2/users/lisem/eniclust/wes_wgs/data/Hartwig/Facets_qc_summary_best_fit_hrd_corrected.txt", header = T, sep ="\t")
qc_file_hmf = qc_file_hmf[which(grepl("_1800_1000_nd20",qc_file_hmf$sample)),]
qc_file_hmf$sample = sapply(strsplit(qc_file_hmf$sample,"_"),"[[",1)
any(duplicated(qc_file_hmf$sample))

all(colnames(qc_file_tcga) == colnames(qc_file_hmf))
qc_file = rbind(qc_file_tcga, qc_file_hmf)
qc_file[,'Sample'] = qc_file$sample
qc_file[,c('Sample', 'genome_doubled', 'ploidy', 'hrd_loh', 'fraction_cna')]
head(qc_file[,'Sample'])
tail(qc_file[,'Sample'])

## 05
receptors_file_tcga = read.table("/oak/stanford/groups/ccurtis2/users/srinivap/Work/Projects/TCGA/samples_labels.txt", sep = ",", header = T)
any(duplicated(str_sub(receptors_file_tcga$ID,1,12)))
receptors_file_tcga = data.frame(Sample = str_sub(receptors_file_tcga$ID,1,12), ER = receptors_file_tcga$ER, HER2 = receptors_file_tcga$HER2, stringsAsFactors = F)
receptors_file_tcga = receptors_file_tcga[which(!duplicated(receptors_file_tcga$Sample)),]

receptors_file_hmf = read.table("/oak/stanford/groups/ccurtis2/users/srinivap/Work/Projects/Hartwig/metadata.tsv", sep = "\t", header = T)
receptors_file_hmf = data.frame(Sample = receptors_file_hmf$sampleId, ER = ifelse(grepl("ER-positive",receptors_file_hmf$cancerSubtype),1,0), HER2 = ifelse(grepl("HER2-positive",receptors_file_hmf$cancerSubtype),1,0), stringsAsFactors = F)

receptors_file = rbind(receptors_file_tcga, receptors_file_hmf)
head(receptors_file[,'Sample'])
tail(receptors_file[,'Sample'])

## iC10 reference
tcga_ref = read.table("/oak/stanford/groups/ccurtis2/users/lisem/TCGA/iC10/tcga_iC10_labels.txt", sep = ' ', header = T)
hmf_ref = read.table("/oak/stanford/groups/ccurtis2/users/lisem/Hartwig/iC10/hmf_iC10_labels.txt", sep = ' ', header = T)

any(duplicated(tcga_ref$Sample))
any(duplicated(hmf_ref$Sample))

common_samples = unique(intersect(intersect(intersect(intersect(gene_cn_data$Sample, segments_file$Sample), qc_file$Sample), receptors_file$Sample), mutations_file$Sample))
length(common_samples)

mysamples = common_samples[which(common_samples %in% c(tcga_ref$Sample,hmf_ref$Sample))]
length(mysamples)

reference_file = data.frame(Sample = mysamples, ic10 = NA, stringsAsFactors = F)
reference_file$ic10 = sapply(reference_file$Sample, function(x) if(x %in% tcga_ref$Sample){tcga_ref$iC10[which(tcga_ref$Sample == x)]}else{hmf_ref$iC10[which(hmf_ref$Sample == x)]})
any(is.na(reference_file$ic10))

## Loading set of sample
training_testing = read.table('/oak/stanford/groups/ccurtis2/users/lisem/eniclust/wes_wgs/dev/test/rand/rand_59_cv_noLOW_v3/Pipeline_2/training_and_testing_set.txt', header = T)
validation = read.table('/oak/stanford/groups/ccurtis2/users/lisem/eniclust/wes_wgs/dev/test/rand/rand_59_cv_noLOW_v3/Pipeline_2/validation_set.txt', header = T)

mysamples = c(training_testing$Sample, validation$Sample)
length(mysamples[which(grepl("TCGA",mysamples))])
length(mysamples[which(!grepl("TCGA",mysamples))])

## Re-doing reference
reference_file = data.frame(Sample = mysamples, ic10 = NA, stringsAsFactors = F)
reference_file$ic10 = sapply(reference_file$Sample, function(x) if(x %in% tcga_ref$Sample){tcga_ref$iC10[which(tcga_ref$Sample == x)]}else{hmf_ref$iC10[which(hmf_ref$Sample == x)]})
any(is.na(reference_file$ic10))

## Saving: Missing some TCGA samples because of receptors (n=53)
data = receptors_file
length(data$Sample[which(data$Sample %in% mysamples)][which(grepl("TCGA",data$Sample[which(data$Sample %in% mysamples)]))])
length(data$Sample[which(data$Sample %in% mysamples)][which(!grepl("TCGA",data$Sample[which(data$Sample %in% mysamples)]))])

length(unique(segments_file$Sample)[which(unique(segments_file$Sample) %in% mysamples)][which(grepl("TCGA",unique(segments_file$Sample)[which(unique(segments_file$Sample) %in% mysamples)]))])
length(unique(segments_file$Sample)[which(unique(segments_file$Sample) %in% mysamples)][which(!grepl("TCGA",unique(segments_file$Sample)[which(unique(segments_file$Sample) %in% mysamples)]))])

write.table(gene_cn_data[which(gene_cn_data$Sample %in% mysamples),], paste(outdir, "01_gene_level.txt", sep=""), sep="\t", row.names=F, quote=F)
write.table(segments_file[which(segments_file$Sample %in% mysamples),], paste(outdir, "02_segments.txt", sep=""), sep="\t", row.names=F, quote=F)
write.table(qc_file[which(qc_file$Sample %in% mysamples),], paste(outdir, "03_qc.txt", sep=""), sep="\t", row.names=F, quote=F)
write.table(receptors_file[which(receptors_file$Sample %in% mysamples),], paste(outdir, "05_receptors.txt", sep=""), sep="\t", row.names=F, quote=F)
write.table(reference_file[which(reference_file$Sample %in% mysamples),], paste(outdir, "ic10_references.txt", sep=""), sep="\t", row.names=F, quote=F)
