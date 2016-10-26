library(matrixStats)      
library(edgeR)
library(WGCNA)
library(biomaRt)
library(plotrix)

library(foreach)
library(doMC)
registerDoMC()
library(proxy)
library(sgof)
library(multtest)
library(plyr)

getDoParWorkers()
options(cores=14)
getDoParWorkers()

# linux does not work with box so this code has a local data directory

setwd("/lawrencedata/ongoing_analyses/RNASeq023")

source("/lawrencedata/ongoing_analyses/RNASeq023/scripts/functionDefinitions.R")
try(dir.create("analysis/resultsCoexpr"), silent = F)
try( dir.create("analysis/figuresCoexpr"), silent = F)
try(dir.create("analysis/resultsCoSplicEx"), silent = F)
try( dir.create("analysis/figuresCoSplicEx"), silent = F)

# read raw data - can be improved by using read.table from original .txt file
geneReadsRaw=read.table("analysis/RNASeq023_mm10_coverage_splitoption.txt",sep="\t",header=F)

geneNames=geneReadsRaw[,4]

# Change column 4 name to "gene_sym"
names(geneReadsRaw)[4]<-"gene_sym"

# Combine the chromosome number, start location, and gene symbol to create a unique id  
# column for each exon
geneReadsRaw$exon<-paste(geneReadsRaw$V1,geneReadsRaw$V2,geneReadsRaw$V3,geneReadsRaw$gene_sym,sep="_")

# Create a data frame with gene symbol and exon read counts
exon_counts<-geneReadsRaw[,7:59]
exon_counts<-cbind(geneReadsRaw$gene_sym,exon_counts)

# Calculate the total counts for each gene for each sample
gene_counts<-ddply(exon_counts, 1, numcolwise(sum))

# Change the row names of the exon data frame to the exon unique ids (created above)
rownames(exon_counts)<-exon_counts$exon

# Remove the gene symbol column from the exon data frame
exon_counts<-exon_counts[,2:53]

# Change the row names of the gene data frame to the gene symbols
names(gene_counts)[1]<-"gene_sym"
rownames(gene_counts)<-gene_counts$gene_sym

# Remove the gene symbol column from the gene data frame
gene_counts<-gene_counts[,2:53]

# read sample info - sample names need to be inspected and categories extracted differently for each dataset !!!!!!!!!
sampleKey=read.table("data/sample_key_sorted.txt",sep="\t",header=F)
sampleKey[,1]=paste("S", sampleKey[,"V1"], sep = "")

countsSampleNames=read.table("data/bam_files_counts_header_order.txt", sep=",", header=F)

colnames(gene_counts) <- t(countsSampleNames[1,])

colnames(exon_counts) <- t(countsSampleNames[1,])

write.table(gene_counts,"RNASeq023_gene_reads_not_normalized.txt", sep="\t",quote=F,col.names=T,row.names=T)
write.table(exon_counts,"RNASeq023_exon_reads_not_normalized.txt", sep="\t",quote=F,col.names=T,row.names=T)
save.image("RNASeq023_Normalization.RData")

# calculate edgeR normalization factors and normalize the data - use all data
UQnormFactors_exons=calcNormFactors(exon_counts, method=c("upperquartile"))
UQnormFactors_genes=calcNormFactors(gene_counts, method=c("upperquartile"))

effectiveLibrarySizes_exons= UQnormFactors_exons*colSums(exon_counts)
effectiveLibrarySizes_genes= UQnormFactors_genes*colSums(gene_counts)

meanEffLibSize_exons=mean(effectiveLibrarySizes_exons)
meanEffLibSize_genes=mean(effectiveLibrarySizes_genes)

countNormFactor_exons= meanEffLibSize_exons/effectiveLibrarySizes_exons
countNormFactor_genes= meanEffLibSize_genes/effectiveLibrarySizes_genes

normalizedGeneCountsUQ_exons=0* exon_counts
normalizedGeneCountsUQ_genes=0* gene_counts

for (sample in 1:dim(normalizedGeneCountsUQ_exons)[2]){  
  normalizedGeneCountsUQ_exons[,sample]= exon_counts[, sample]* countNormFactor_exons[sample]  
}

for (sample in 1:dim(normalizedGeneCountsUQ_genes)[2]){  
  normalizedGeneCountsUQ_genes[,sample]= gene_counts[, sample]* countNormFactor_genes[sample]  
}

normalizedGeneCountsUQ_exons =round(normalizedGeneCountsUQ_exons)
normalizedGeneCountsUQ_genes =round(normalizedGeneCountsUQ_genes)

write.table(normalizedGeneCountsUQ_exons,"analysis/RNASeq023_exon_reads_UQNormalized.txt", sep="\t",quote=F,col.names=T,row.names=T)
write.table(normalizedGeneCountsUQ_genes,"analysis/RNASeq023_gene_reads_UQNormalized.txt", sep="\t",quote=F,col.names=T,row.names=T)

save.image("RNASeq023_Normalization.RData")

