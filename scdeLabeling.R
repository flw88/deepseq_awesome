#!/usr/bin/Rscript

## setwd and load files
print("Loading files...")

args=commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
	stop("This script requires one argument (the pca input file)")
} else {
	input_DE = read.table(args[1],header=TRUE)
	input_DE_name = basename(args[1])
}

setwd("/cygdrive/c/Users/eflyn/Google Drive/Deep Seq/DeepSeq_Project")
refseq<-read.table("Data/RefSeq/RefSeq_Genes_unique.txt",
	header=TRUE,sep='\t')
#fpkm<-read.csv("datasets/scRNASeq/fpkm_table.csv")
#pheno<-read.csv("datasets/scRNASeq/columns-cells.csv")
#genes<-read.csv("datasets/scRNASeq/rows-genes.csv")


merge(input_DE,refseq,by.x='gene_symbol',by.y='name2') -> DE_labeled
write.table(DE_labeled,
	file=paste("Results/scde/",gsub('tsv','refseq.tsv',input_DE_name),sep=''),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')


colnames(DE_labeled)[11] <- "gID"
DE_labeled[c("gID","Z")] -> Z_table
write.table(Z_table,
	file=paste("Results/iGET_inputfiles/",
		gsub('tsv','refseq',input_DE_name),
		".Z.txt",
		sep=''),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

Z_table -> Z_table_bin
Z_table_bin$Z_bin <- 0
Z_table_bin[row.names(subset(Z_table_bin,
	Z > -2 & Z <= 0)),]$Z_bin <- 1
Z_table_bin[row.names(subset(Z_table_bin,
	Z > 0 & Z < 2)),]$Z_bin <- 2
Z_table_bin[row.names(subset(Z_table_bin,
	Z >= 2 & Z < 3)),]$Z_bin <- 3
Z_table_bin[row.names(subset(Z_table_bin,
	Z >= 3 & Z < 4)),]$Z_bin <- 4
Z_table_bin[row.names(subset(Z_table_bin,
	Z >= 4)),]$Z_bin <- 5
write.table(Z_table_bin[,-2],
	file=paste("Results/iGET_inputfiles/",
		gsub('tsv','refseq',input_DE_name),
		".Z_bin.txt",
		sep=''),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')


