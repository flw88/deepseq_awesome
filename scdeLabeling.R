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


colnames(DE_labeled)[12] <- "gID"
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



######################## With gtf multiple transcripts #########################

setwd("/Users/eliseflynn/Google Drive/Deep Seq/DeepSeq_Project")
gtf_NM <- read.table(file="Data/rsem_GRCm38.p3.NMtranscript.txt",
	sep='\t',header=TRUE)

merge(input_DE,gtf_NM,by="gene_symbol") -> DE_labeled
write.table(DE_labeled,
	file=paste("Results/scde/",gsub('tsv','gtf.tsv',input_DE_name),sep=''),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

colnames(DE_labeled)[12] <- "gID"
DE_labeled[c("gID","Z")] -> Z_table
write.table(Z_table,
	file=paste("Results/iGET_inputfiles/",
		gsub('tsv','gtf',input_DE_name),
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
		gsub('tsv','gtf',input_DE_name),
		".Z_bin.txt",
		sep=''),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')



subset(Z_table_bin, Z<= -2 | Z>= 2) -> Z_table_bin_red
subset(Z_table_bin, Z> -2 & Z< 2) -> uninteresting
set.seed(88)
rbind(uninteresting[sample(nrow(uninteresting),10000), ],
	Z_table_bin_red) -> Z_table_bin_red
write.table(Z_table_bin_red[,-2],
	file=paste("Results/iGET_inputfiles/",
		gsub('tsv','gtf',input_DE_name),
		".Z_bin.reduced.txt",
		sep=''),
	col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')



