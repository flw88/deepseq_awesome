#!/usr/bin/Rscript


## setwd and load files
print("Loading files...")
setwd("/Users/eliseflynn/Google Drive/Deep Seq/DeepSeq_Project")
gtf<-read.table("Data/rsem_GRCm38.p3.gtf",
	header=FALSE,sep='\t')
#fpkm<-read.csv("datasets/scRNASeq/fpkm_table.csv")
#pheno<-read.csv("datasets/scRNASeq/columns-cells.csv")
#genes<-read.csv("datasets/scRNASeq/rows-genes.csv")


## reformat last column of gtf
gtf_tr <- data.frame(do.call(rbind, strsplit(as.character(gtf[,9]),"; ")))
names(gtf_tr) <- c("gene_id","gene_symbol","transcript_id")
gtf_tr$gene_id <- gsub('gene_id ','',gtf_tr$gene_id)
gtf_tr$gene_symbol <- gsub('gene_symbol ','',gtf_tr$gene_symbol)
gtf_tr$transcript_id <- gsub('transcript_id ','',gtf_tr$transcript_id)

gtf_tr_unique <- unique(gtf_tr)

subset(gtf_tr_unique,grepl("^[XM][A-Z]_",transcript_id)) -> gtf_NM
write.table(gtf_NM,file="Data/rsem_GRCm38.p3.transcript.txt",
	sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)