#!/usr/bin/Rscript


## setwd and load files
print("Loading files...")
setwd("/Users/eliseflynn/Google Drive/Deep Seq/DeepSeq_Project")
refseq<-read.table("Data/RefSeq/RefSeq_Genes.txt",
	header=TRUE,sep='\t')
#fpkm<-read.csv("datasets/scRNASeq/fpkm_table.csv")
#pheno<-read.csv("datasets/scRNASeq/columns-cells.csv")
#genes<-read.csv("datasets/scRNASeq/rows-genes.csv")


## combine refseq table to one NM (smallest) per gene
refseq_NM <- refseq[c("name","name2")]
unique(refseq_NM) -> refseq_NM
refseq_NM$name.integer <- as.integer(gsub('N[A-Z]_','',refseq_NM$name))
aggregate(refseq_NM$name.integer, 
	by=list(refseq_NM$name2), 
	min) -> refseq_collapsed


## add additional info to NM table
merge(refseq_collapsed, refseq_NM, 
	by.x = c("Group.1","x"), 
	by.y=c("name2","name.integer")) -> refseq_collapsed_labeled
names(refseq_collapsed_labeled) <- c('name2','name.integer','name')

refseq[!duplicated(refseq$name),] -> refseq_unique
merge(refseq_collapsed_labeled,refseq_unique, 
	by=c("name","name2"))[,-3] -> refseq_collapsed_full


## write final table to file
write.table(refseq_collapsed_full,
	file="Data/RefSeq/RefSeq_Genes_unique.txt",
	col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)