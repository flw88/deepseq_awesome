#!/usr/bin/Rscript



## setwd and load files
print("Loading files...")
setwd("/home/local/users/eflynn/deepseq/")
fpkm<-read.csv("datasets/scRNASeq/fpkm_table.csv")
pheno<-read.csv("datasets/scRNASeq/columns-cells.csv")
genes<-read.csv("datasets/scRNASeq/rows-genes.csv")
het_genes<-read.table("results/DE/union_of_batches_het_genes.txt")

## reformat fpkm df
fpkm[,1] -> row.names(fpkm)		     ## rename rows
fpkm[,-1] -> fpkm			       ## remove gene id column
gsub("^X","",names(fpkm)) -> names(fpkm)	## remove leading X from cell IDs
#fpkm <- t(fpkm)				## transform so rows are cells, columns are genes (nvm, do later)

## Compute PCs of all genes, write to file
print("Computing principal components on all genes")
#fpkm_pr<-prcomp(t(fpkm),center=T,scale.=F)
fpkm_pr<-prcomp(fpkm,center=T,scale.=F)
unclass(fpkm_pr$rotation) -> fpkm_pr_rot
fpkm_pr_rot[1:5,1:5]
write.table(fpkm_pr_rot,file="results/batch_effect/fpkm_pr.txt",quote=FALSE,
	sep="\t", row.names=TRUE, col.names=TRUE)


#print("Computing principal components on genes fpkm > 0")
#
### Isolate genes w/ fpkm > 0
#row_sub = apply(fpkm, 1, function(row) all(row!=0))
#fpkm[row_sub,] -> sub_fpkm
#
#
### Compute components and plot
##sub_fpkm_pr<-prcomp(t(sub_fpkm),center=T,scale.=F)
#sub_fpkm_pr<-prcomp(sub_fpkm,center=T,scale.=F)
#unclass(sub_fpkm_pr$rotation) -> sub_fpkm_pr_rot
#sub_fpkm_pr_rot[1:5,1:5]
#write.table(sub_fpkm_pr_rot,file="results/batch_effect/sub_fpkm_pr.txt",quote=FALSE,
#	sep="\t", row.names=TRUE, col.names=TRUE)

print("Computing principal components on het genes")

## Isolate genes
rownames(fpkm) -> fpkm$gene_id
subset(fpkm, gene_in %in% het_genes$V1) -> het_fpkm
het_fpkm[,-ncol(het_fpkm)] -> het_fpkm

## Compute PCs, write to file
fpkm_het_pr<-prcomp(het_fpkm,center=T,scale.=F)
unclass(fpkm_het_pr$rotation) -> fpkm_het_pr_rot
fpkm_het_pr_rot[1:5,1:5]
write.table(fpkm_het_pr_rot,file="results/batch_effect/fpkm_hetgenes_pr.txt",quote=FALSE,
        sep="\t", row.names=TRUE, col.names=TRUE)


