#!/usr/bin/Rscript



## setwd and load files
print("Loading files...")
setwd("/home/local/users/eflynn/deepseq/")

fpkm<-read.csv("datasets/scRNASeq/fpkm_table.csv")
pheno<-read.csv("datasets/scRNASeq/columns-cells.csv")
#genes<-read.csv("datasets/scRNASeq/rows-genes.csv")

## reformat fpkm df
fpkm[,1] -> row.names(fpkm)                     ## rename rows
fpkm[,-1] -> fpkm                               ## remove gene id column
gsub("^X","",names(fpkm)) -> names(fpkm)        ## remove leading X from cell IDs
fpkm_t <- data.frame(t(fpkm)                    ## transform so rows are cells, columns are genes

names(fpkm_t) <- gsub("^X","",names(fpkm_t))
rownames(fpkm_t)-> fpkm_t$rnaseq_profile_id


## create Cre line data tables
subset(fpkm_t, rnaseq_profile_id %in% subset(pheno, 
	pheno$genotype_driver == 'Snap25-IRES2-Cre')$rnaseq_profile_id) -> fpkm_Snap25
subset(fpkm_t, rnaseq_profile_id %in% subset(pheno, 
	pheno$genotype_driver == 'Slc32a1-IRES-Cre')$rnaseq_profile_id) -> fpkm_Slc32
subset(fpkm_t, rnaseq_profile_id %in% subset(pheno, 
	pheno$genotype_driver == 'Slc17a6-IRES-Cre')$rnaseq_profile_id) -> fpkm_Slc17
subset(fpkm_t, rnaseq_profile_id %in% subset(pheno, 
	pheno$genotype_driver == 'Gad2-IRES-Cre')$rnaseq_profile_id) -> fpkm_Gad2


## run pca for Snap25 Cre line
print("Computing principal components on Snap25")
fpkm_pr<-prcomp(fpkm_Snap25[,-ncol(fpkm_Snap25)],center=T,scale.=F)
unclass(fpkm_pr$rotation) -> fpkm_pr_rot
fpkm_pr_rot[1:5,1:5]
write.table(fpkm_pr_rot,file="results/batch_effect/fpkm_pr_Snap25.txt",quote=FALSE,
        sep="\t", row.names=TRUE, col.names=TRUE)


## run pca for Slc32 Cre line
print("Computing principal components on Slc32")
fpkm_pr<-prcomp(fpkm_Slc32[,-ncol(fpkm_Slc32)],center=T,scale.=F)
unclass(fpkm_pr$rotation) -> fpkm_pr_rot
fpkm_pr_rot[1:5,1:5]
write.table(fpkm_pr_rot,file="results/batch_effect/fpkm_pr_Slc32.txt",quote=FALSE,
        sep="\t", row.names=TRUE, col.names=TRUE)


## run pca for Slc17 Cre line
print("Computing principal components on Slc17")
fpkm_pr<-prcomp(fpkm_Slc17[,-ncol(fpkm_Slc17)],center=T,scale.=F)
unclass(fpkm_pr$rotation) -> fpkm_pr_rot
fpkm_pr_rot[1:5,1:5]
write.table(fpkm_pr_rot,file="results/batch_effect/fpkm_pr_Slc17.txt",quote=FALSE,
        sep="\t", row.names=TRUE, col.names=TRUE)


## run pca for Gad2 Cre line
print("Computing principal components on Gad2")
fpkm_pr<-prcomp(fpkm_Gad2[,-ncol(fpkm_Gad2)],center=T,scale.=F)
unclass(fpkm_pr$rotation) -> fpkm_pr_rot
fpkm_pr_rot[1:5,1:5]
write.table(fpkm_pr_rot,file="results/batch_effect/fpkm_pr_Gad2.txt",quote=FALSE,
        sep="\t", row.names=TRUE, col.names=TRUE)
