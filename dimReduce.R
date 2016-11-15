etwd and load files
setwd("/results/G4017_Users/edf2126/team_awesome/")
fpkm<-read.csv("datasets/scRNASeq/fpkm_table.csv")
pheno<-read.csv("datasets/scRNASeq/columns-cells.csv")
genes<-read.csv("datasets/scRNASeq/rows-genes.csv")


## reformat fpkm df
fpkm[,1] -> row.names(fpkm)                     ## rename rows
fpkm[,-1] -> fpkm                                       ## remove gene id column
gsub("^X","",names(fpkm)) -> names(fpkm)        ## remove leading X from cell IDs
#fpkm <- t(fpkm)                                        ## transform so rows are cells, columns are genes (nvm, do later)

## Compute PCs of all genes, write to file
print("Computing principal components on all genes")
fpkm_pr<-prcomp(t(fpkm),center=T,scale.=F)
unclass(fpkm_pr$rotation) -> fpkm_pr_rot
fpkm_pr_rot[1:5,1:5]
write.table(fpkm_pr_rot,file="results/batch_effect/fpkm_pr.txt",quote=FALSE,
        sep="\t", row.names=TRUE, col.names=TRUE)


## Isolate genes w/ fpkm > 0
row_sub = apply(fpkm, 1, function(row) all(row!=0))
fpkm[row_sub,] -> sub_fpkm


## Compute components and plot
print("Computing principal components on genes fpkm > 0")
sub_fpkm_pr<-prcomp(t(sub_fpkm),center=T,scale.=F)
unclass(sub_fpkm_pr$rotation) -> sub_fpkm_pr_rot
sub_fpkm_rot[1:5,1:5]
write.table(sub_fpkm_pr_rot,file="results/batch_effect/sub_fpkm_pr.txt",quote=FALSE,
        sep="\t", row.names=TRUE, col.names=TRUE)
