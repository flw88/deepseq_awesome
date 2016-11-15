#!/usr/bin/Rscript



## setwd and load files
print("Loading files...")
setwd("/home/local/users/eflynn/deepseq/")
fpkm<-read.csv("datasets/scRNASeq/fpkm_table.csv")
pheno<-read.csv("datasets/scRNASeq/columns-cells.csv")
genes<-read.csv("datasets/scRNASeq/rows-genes.csv")


## reformat fpkm df
fpkm[,1] -> row.names(fpkm)                     ## rename rows
fpkm[,-1] -> fpkm                               ## remove gene id column
gsub("^X","",names(fpkm)) -> names(fpkm)        ## remove leading X from cell IDs
#fpkm <- t(fpkm)                                ## transform so rows are cells, columns are genes (nvm, do later)

