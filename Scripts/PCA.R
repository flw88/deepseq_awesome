#PURPOSE: Do basic PCA as first pass to seeing data structure

#####setwd and load packages/libraries#####
library(ggfortify)
setwd("~/Google\ Drive/DeepSeq_Project")
#####load files#####
fpkm<-read.csv("fpkm_table.csv")
pheno<-read.csv("columns-cells.csv")
genes<-read.csv("rows-genes.csv")
#get attributes and re-format df
tmp<-colnames(fpkm)[-c(1)]
cells<-gsub("X","",tmp)
cells_class<-pheno[which(pheno$rnaseq_profile_id %in% cells),"broad_class"]
rownames(fpkm)<-fpkm$gene_id...rnaseq_profile_id
#realy crude but only want genes with fpkm>0
r<-apply(fpkm[,-c(1)],1,function(x){all(x!=0)})
sub_fpkm<-t(fpkm[,-c(1)])
rownames(sub_fpkm)<-cells
#compute components and plot
pr<-prcomp(sub_fpkm,center=T,scale.=F)
plot(pr, type = "l")
summary(pr)$importance[,1:2]
pdf("~/Google\ Drive/DeepSeq_Project/first_pass_pca.pdf",width=7,height=4)
autoplot(pr,data=pheno,shape="broad_class",
         colour="cell_prep_sample_id",size=1)+
  ggtitle("Spatial map of LGd neurons-first pass")+
  theme(
    title=element_text(face=2,size=8),
    axis.text=element_text(face=2),
    axis.title=element_text(face=2,size=8),
    legend.text=element_text(face=2,size=7)
  )
dev.off()

install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
pdf("~/Google\ Drive/DeepSeq_Project/first_pass_biplot_pcs23.pdf",width=7,height=4)
ggbiplot(pr,choices=2:3,
         var.axes=F,groups =pheno$broad_class,
         size=1)+
  ggtitle("Biplot PCs 1 and 2 of scLGd neurons")+
  theme(
    title=element_text(face=2,size=8),
    axis.text=element_text(face=2),
    axis.title=element_text(face=2,size=8),
    legend.text=element_text(face=2,size=7)
  )
dev.off()
