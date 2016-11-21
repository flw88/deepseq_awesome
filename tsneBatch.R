#!/usr/bin/Rscript

library(ggplot2)
library(tsne)


## setwd and load files
print("Loading files...")
setwd("/home/local/users/eflynn/deepseq/")

args=commandArgs(tralingOnly=TRUE)
if (length(args)!=1) {
	stop("This script requires one argument (the pca input file)")
} else {
	input_pca = read.table(args[1],header=TRUE)
}

#fpkm<-read.csv("datasets/scRNASeq/fpkm_table.csv")
pheno<-read.csv("datasets/scRNASeq/columns-cells.csv")
#genes<-read.csv("datasets/scRNASeq/rows-genes.csv")


## label rnaseq details
row.names(input_pca) -> input_pca$rnaseq_profile_id
merge(input_pca,pheno,by="rnaseq_profile_id") -> pca_labeled


## plot scatter plots of PCA
ggplot(pca_labeled,aes(PC1,PC2,color=as.character(donor_id))) + geom_point()
ggsave(paste("results/batch_effect/PDFs/",gsub('.txt','',args[1]),"_PC1PC2.pdf",sep=''),
	plot=last_plot(), width=20,height=10)
ggplot(pca_labeled,aes(PC2,PC3,color=as.character(donor_id))) + geom_point()
ggsave(paste("results/batch_effect/PDFs/",gsub('.txt','',args[1]),"_PC2PC3.pdf",sep=''),
	plot=last_plot(), width=20,height=10)


## run tSNE
colors=rainbow(length(unique(pca_labeled$donor_id)))
names(colors) = unique(pca_labeled$donor_id)
ecb = function(x,y) { plot(x,t='n'); 
	text(x,labels=pca_labeled$broad_class
		col=colors[as.character(pca_labeled$donor_id)])}

for ( PCnum in c(5,10,20,30) ) {
print(paste("Running at PC level ",PCnum,sep=''))
	for ( perplexnum in c(10,20,50) ){
	print(paste("Running at perplexity ",perplexnum,sep=''))
		for (i in 1:10) {
			pdf(paste("results/batch_effect/PDFs/",
				gsub('.txt','',args[1]),
				"_tSNE_PC",PCnum,
				"_perplex",perplexnum,
				'.PDF',sep=''))
			tsne_run = tsne(pca_labeled[2:2+PCnum],
				epoch_callback = ecb
				perplexity = perplexnum)
			dev.off()
		}

	}

}

