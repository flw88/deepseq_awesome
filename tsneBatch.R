#!/usr/bin/Rscript

library(ggplot2)
library(tsne)


## setwd and load files
print("Loading files...")

args=commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
	stop("This script requires one argument (the pca input file)")
} else {
	input_pca = read.table(args[1],header=TRUE)
	input_pca_name = basename(args[1])
}


setwd("/home/local/users/eflynn/deepseq/")
#fpkm<-read.csv("datasets/scRNASeq/fpkm_table.csv")
pheno<-read.csv("datasets/scRNASeq/columns-cells.csv")
#genes<-read.csv("datasets/scRNASeq/rows-genes.csv")


## label rnaseq details
row.names(input_pca) -> input_pca$rnaseq_profile_id
merge(input_pca,pheno,by="rnaseq_profile_id") -> pca_labeled


## plot scatter plots of PCA
print("Plotting PCAs")
ggplot(pca_labeled,aes(PC1,PC2,color=as.character(donor_id))) + 
	geom_point() -> temp

png(paste("results/batch_effect/PDFs/",gsub('.txt','',input_pca_name),"_PC1PC2.png",sep=''),
	height=2*300,width=2.5*300)
temp
dev.off()

pdf(paste("results/batch_effect/PDFs/",gsub('.txt','',input_pca_name),"_PC1PC2.pdf",sep=''),
	height=10,width=10)
temp
dev.off()
#ggsave(paste("results/batch_effect/PDFs/",gsub('.txt','',input_pca_name),"_PC1PC2.pdf",sep=''),
#	plot=last_plot(), width=20,height=10)

ggplot(pca_labeled,aes(PC2,PC3,color=as.character(donor_id))) + 
	geom_point() -> temp

png(paste("results/batch_effect/PDFs/",gsub('.txt','',input_pca_name),"_PC2PC3.png",sep=''),
	height=2*300,width=2.5*300)
temp
dev.off()

pdf(paste("results/batch_effect/PDFs/",gsub('.txt','',input_pca_name),"_PC2PC3.pdf",sep=''),
	height=10,width=10)
temp
dev.off()

#ggsave(paste("results/batch_effect/PDFs/",gsub('.txt','',input_pca_name),"_PC2PC3.pdf",sep=''),
#	plot=last_plot(), width=20,height=10)


## run tSNE
print("Running tSNE")
colors=rainbow(length(unique(pca_labeled$donor_id)))
names(colors) = unique(pca_labeled$donor_id)
#ecb = function(x,y) { 
#	png(paste("results/batch_effect/PDFs/",
#		gsub('.txt','',input_pca_name),
#		"_tSNE_PC",PCnum,
#		"_perplex",perplexnum,
#		"_",i,
#		'.png',sep=''),
#		height=2*300,width=2.5*300)
#	plot(x,t='n');
#	text(x,labels=pca_labeled$broad_class,
#		col=colors[as.character(pca_labeled$donor_id)])
#	dev.off()
#}

for ( PCnum in c(5,10,20,30) ) {
	print(paste("Running at PC level ",PCnum,sep=''))
	last_col=2+PCnum

	for ( perplexnum in c(10,20,50) ) {
		print(paste("Running at perplexity ",perplexnum,sep=''))

		for (i in 1:5) {

			ecb = function(x,y) { 
				png(paste("results/batch_effect/PDFs/",
					gsub('.txt','',input_pca_name),
					"_tSNE_PC",PCnum,
					"_perplex",perplexnum,
					"_",i,
					'.png',sep=''),
					height=2*300,width=2*300)
				plot(x,t='n');
				text(x,labels=pca_labeled$broad_class,
					col=colors[as.character(pca_labeled$donor_id)])
				dev.off()
			
				pdf(paste("results/batch_effect/PDFs/",	
					gsub('.txt','',input_pca_name),
					"_tSNE_PC",PCnum,
					"_perplex",perplexnum,
					"_",i,
					'.pdf',sep=''),height=10,width=10)
				plot(x,t='n');
				text(x,labels=pca_labeled$broad_class,
					col=colors[as.character(pca_labeled$donor_id)])
				dev.off()
			}

			tsne_run = tsne(pca_labeled[2:last_col],
				epoch_callback = ecb,
				epoch = 1000,
				perplexity = perplexnum)

#			png(paste("results/batch_effect/PDFs/",gsub('.txt','',input_pca_name),"_PC2PC3_2.png",sep=''),
#				height=2*300,width=2.5*300)
#			temp
#			dev.off()

			#plot(tsne_run$Y,labels=pca_labeled$broad_class,
			#	col=colors[as.character(pca_labeled$donor_id)])

#			png(paste("results/batch_effect/PDFs/",
#				gsub('.txt','',input_pca_name),
#				"_tSNE_PC",PCnum,
#				"_perplex",perplexnum,
#				"_",i,
#				'.png',sep=''),height=2*300,width=2.5*300)
#			temp
#			dev.off()
	
		}

	}

}

