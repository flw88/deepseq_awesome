#!/usr/bin/Rscript

library(ggplot2)
library(tsne)


## setwd and load files
print("Loading files...")
setwd("/cygdrive/c/Users/eflyn/Google Drive/Deep Seq/DeepSeq_Project/")
input_pca_all <- read.table("Results/batch_effect/input/fpkm_hetgenes_pr.txt",
	header=TRUE,sep='\t')
input_pca_Snap25 <- read.table("Results/batch_effect/input/fpkm_hetgenes_pr_Snap25.txt",
	header=TRUE,sep='\t')
pheno<-read.csv("Data/columns-cells.csv")


######################## All RNASeqs ########################
row.names(input_pca_all) -> input_pca_all$rnaseq_profile_id
merge(input_pca_all,pheno,by="rnaseq_profile_id") -> pca_labeled

print("Running tSNE")
colors_donor = rainbow(length(unique(pca_labeled$donor_id)))
names(colors_donor) = unique(pca_labeled$donor_id)

colors_broad = rainbow(length(unique(pca_labeled$broad_class)))
names(colors_broad) = unique(pca_labeled$broad_class)


for ( PCnum in c(5,30) ) {
	print(paste("Running at PC level ",PCnum,sep=''))
	last_col=2+PCnum

	for ( perplexnum in c(10,20,50) ) {
		print(paste("Running at perplexity ",perplexnum,sep=''))

		for (i in 1:5) {

			tsne_run = tsne(pca_labeled[,2:last_col],
				#epoch_callback = ecb,
				#epoch = 1000,
				perplexity = perplexnum,
				max_iter=1000)


			tsne_df <- data.frame(tsne_run)
			names(tsne_df) <- c('X1','X2')
			cbind(tsne_df,pca_labeled) -> labeled


			## labeled by donor ID
			pdf(paste("Results/batch_effect/PDFs/",
				"fpkm_hetgenes_pr_tSNE",
				"_PC",PCnum,
					"_perplex",perplexnum,
				"_size4",
					"_",i,
					'_FINAL.pdf',sep=''),height=10,width=10)
			print(ggplot(labeled,aes(X1,X2,
				color=as.character(donor_id))) +
				geom_point(size=4) + 
				theme_classic())
			dev.off()

			png(paste("Results/batch_effect/PDFs/",
				"fpkm_hetgenes_pr_tSNE",
				"_PC",PCnum,
				"_perplex",perplexnum,
				"_size4",
				"_",i,
				'_FINAL.png',sep=''),
				height=2*300,width=2*300)
			print(ggplot(labeled,aes(X1,X2,
				color=as.character(donor_id))) +
				geom_point(size=4) + 
				theme_classic())
			dev.off()


			## labeled by broad class
			pdf(paste("Results/batch_effect/PDFs/",
				"fpkm_hetgenes_pr_tSNE_broadclass",
				"_PC",PCnum,
				"_perplex",perplexnum,
				"_size4",
				"_",i,
				'_FINAL.pdf',sep=''),height=10,width=10)
			print(ggplot(labeled,aes(X1,X2,
				color=as.character(broad_class))) +
				geom_point(size=4) + 
				theme_classic())
			dev.off()

			png(paste("Results/batch_effect/PDFs/",
				"fpkm_hetgenes_pr_tSNE_broadclass",
				"_PC",PCnum,
				"_perplex",perplexnum,
				"_size4",
				"_",i,
				'_FINAL.png',sep=''),
				height=2*300,width=2*300)
			print(ggplot(labeled,aes(X1,X2,
				color=as.character(broad_class))) +
				geom_point(size=4) + 
				theme_classic())
			dev.off()
	
		}

	}

}




######################## Snap25 RNASeqs ########################
row.names(input_pca_Snap25) -> input_pca_Snap25$rnaseq_profile_id
merge(input_pca_Snap25,pheno,by="rnaseq_profile_id") -> pca_labeled

print("Running tSNE")
colors_donor = rainbow(length(unique(pca_labeled$donor_id)))
names(colors_donor) = unique(pca_labeled$donor_id)

colors_broad = rainbow(length(unique(pca_labeled$broad_class)))
names(colors_broad) = unique(pca_labeled$broad_class)


for ( PCnum in c(5,30) ) {
	print(paste("Running at PC level ",PCnum,sep=''))
	last_col=2+PCnum

	for ( perplexnum in c(10,20,50) ) {
		print(paste("Running at perplexity ",perplexnum,sep=''))

		for (i in 1:2) {

			tsne_run = tsne(pca_labeled[,2:last_col],
				#epoch_callback = ecb,
				#epoch = 1000,
				perplexity = perplexnum,
				max_iter = 1000)


			tsne_df <- data.frame(tsne_run)
			names(tsne_df) <- c('X1','X2')
			cbind(tsne_df,pca_labeled) -> labeled


			## labeled by donor ID
			pdf(paste("Results/batch_effect/PDFs/",
				"fpkm_hetgenes_Snap25_pr_tSNE",
				"_PC",PCnum,
					"_perplex",perplexnum,
				"_size4",
					"_",i,
					'_FINAL.pdf',sep=''),height=10,width=10)
			print(ggplot(labeled,aes(X1,X2,
				color=as.character(donor_id))) +
				geom_point(size=4) + 
				theme_classic())
			dev.off()

			png(paste("Results/batch_effect/PDFs/",
				"fpkm_hetgenes_Snap25_pr_tSNE",
				"_PC",PCnum,
				"_perplex",perplexnum,
				"_size4",
				"_",i,
				'_FINAL.png',sep=''),
				height=2*300,width=2*300)
			print(ggplot(labeled,aes(X1,X2,
				color=as.character(donor_id))) +
				geom_point(size=4) + 
				theme_classic())
			dev.off()


			## labeled by broad class
			pdf(paste("Results/batch_effect/PDFs/",
				"fpkm_hetgenes_Snap25_pr_tSNE_broadclass",
				"_PC",PCnum,
					"_perplex",perplexnum,
				"_size4",
					"_",i,
					'_FINAL.pdf',sep=''),height=10,width=10)
			print(ggplot(labeled,aes(X1,X2,
				color=as.character(broad_class))) +
				geom_point(size=4) + 
				theme_classic())
			dev.off()

			png(paste("Results/batch_effect/PDFs/",
				"fpkm_hetgenes_Snap25_pr_tSNE_broadclass",
				"_PC",PCnum,
				"_perplex",perplexnum,
				"_size4",
				"_",i,
				'_FINAL.png',sep=''),
				height=2*300,width=2*300)
			print(ggplot(labeled,aes(X1,X2,
				color=as.character(broad_class))) +
				geom_point(size=4) + 
				theme_classic())
			dev.off()
	
		}

	}

}
