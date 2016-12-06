#!/usr/bin/Rscript

library(ggplot2)

## setwd and load files
print("Loading files...")


#setwd("/Users/eliseflynn/Google Drive/Deep Seq/DeepSeq_Project")
setwd("/cygdrive/c/Users/eflyn/Google Drive/Deep Seq/DeepSeq_Project/")
pheno<-read.csv("Data/columns-cells.csv")

broad_full <- pheno[c('donor_id','broad_class')]
data.frame(table(broad_full)) -> broad_long
reshape(broad_long, 
	idvar='donor_id', timevar='broad_class', 
	direction='wide') -> broad_wide 
broad_wide[,1] -> row.names(broad_wide)
broad_wide[,-1] -> broad_wide


fisher.test(broad_wide, simulate.p.value=T) -> fisher
print(fisher)
chisq.test(broad_wide) -> chi
print(chi)


pdf("Results/batch_effect/PDFs/all_donorVbroad.pdf")
ggplot(broad_full,aes(as.character(donor_id),
	fill=broad_class)) + 
	geom_bar() +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png("Results/batch_effect/PDFs/all_donorVbroad.png")
ggplot(broad_full,aes(as.character(donor_id),
	fill=broad_class)) + 
	geom_bar() +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


broad_full$category <- "Other"
broad_full[row.names(subset(broad_full,
	grepl('Gad2',broad_class))),]$category <- 'GABAergic'
broad_full[row.names(subset(broad_full,
	grepl('Slc17a6',broad_class))),]$category <- 'Glutamatergic'

png("Results/batch_effect/PDFs/all_donorVcat.png")
ggplot(broad_full,aes(as.character(donor_id),
fill=category)) + 
	geom_bar() +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf("Results/batch_effect/PDFs/all_donorVcat.pdf")
ggplot(broad_full,aes(as.character(donor_id),
	fill=category)) + 
	geom_bar() +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()	



for (cre in unique(pheno$genotype_driver)) {
	print(cre)
	cre_pheno <- subset(pheno, genotype_driver == cre)

	print("Broad class")
	broad_full <- cre_pheno[c('donor_id','broad_class')]
	data.frame(table(broad_full)) -> broad_long
	reshape(broad_long, 
		idvar='donor_id', timevar='broad_class', 
		direction='wide') -> broad_wide 
	broad_wide[,1] -> row.names(broad_wide)
	broad_wide[,-1] -> broad_wide

	fisher.test(broad_wide, simulate.p.value=T) -> fisher
	print(fisher)
	chisq.test(broad_wide) -> chi
	print(chi)

	pdf(paste("Results/batch_effect/PDFs/all_donorVbroad_",
		cre,".pdf",sep=''))
	ggplot(broad_full,aes(as.character(donor_id),
		fill=broad_class)) + 
		geom_bar() +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
	dev.off()

	png(paste("Results/batch_effect/PDFs/all_donorVbroad_",
		cre,".png",sep=''))
	ggplot(broad_full,aes(as.character(donor_id),
		fill=broad_class)) + 
		geom_bar() +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
	dev.off()


	broad_full$category <- "Other"
	try(broad_full[row.names(subset(broad_full,
		grepl('Gad2',broad_class))),]$category <- 'GABAergic')
	try(broad_full[row.names(subset(broad_full,
		grepl('Slc17a6',broad_class))),]$category <- 'Glutamatergic')
	data.frame(table(broad_full[,c('donor_id','category')])) -> broad_long
	reshape(broad_long, 
		idvar='donor_id', timevar='category', 
		direction='wide') -> broad_wide 
	broad_wide[,1] -> row.names(broad_wide)
	broad_wide[,-1] -> broad_wide

	fisher.test(broad_wide, simulate.p.value=T) -> fisher
	print(fisher)
	chisq.test(broad_wide) -> chi
	print(chi)

	png(paste("Results/batch_effect/PDFs/all_donorVcat_",
		cre,".png",sep=''))
	ggplot(broad_full,aes(as.character(donor_id),
		fill=category)) + 
		geom_bar() +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
	dev.off()

	pdf(paste("Results/batch_effect/PDFs/all_donorVcat_",
		cre,".pdf",sep=''))
	ggplot(broad_full,aes(as.character(donor_id),
		fill=category)) + 
		geom_bar() +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
	dev.off()	

}
#	for ( i in 1:nrow(broad_wide) ) {
#		for ( j in i:nrow(broad_wide) ) {
#			if ( i != j ) {
#				print(row.names(broad_wide)[i])
#				print(row.names(broad_wide)[j])
#				rows = c(i,j)
#				fisher.test(broad_wide[rows,]) -> fisher
#				print(fisher)
#				chisq.test(broad_wide[rows,]) -> chi
#				print(chi)
#			}
#		}
#	}

#	print("Sub class")
#	data.frame(table(cre_pheno[c('donor_id','subclass')])) -> sub_long
#	reshape(sub_long, 
#		idvar='donor_id', timevar='subclass', 
#		direction='wide') -> sub_wide 
#	sub_wide[,1] -> row.names(sub_wide)
#	sub_wide[,-1] -> sub_wide
#	#try(fisher.test(broad_wide) -> fisher; print(fisher))
#	chisq.test(sub_wide) -> chi
#	print(chi)
#
#	ggplot()
#	
#}
#long_broad <- pheno[c('donor_id','broad_class')]




cre="Snap25-IRES2-Cre"
cre_pheno <- subset(pheno, genotype_driver == cre)
print("Broad class")
data.frame(table(cre_pheno[c('donor_id','broad_class')])) -> broad_long
reshape(broad_long, 
	idvar='donor_id', timevar='broad_class', 
	direction='wide') -> broad_wide 
broad_wide[,1] -> row.names(broad_wide)
broad_wide[,-1] -> broad_wide
#try(fisher.test(broad_wide) -> fisher; print(fisher))
chisq.test(broad_wide) -> chi
print(chi)


for ( i in 1:nrow(broad_wide) ) {
	for ( j in i:nrow(broad_wide) ) {
		if ( i != j ) {
			print(row.names(broad_wide)[i])
			print(row.names(broad_wide)[j])
			rows = c(i,j)
			fisher.test(broad_wide[rows,]) -> fisher
			print(fisher)
			chisq.test(broad_wide[rows,]) -> chi
			print(chi)
		}
	}
}

## NOTE:

#[1] "Snap25-IRES2-Cre"
#[1] "Broad class"
#
#        Fisher's Exact Test for Count Data with simulated p-value (based on
#        2000 replicates)
#
#data:  broad_wide
#p-value = 0.001499
#alternative hypothesis: two.sided
#
#
#        Pearson's Chi-squared test
#
#data:  broad_wide
#X-squared = 71.638, df = 35, p-value = 0.0002548


#[1] "484742351"
#[1] "487244822"
#
#        Fisher's Exact Test for Count Data
#
#data:  broad_wide[rows, ]
#p-value = 0.006241
#alternative hypothesis: two.sided
#
#
#        Pearson's Chi-squared test
#
#data:  broad_wide[rows, ]
#X-squared = NaN, df = 5, p-value = NA


#[1] "487244822"
#[1] "487244826"
#
#        Fisher's Exact Test for Count Data
#
#data:  broad_wide[rows, ]
#p-value = 0.001791
#alternative hypothesis: two.sided
#
#
#        Pearson's Chi-squared test
#
#data:  broad_wide[rows, ]
#X-squared = NaN, df = 5, p-value = NA


#[1] "487244822"
#[1] "488448320"
#
#        Fisher's Exact Test for Count Data
#
#data:  broad_wide[rows, ]
#p-value = 0.001619
#alternative hypothesis: two.sided
#
#
#        Pearson's Chi-squared test
#
#data:  broad_wide[rows, ]
#X-squared = NaN, df = 5, p-value = NA
