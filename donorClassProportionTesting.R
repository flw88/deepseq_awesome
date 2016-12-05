#!/usr/bin/Rscript


## setwd and load files
print("Loading files...")


setwd("/Users/eliseflynn/Google Drive/Deep Seq/DeepSeq_Project")
pheno<-read.csv("Data/columns-cells.csv")

data.frame(table(pheno[c('donor_id','broad_class')])) -> broad_long
	reshape(broad_long, 
		idvar='donor_id', timevar='broad_class', 
		direction='wide') -> broad_wide 
	broad_wide[,1] -> row.names(broad_wide)
	broad_wide[,-1] -> broad_wide
	#try(fisher.test(broad_wide) -> fisher; print(fisher))
	chisq.test(broad_wide) -> chi
	print(chi)


for (cre in unique(pheno$genotype_driver)) {
	print(cre)
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

	print("Sub class")
	data.frame(table(cre_pheno[c('donor_id','subclass')])) -> sub_long
	reshape(sub_long, 
		idvar='donor_id', timevar='subclass', 
		direction='wide') -> sub_wide 
	sub_wide[,1] -> row.names(sub_wide)
	sub_wide[,-1] -> sub_wide
	#try(fisher.test(broad_wide) -> fisher; print(fisher))
	chisq.test(sub_wide) -> chi
	print(chi)

	
}
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
