#PURPOSE: For doing exploratorey analyses of the scRNASeq

#####SETWD#####
setwd("~/Google\ Drive/DeepSeq_Project/Results")
#####load files#####
setwd("~/Google\ Drive/DeepSeq_Project/Data")
fpkm<-read.csv("fpkm_table.csv")
pheno<-read.csv("columns-cells.csv")
genes<-read.csv("rows-genes.csv")
gtf<-read.table("rsem_GRCm38.p3.gtf",sep="\t",stringsAsFactors = F)
counts<-read.csv("count_table.fpkm2count.csv")
#####load libraries#####
library(genefilter)
library(RColorBrewer)
library(gplots)
#####get attributes and re-format df#####
tmp<-colnames(fpkm)[-c(1)]
cells<-gsub("X","",tmp)
cells_class<-pheno[which(pheno$rnaseq_profile_id %in% cells),"broad_class"]
rownames(fpkm)<-fpkm$gene_id...rnaseq_profile_id
#####subsetting fpkm#####
#r<-apply(fpkm[,-c(1)],1,function(x){all(x!=0)})
sub_fpkm<-t(fpkm[,-c(1)])
rownames(sub_fpkm)<-cells
#####agg by pheno#####
agg_class<-aggregate(sub_fpkm,by=list(pheno$broad_class),FUN=mean)
agg_driver<-aggregate(sub_fpkm,by=list(pheno$genotype_driver),FUN=mean)
sub_agg<-agg_class[,-c(1)]
rownames(sub_agg)<-agg_class$Group.1
#####scRNA variabiliity ppt#####
rownames(counts)<-colnames(sub_fpkm)
mod_counts<-t(counts)[-c(1),]
tmp<-mod_counts
#a look at the cell's quantified RNA
quants<-apply(tmp,1,function(x){quantile(x,probs=seq(0.7,1,0.05))})
#for each quantile of expression across the cells,
#what is the max gene value?
max(unlist(lapply(1:nrow(quants),function(x){max(quants[x,])})))
#what is the min gene value?
lapply(1:nrow(quants),function(x){min(quants[x,])})
#what is the mean gene count across cells?
lapply(1:nrow(quants),function(x){mean(quants[x,])})
#for each quantile of expression across the cells,
#what is the percent drop out rate?
lapply(1:nrow(quants),function(x){length(which(quants[x,]==0))/length(quants[x,])})
#how do gene levels differ between cells?
pdf(paste0(getwd(),"/Example_gene_variability.pdf"),width=5,height=5)
smoothScatter(log10(tmp[1,]+1),log10(tmp[2,]+1), pch=1,xlab="cell #1 log10 FPKM",ylab="cell #2 log10 FPKM",main="Example gene levels between Gad2_Chrna6 cells")
dev.off()
#cell-to-cell variability expanded
n=5
bclasses<-unique(pheno$broad_class)
inds<-which(pheno$broad_class %in% bclasses[1])
raninds<-sample(inds,n)
pdf(paste0(getwd(),"/expanded_cell_variability_scatterplot.pdf"),width=20,height=20)
par(mfrow=c(n,n))
for(i in raninds){
  for(j in raninds){
    if(i==j){
      hist(log10(tmp[i,]+1),breaks=50,xlim=c(0,1),xlab="expression bin",main="")
    }else{
          smoothScatter(log10(tmp[i,]+1),log10(tmp[j,]+1), pch=1,xlab="",ylab="")
    }
  }
}
dev.off()
#CV v sample mean
cov_allcells<-lapply(1:nrow(tmp),function(x){sd(tmp[x,])/mean(tmp[x,])^2})
smean_allcells<-lapply(1:nrow(tmp),function(x){mean(tmp[x,])})
x<-unlist(smean_allcells)
y<-unlist(cov_allcells)
pdf(paste0(getwd(),"/CV2_vs_mean_all_cells.pdf"),width=5,height=5)
smoothScatter(log10(x),log(y),ylab="CV^2 (log10)",xlab="Mean (log10)",main="Density of gene variation across cells")
dev.off()
#example variability in select batches
tmp_df<-batch_dfs[[which.max(cors)]]
cov<-lapply(1:nrow(tmp_df),function(x){(sd(tmp_df[x,])/mean(tmp_df[x,]))^2})
smean<-lapply(1:nrow(tmp_df),function(x){mean(tmp_df[x,])})
x<-unlist(smean)
y<-unlist(cov)
cors<-c(cors,cor(x,y))
pdf(paste0(getwd(),"/high_cor_cov_vs_mean_gene_density_example.pdf",width=5,height=5))
plot(log10(x),log(y),pch=19,ylab="CV^2 (log10)",xlab="sample mean (log10)",main=paste0("Density of gene variation\nacross ",names(batch_dfs)[prep]," prep id; Cor ",round(cor(x,y),2)))
dev.off()
tmp_df<-batch_dfs[[which.min(cors)]]
cov<-lapply(1:nrow(tmp_df),function(x){(sd(tmp_df[x,])/mean(tmp_df[x,]))^2})
smean<-lapply(1:nrow(tmp_df),function(x){mean(tmp_df[x,])})
x<-unlist(smean)
y<-unlist(cov)
cors<-c(cors,cor(x,y))
pdf(paste0(getwd(),"/low_cor_cov_vs_mean_gene_density_example.pdf"),width=5,height=5)
plot(log10(x),log(y),pch=19,ylab="CV^2 (log10)",xlab="sample mean (log10)",main=paste0("Density of gene variation\nacross ",names(batch_dfs)[prep]," prep id; Cor ",round(cor(x,y),2)))
dev.off()
#####Figures#####
tmp<-mod_counts
pdf(paste0(getwd(),"/Figure_example_gene_variability.pdf"),width=6,height=5)
par(cex=1,font.main=2,font.lab=2)
plot(log10(tmp[1,]+1),log10(tmp[2,]+1), pch=1,xlab="cell #1 log10 counts",ylab="cell #2 log10 counts",main="Example gene variability between Gad2_Chrna6 cells")
dev.off()

tmp<-mod_counts
LCounts <- log10(t(tmp)+1)
Lmeans <- rowMeans( LCounts )
Lvars <- unlist(apply(LCounts,1,var))
Lcv2 <- Lvars / Lmeans^2
LogNcountsList = list()
#75 percentile of all cell means
cutoff<-quantile(Lmeans,prob=seq(0,1,0.05))["75%"]
useForFitL = Lmeans>cutoff #number of het genes tapers after .7
LogNcountsList$mean = Lmeans[useForFitL]
LogNcountsList$cv2 = Lcv2[useForFitL]
#predefined function plot(LogNcountsList$mean,10*10^(-2*LogNcountsList$mean))
#fitting function to cv2 starting from more or less arb. params
fit_loglin = nls(cv2 ~ a* 10^(-k*mean), LogNcountsList,start=c(a=10,k=2))
#defining estimated function params
a<-coefficients(fit_loglin)["a"]
k<-coefficients(fit_loglin)["k"]
#defining threshold and labeling het genes
#het genes have to have above the mean threshold and less than Lcv2
is_het = (a *10^(-k*Lmeans) < Lcv2) &  Lmeans>cutoff
#fitted line defined
LogVar_tech_logfit <- a *10^(-k*Lmeans)*Lmeans^2
#plot mean/cv2 relationship and variable genes
#making y axis log scale so we can see the variability better
pdf(paste0(getwd(),"/Figure_techvbiol_variable_genes.pdf"),height=5,width=5)
par(cex=1,font.main=2,font.lab=2)
plot( Lmeans, Lcv2, log="y", col=1+is_het,xlab='mean',ylab='CV2 (log10)',ylim=c(1e-3,1e2))  
xg <- seq( 0, 4.5, length.out=100 )
#line that defines "technical variance barrier" based on estimated params
lines( xg, a *10^(-k*xg ),lwd=2,col='blue' )
legend('top',c('Variable genes'),pch=c(1),col=c('red'),cex=0.9)
dev.off()
##########exploring##########
#####example scatterplots#####
sub_agg<-sapply(sub_agg,as.numeric)
sub_fpkm<-t(fpkm[,-c(1)])
cov<-lapply(1:nrow(sub_fpkm),function(x){(sd(sub_fpkm[x,])/mean(sub_fpkm[x,]))^2})
smean<-lapply(1:nrow(sub_fpkm),function(x){mean(sub_fpkm[x,])})
c<-unlist(cov)
s<-unlist(smean)
x<-s[order(s)]
y<-c[order(s)]
smoothScatter(log10(x),log(y),ylab="CV^2 (log10)",xlab="sample mean (log10)",main="Density of gene variation across cells")
smoothScatter(log10(agg[1,-c(1)]),log10(agg[2,-c(1)]))
smoothScatter(log10(sub_agg[1,]),log10(sub_agg[3,]), pch=1)
#####fig 1a like in kharenchenko et al.#####
smoothScatter(log10(sub_fpkm[1,]+1),log10(sub_fpkm[2,]+1), pch=1)
#do the above for the first 7 cells for each broad class 
#to see how much drop outs effects cell-to-cell variability
plot(density(log1p(sub_fpkm[1,])),xlim=c(-1,6))
#####get prep batches#####
tmp<-mod_counts
batches<-unique(pheno$cell_prep_sample_id)
batch_rows<-NULL
batch_dfs<-NULL
for(i in 1:length(batches)){
  rows<-rownames(subset(pheno, cell_prep_sample_id==batches[i]))
  batch_rows[[i]]<-as.numeric(rows)
  batch_dfs[[i]]<-tmp[as.numeric(rows),]
}
names(batch_rows)<-batches
names(batch_dfs)<-batches
cors<-c()
for(i in 1:length(batches)){
  prep<-i
  tmp_df<-batch_dfs[[prep]]
  cov<-lapply(1:nrow(tmp_df),function(x){(sd(tmp_df[x,])/mean(tmp_df[x,]))^2})
  smean<-lapply(1:nrow(tmp_df),function(x){mean(tmp_df[x,])})
  x<-unlist(smean)
  y<-unlist(cov)
  cors<-c(cors,cor(x,y))
  plot(log10(x),log(y),pch=19,ylab="CV^2 (log10)",xlab="sample mean (log10)",main=paste0("Density of gene variation\nacross ",names(batch_dfs)[prep]," prep id; Cor ",round(cor(x,y),2)))
}
#
#####get variability by classes#####
tmp<-mod_counts
classes<-unique(pheno$broad_class)
class_rows<-NULL
class_dfs<-NULL
for(i in 1:length(classes)){
  rows<-rownames(subset(pheno, broad_class==classes[i]))
  class_rows[[i]]<-as.numeric(rows)
  class_dfs[[i]]<-tmp[as.numeric(rows),]
}
names(class_rows)<-classes
names(class_dfs)<-classes
cors<-c()
for(i in 1:length(classes)){
  class<-i
  tmp_df<-class_dfs[[class]]
  cov<-lapply(1:nrow(tmp_df),function(x){(sd(tmp_df[x,])/mean(tmp_df[x,]))^2})
  smean<-lapply(1:nrow(tmp_df),function(x){mean(tmp_df[x,])})
  x<-unlist(smean)
  y<-unlist(cov)
  cors<-c(cors,cor(x,y))
  smoothScatter(log10(x),log(y),ylab="CV^2 (log10)",xlab="sample mean (log10)",main=paste0("Density of gene variation\nacross ",names(class_dfs)[class]," class; Cor ",round(cor(x,y),2)))
}
#####across batches technical vs biologically variable genes#####
#from https://github.com/PMBio/scLVM/blob/master/R/scripts/transform_counts_demo_no_spikeins.Rmd
##
##Notes##
##determining nonlinear (weighted) least-squares estimates of the parameters
#of a nonlinear model.
##
##
tmp<-mod_counts
LCounts <- log10(t(tmp)+1)
Lmeans <- rowMeans( LCounts )
Lvars <- unlist(apply(LCounts,1,var))
Lcv2 <- Lvars / Lmeans^2
LogNcountsList = list()
#75 percentile of all cell means
cutoff<-quantile(Lmeans,prob=seq(0,1,0.05))["75%"]
useForFitL = Lmeans>cutoff #number of het genes tapers after .7
LogNcountsList$mean = Lmeans[useForFitL]
LogNcountsList$cv2 = Lcv2[useForFitL]
#predefined function plot(LogNcountsList$mean,10*10^(-2*LogNcountsList$mean))
#fitting function to cv2 starting from more or less arb. params
fit_loglin = nls(cv2 ~ a* 10^(-k*mean), LogNcountsList,start=c(a=10,k=2))
#defining estimated function params
a<-coefficients(fit_loglin)["a"]
k<-coefficients(fit_loglin)["k"]
#defining threshold and labeling het genes
#het genes have to have above the mean threshold and less than Lcv2
is_het = (a *10^(-k*Lmeans) < Lcv2) &  Lmeans>cutoff
#fitted line defined
LogVar_tech_logfit <- a *10^(-k*Lmeans)*Lmeans^2
#plot mean/cv2 relationship and variable genes
#making y axis log scale so we can see the variability better
plot( Lmeans, Lcv2, log="y", col=1+is_het,xlab='meansLog',ylab='cv2Log',ylim=c(1e-3,1e2))  
xg <- seq( 0, 4.5, length.out=100 )
#line that defines "technical variance barrier" based on estimated params
lines( xg, a *10^(-k*xg ),lwd=2,col='blue' )
legend('top',c('Variable genes'),pch=c(1),col=c('red'),cex=0.8)
across_batches_het_genes<-names(which(is_het))
write.table(across_batches_het_genes,file=paste0(getwd(),"/across_batches_het_genes.txt"),quote=F,row.names=F,col.names=F)
#####batch-specific biologically variable genes#####
#from https://github.com/PMBio/scLVM/blob/master/R/scripts/transform_counts_demo_no_spikeins.Rmd
##
##Notes##
##determining nonlinear (weighted) least-squares estimates of the parameters
#of a nonlinear model.
##
##
het_genes<-NULL
for (i in 1:length(batch_dfs)){
  LCounts <- log10(t(batch_dfs[[i]])+1)
  Lmeans <- rowMeans( LCounts )
  Lvars <- unlist(apply(LCounts,1,var))
  Lcv2 <- Lvars / Lmeans^2
  LogNcountsList = list()
  #85 percentile of all cell means
  cutoff<-quantile(Lmeans,prob=seq(0,1,0.05))["75%"]
  useForFitL = Lmeans>cutoff #number of het genes tapers after .7
  LogNcountsList$mean = Lmeans[useForFitL]
  LogNcountsList$cv2 = Lcv2[useForFitL]
  #predefined function plot(LogNcountsList$mean,10*10^(-2*LogNcountsList$mean))
  #fitting function to cv2 starting from more or less arb. params
  fit_loglin = nls(cv2 ~ a* 10^(-k*mean), LogNcountsList,start=c(a=10,k=2))
  #defining estimated function params
  a<-coefficients(fit_loglin)["a"]
  k<-coefficients(fit_loglin)["k"]
  #defining threshold and labeling het genes
  #het genes have to have above the mean threshold and less than Lcv2
  is_het = (a *10^(-k*Lmeans) < Lcv2) &  Lmeans>cutoff
  #fitted line defined
  LogVar_tech_logfit <- a *10^(-k*Lmeans)*Lmeans^2
  # #plot mean/cv2 relationship and variable genes
  # #making y axis log scale so we can see the variability better
  # plot( Lmeans, Lcv2, log="y", col=1+is_het,xlab='meansLog',ylab='cv2Log',ylim=c(1e-3,1e2))  
  # xg <- seq( 0, 4.5, length.out=100 )
  # #line that defines "technical variance barrier" based on estimated params
  # lines( xg, a *10^(-k*xg ),lwd=2,col='blue' )
  # legend('top',c('Variable genes'),pch=c(1),col=c('red'),cex=0.8)
  het_genes[[names(batch_dfs)[i]]]<-names(which(is_het))
}
union_of_batches_het_genes<-unique(unlist(het_genes))
intersectList<-function(l=list,len=length(list)){
  intersect(l[[1]],
            intersectListAgain(l,len)<-function(l,x){
               intersect(l[[2]],l[[len]])
            }
              )
}
intersection_of_batches_het_genes<-unique(unlist(het_genes))
write.table(unique(unlist(union_of_batches_het_genes)),file=paste0(getwd(),"/union_of_batches_het_genes.txt"),quote=F,row.names=F,col.names=F)
#####hierarchical clustering#####
tmp<-mod_counts
het_genes<-as.character(read.table(paste0(getwd(),"/Results/het_genes.txt"),stringsAsFactors=F)[,1])
m<-tmp[,het_genes]
dist<-dist(m)
hc<-hclust(dist,method="ward.D")
p<-pheno$broad_class
cols=rainbow(length(unique(p)))
cut<-length(unique(pheno$broad_class))
cutTier3<-cutree(hc,cut)
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- cols[cutTier3[which(names(cutTier3) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}
clusDendro = dendrapply(as.dendrogram(hc), colLab)
pdf("~/Google\ Drive/DeepSeq_Project/Results/broad_class_dendrogram.pdf",width=50,height=15)
plot(clusDendro,main="Ward Dendrogram of LGd")
legend(-80,11e6,legend=levels(p),xpd=T,text.width=1.8,inset=.1,horiz=F,col=cols,pch=15,cex=1,bty="n")
dev.off()
#####heatmap of top varying genes#####
tmp<-t(sub_fpkm)
#bvargenes<-read.table(paste0(getwd(),"/union_of_batches_het_genes.txt"),stringsAsFactors = F)[,1]
topVarGenes <- head(order(-rowVars(tmp)),500)
#topVarGenes<-rownames(tmp) %in% bvargenes
hm<-tmp[topVarGenes,]
hm<-hm-colMeans(hm)
class<-as.factor(pheno$broad_class)
dist_hm<-dist(hm)
hc_hm<-hclust(dist_hm,method="ward.D")
p_hm<-class[hc_hm$order]
cols=rainbow(length(unique(p_hm)))
sidecols <- cols[class]
pdf("~/Google\ Drive/DeepSeq_Project/Results/heatmap_top500VaryingGenes.pdf",width=20,height=15)
par(cex.main=1)
heatmap.2(hm, trace="none",
          Rowv=F,ColSideColors=sidecols,
          labRow=rownames(hm),key=T,
          mar=c(10,2), scale="row",labCol=F,
          main="LGd scRNASeq hierarchical clustering with top 500 varying genes",
          cexRow=1,margins=c(4,27),
          dendrogram="column")
legend(0,15000,legend=levels(p_hm),xpd=T,text.width=1.8,inset=.1,horiz=F,col=cols,pch=15,cex=1,bty="n")
dev.off()
#####save image#####
setwd("~/Google\ Drive/DeepSeq_Project/RDa")
save.image(paste0(getwd(),"/exploratoryAnalyses.RDa"))
