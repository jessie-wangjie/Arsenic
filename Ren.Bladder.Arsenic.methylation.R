library(minfi)
library(reshape2)

## parse annotation
load("~/data/ren/methylation/analysis/GRSet.funnorm.RDat")
ann<-getAnnotation(GRSet.funnorm)
gene_name<-strsplit(ann$UCSC_RefGene_Name,";")
gene_group<-strsplit(ann$UCSC_RefGene_Group,";")
ann.450k<-data.frame(Name=rep(ann$Name,sapply(gene_name,length)),UCSC_RefGene_Name=unlist(gene_name),UCSC_RefGene_Group=unlist(gene_group))
ann.450k<-merge(ann.450k,ann[,c("chr","pos","strand","Islands_Name","Relation_to_Island","Regulatory_Feature_Group")],by.x=1,by.y=0,all=T)

## Beta
GRSet.funnorm<-GRSet.funnorm[,grep("JQ",pData(GRSet.funnorm)$Sample_Name)]
m_beta<-getBeta(GRSet.funnorm)
m_beta<-merge(ann[,c("chr","pos","strand","Islands_Name","Relation_to_Island","Regulatory_Feature_Group")],m_beta,by.x=0,by.y=0)

## MDS plot
pdf("Ren.Bladder.Arsenic.mds.pdf")
mdsPlot(as.matrix(m_beta[,-c(1:7)]),sampNames=pData(GRSet.funnorm)$Sample_Name)
dev.off();

## dmp, no replicate failed
# m_beta<-getBeta(GRSet.funnorm,betaThreshold=0.001)
# dmp<-dmpFinder(m_beta,pheno=pData(GRSet.funnorm)$Family,type="categorical")

## dmr
# dmr<-bumphunter(GRSet.funnorm[,1:2],design=model.matrix(~pData(GRSet.funnorm)$Family[1:2]),cutoff=0.2,B=0,type="Beta")
shortcpg<-cpgCollapse(GRSet.funnorm,what="Beta",returnBlockInfo=FALSE)
block<-blockFinder(shortcpg,design=model.matrix(~pData(GRSet.funnorm)$Family),what="Beta",cutoff=0.0001,smooth=FALSE)


## Sample Distance
sampleDists<-dist(t(getBeta(GRSet.funnorm)))  
library(gplots)
library('RColorBrewer')
sampleDistMatrix<-as.matrix(sampleDists)
colnames(sampleDistMatrix)<-pData(GRSet.funnorm)[colnames(sampleDistMatrix),"Sample_Name"]
rownames(sampleDistMatrix)<-pData(GRSet.funnorm)[rownames(sampleDistMatrix),"Sample_Name"]
colors<-colorRampPalette(rev(brewer.pal(9,'Blues')))(255)
hc<-hclust(sampleDists)
pdf("sampleDists.pdf")
heatmap.2(sampleDistMatrix,Rowv=as.dendrogram(hc),symm=TRUE,trace='none',col=colors,margins=c(10,10),cexRow=0.8,cexCol=0.8)
#heatmap.2(getBeta(GRSet.funnorm),distfun=function(c) as.dist(1-c),symm=TRUE,trace='none',col=colors,margins=c(10,10),cexRow=0.8,cexCol=0.8)
dev.off();


pdf("hist.pdf",width=12)
par(mfrow=c(2,4))
hist(getBeta(GRSet.funnorm)[,1],col="grey50", border="white",xlab="Beta",main=pData(GRSet.funnorm)[1,"Sample_Name"],ylim=c(0,100000))
hist(getBeta(GRSet.funnorm)[,2],col="grey50", border="white",xlab="Beta",main=pData(GRSet.funnorm)[2,"Sample_Name"],ylim=c(0,100000))
hist(getBeta(GRSet.funnorm)[,3],col="grey50", border="white",xlab="Beta",main=pData(GRSet.funnorm)[3,"Sample_Name"],ylim=c(0,100000))
hist(getBeta(GRSet.funnorm)[,4],col="grey50", border="white",xlab="Beta",main=pData(GRSet.funnorm)[4,"Sample_Name"],ylim=c(0,100000))
hist(getBeta(GRSet.funnorm)[,5],col="grey50", border="white",xlab="Beta",main=pData(GRSet.funnorm)[5,"Sample_Name"],ylim=c(0,100000))
hist(getBeta(GRSet.funnorm)[,6],col="grey50", border="white",xlab="Beta",main=pData(GRSet.funnorm)[6,"Sample_Name"],ylim=c(0,100000))
hist(getBeta(GRSet.funnorm)[,7],col="grey50", border="white",xlab="Beta",main=pData(GRSet.funnorm)[7,"Sample_Name"],ylim=c(0,100000))
hist(getBeta(GRSet.funnorm)[,8],col="grey50", border="white",xlab="Beta",main=pData(GRSet.funnorm)[8,"Sample_Name"],ylim=c(0,100000))
dev.off();
