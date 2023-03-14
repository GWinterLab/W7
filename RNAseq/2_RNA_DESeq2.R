##############################
# DESeq2 analysis
##############################

rm(list=ls())

library("DESeq2")
library("edgeR")
require(ggplot2)
library("EDASeq")

list.files()

coldata=as.data.frame(read.table("coldata_K.txt",sep="\t", header=TRUE, row.names=1, dec="."))

infile=as.data.frame(read.table("count_table_K.txt",sep="\t", header=TRUE, row.names=1, dec=".", fill = TRUE))



head(infile)
dim(infile)


COUNTS <- infile
head(COUNTS)

samples <- colnames(COUNTS)
genes <- rownames(COUNTS)
COUNTS <- apply(as.matrix(COUNTS),2,as.numeric)
rownames(COUNTS) <- genes
head(COUNTS)

Data <- newSeqExpressionSet(counts=COUNTS)



## Filtering absent genes

Filter=0

detectionLimit <- Filter
filter <- rowSums(cpm(exprs(Data))> detectionLimit) >= 1 # cpm = function from edgeR and means 'counts-per-million'
nGenes <- length(filter)
nAbove <- sum(filter)
PpData <- newSeqExpressionSet(counts=exprs(Data)[filter,],featureData=fData(Data)[filter,],phenoData=pData(Data))
head(exprs(PpData))
dim(PpData)

data<-exprs(PpData)

# load the data with raw counts, where genes with low counts are filtered out (genes in rows, samples in colums)
infile<-data
dim(infile)

head(infile)







##################
# Run the DEseq
##################
dds <- DESeqDataSetFromMatrix(countData = infile, colData = coldata, design = ~ condition)

dds <- DESeq(dds)





##################################################################################################
# Data normalization: Variance stabilizing Transformation
##################################################################################################
reg_data=paste("star_RNA_noSpikeIn_Filter_",Filter,"_",sep="")

vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=TRUE) #default blind=TRUE
vstMat<-assay(vsd)
dim(vstMat)



write.table(vstMat, file=paste(reg_data,".DESeq2normalized-vst-model.blindT.txt",sep=""),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)



###############
# PCA plot
###############
k=2

pdf(paste(reg_data,"_DESeq2_vstMat_BLIND_fitType_parametric_EDAseq_plotPCA_k_",k,".pdf",sep=""))

EDASeq::plotPCA(vstMat,k=k,labels=T,isLog=T,col=colors_selected)
dev.off()






####################################
# DE analysis
####################################

tissue1="W7"

tissue2="DMSO"

print(paste(tissue1,tissue2,sep="_vs_"))



res <- DESeq2::results(dds, independentFiltering=T, cooksCutoff=T,contrast=c("condition",tissue1,tissue2))

head(res)
A<-as.data.frame(res)

write.table(A, file=paste(reg_data,"DESeq2_RNAseq_CooksT_FiltT.",paste(tissue1,tissue2,sep="_vs_"),".Filt_",Filter,".txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)

######
# meta-scores based on baseMean, log2FoldChange and pvalue

DN<-cbind(A[which(A[,2]<0),],log10(A[which(A[,2]<0),6]))
UP<-cbind(A[which(A[,2]>=0),],-log10(A[which(A[,2]>=0),6]))
colnames(DN)<-cbind(t(colnames(A)),"signLog10")
colnames(UP)<-cbind(t(colnames(A)),"signLog10")


order_log10adjP<-rbind(UP[order(UP[,7],decreasing=T),],DN[order(DN[,7],decreasing=T),])
head(order_log10adjP)

A_adapted <- cbind(order_log10adjP, order_log10adjP[,1], abs(order_log10adjP[,2]), -order_log10adjP[,5])

head(A_adapted)



meta_score<-apply(apply(A_adapted,2,rank,na.last=F, ties.method="max")[,c(8,9,10)],1,min)
head(meta_score)

A_final<-cbind(order_log10adjP,meta_score)
head(A_final)


DN_final<-cbind(A_final[which(A_final[,2]<0),],-(A_final[which(A_final[,2]<0),8]))
UP_final<-cbind(A_final[which(A_final[,2]>=0),],A_final[which(A_final[,2]>=0),8])
colnames(DN_final)<-cbind(t(colnames(A_final)),"signed_meta_score")
colnames(UP_final)<-cbind(t(colnames(A_final)),"signed_meta_score")

order_meta_score<-rbind(UP_final[order(UP_final[,9],decreasing=T),],DN_final[order(DN_final[,9],decreasing=T),])


write.table(order_meta_score, file=paste("DESeq2_RNAseq_CooksT_FiltT.",paste(tissue1,tissue2,sep="_vs_"),".meta_score.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)



####




####################
# DESeq2::plotMA(res)
####################

pdf(paste(reg_data,"_plotMA.",tissue1,"_vs_",tissue2,".padj05.Filt_",Filter,".pdf",sep=""))
DESeq2::plotMA(res,alpha = 0.05,main=c(paste(tissue1,tissue2,sep=" vs "),("(alpha <= 0.05)")))
dev.off()



####################
# Volcano plot
####################

library(extrafont)
loadfonts(device = "pdf")
fonts()

padj_threshold_plot=0.05
pdf(paste(reg_data,"_Volcano_RNA_samples_plot_log2FC1_vs_adjP.",tissue1,"_vs_",tissue2,".",padj_threshold_plot,".Filt_",Filter,".pdf",sep=""))

par(family = "Arial",cex=1.3)
plot(res$log2FoldChange,-log10(res$padj)
     ,pch=16
     ,main=c(paste(tissue1, "vs" , tissue2,"adj. P-value <= 0.05 & |log2FC| >= 1"))
     ,cex=0.35,xlab="Log2 Fold Change", ylab="-log10(adj. P-value)"
     ,family="Arial",col="#B4B4B4")
color_dn="#651DA1" 
color_up="#651DA1"
points(res[which((res$log2FoldChange)>=1 & res$padj<=0.05),]$log2FoldChange,-log10(res[which((res$log2FoldChange)>=1 & res$padj<=0.05),]$padj)
       ,pch=16,cex=0.35
      ,col=color_up)

points(res[which((res$log2FoldChange)<=-1 & res$padj<=0.05),]$log2FoldChange,-log10(res[which((res$log2FoldChange)<=-1 & res$padj<=0.05),]$padj)
       ,pch=16,cex=0.35
      ,col=color_dn)

dev.off()

}

