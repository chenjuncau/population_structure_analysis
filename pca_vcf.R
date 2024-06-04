# Jun Chen
# 3/16/2023

# vcftools --vcf 1.vcf --maf 0.05 --max-missing 0.9 --relatedness
vcftools --vcf MAF0.2_DP10_Corrected_RD_171479_SNPs_COB_Samples.recode.renamed.annotate_NAME_CORRECTED_2.vcf --relatedness


# vcftools --vcf MAF0.2_DP10_Corrected_COB_183401_SNPs_annotated_renamed_juncopy.vcf --relatedness



# pca graph
setwd("V:\\31\\Jun\\rapidgenome\\06092023\\")
# setwd("V:\\31\\Jun\\rapidgenome\\3152023\\")


#sampleName.code.txt
library(reshape2)
# checking
# setwd("//cvirdlux05/work/31/Jun/biointerval/bgalGal1/newdata2")
# setwd("C:\\jchen\\google_computer\\gg7")


a_vs_g <- read.table("out.relatedness",head=TRUE,stringsAsFactors=F)# by vcftools
Matrix.value <- dcast(a_vs_g[,1:3], INDV1 ~ INDV2)  # it is not symmetric matrix .
# write.table(Matrix.value, file = "G.csv", sep = ",", col.names = TRUE,row.names=FALSE)  # create the plot in excel
Matrix.value[is.na(Matrix.value)] <- 0
Matrix.value <- as.matrix(Matrix.value)
Matrix.value2 <- Matrix.value[,-1]
rownames(Matrix.value2) <- Matrix.value[,1]
Matrix.value3=t(Matrix.value2)
Matrix.value2[lower.tri(Matrix.value2, diag = FALSE)]=Matrix.value3[lower.tri(Matrix.value3, diag = FALSE)]
write.table(Matrix.value2, file = "G.csv", sep = ",", col.names = TRUE,row.names=TRUE)  # create the plot in excel


d <- dist(t(as.matrix(Matrix.value2)))  
mode(d) <- 'numeric' 
# pdf("clusterplot_SEQ_value.pdf",width=150,height=25)    #complete  50 25   100 25
# pdf("clusterplot_SEQ_value.pdf",width=1000,height=25)    #complete  50 25   100 25
# plot(hclust(d,"ave"),cex = 1)
# dev.off()
pdf("clusterplot_SEQ_value.pdf",width=400,height=10)    #complete  50 25   100 25
# plot(hclust(d,"ave"),cex = 0.10)
plot(hclust(d,"ave"),cex = 0.05)   # save by right click using Rstudio graph click. 
dev.off() 
 
# install.packages("pvclust")

# MySNP=Matrix.value2
mode(Matrix.value2) <- 'numeric' 
# 4 
library(pvclust)
result <- pvclust(t(Matrix.value2), method.dist="cor", method.hclust="average", nboot=1000)
plot(result)

pdf("clusterplot_SEQ_value_bootstrap.pdf",width=200,height=8)
plot(result,cex=0.8)
dev.off()
  
  

# PCA

cps.full <- cmdscale(d, eig = T, k = 10)
names(cps.full)
cps <- cps.full$points
plot(cps[,1], cps[,2], pch = 1)

# pca
library(ggplot2)
pca=as.data.frame(cps)
# cobb_lines=factor(substring(rownames(pca),1,6))
cobb_lines=factor(substring(rownames(pca),1,3))
qplot(x=V1, y=V2, data=pca, colour=cobb_lines,xlab="PC1",ylab="PC2") + theme(legend.position="none")
write.table(pca,file="pca_2020.csv",append=FALSE,quote=FALSE,sep=",",row.names=TRUE,col.names=NA)

qplot(x=V1, y=V2, data=pca, colour=cobb_lines,xlab="PC1",ylab="PC2") 

qplot(x=V2, y=V3, data=pca, colour=cobb_lines,xlab="PC2",ylab="PC3") 

qplot(x=V1, y=V3, data=pca, colour=cobb_lines,xlab="PC1",ylab="PC3") 


# level plot
require(lattice)
pdf("levelplot_G.pdf",width=150,height=150)    #complete  50 50 100 100 
levelplot(Matrix.value2,scales=list(x=list(rot=90)),, main="G Matrix",col.regions = terrain.colors(100))

dev.off()

# library(ggplot2)
pdf("heatmap_G.pdf",width=100,height=100)    #complete
heatmap <- qplot(x=Var1, y=Var2, data=melt(Matrix.value2), geom="tile",fill=value) 
heatmap
dev.off()




# L12_SH_MG180_180320974
# L12_GM_MG191_penCrop_positive_191614771



library(scatterplot3d)
# jpeg('pca-3D.png', quality = 100, bg = "white", res = 200, width = 12, height = 12)
scatterplot3d(pca$V1,pca$V2,pca$V3, color=as.numeric(as.factor(as.character(cobb_lines))),xlab = "PC1", ylab = "PC2", zlab = "PC3",main="PCA")
legend("topleft", inset=.05,   pch = 19,  yjust=0,bty="n", cex=0.8,     title="COBB_LINES",   legend = levels(as.factor(as.character(cobb_lines))), col=seq_along(levels(as.factor(as.character(cobb_lines)))))    
# dev.off()



# path="Z:\\pca\\nonGS"
# setwd(path)
setwd("V:\\rapidgenome\\06092023")

library(ggplot2)
pca <- read.table("pca_2020.csv",header=TRUE,sep=",")
#group_info <- read.table("group_code_2019.csv",header=FALSE,sep=",")    # =MID(B1,1,5)
group_info <- read.table("sampleName.code.txt",header=TRUE,sep="\t")    # =MID(B1,1,5)
# colnames(group_info)=c("Rname","Rname2","code")	
colnames(group_info)=c("Rname","blupid")
group_info$code=factor(substring(group_info$blupid,1,3))	
#colnames(group_info)=c("Rname","blupid","code")		
pca$Rname=pca$X
pca2=merge(pca,group_info,by="Rname")
Groups=factor(pca2$code)
cobb_lines=factor(pca2$code)
qplot(x=V1, y=V2, data=pca2, colour=Groups,xlab="PC1",ylab="PC2") 
ggsave(file = "pca2020.pdf",,width=10,height=10)

		
	
		
		
