 l *.zip | wc -l
 for z in *.zip; do unzip $z; done
 l *.txt | wc -l
 rm -f *.zip

  
 
path="/work/31/Jun/rapidgenome/06092023"
setwd(path)
library(ggplot2)
snp_60k <- read.table("/biotech/genomic_selection/MAP/PT_map_2855_RAPID_FS9602_galGal7.txt",header=FALSE,sep="\t")
out.file<-""
out.filename<-NULL
file.names <- dir(path, pattern ="*MatrixAB.txt")
for(i in 1:length(file.names)){
file <- read.table(file=file.names[i],head=TRUE,sep="\t",check.names=F,colClasses="character",skip=0,row.names=1,na.strings = c("--"))

num_4 <- match(as.character(unlist(snp_60k[,1])),rownames(file))
num_5 <- which(num_4>0)
length(num_5)
Mydata <- file[c(num_4),]
Mydata2 <- data.frame(colnames(Mydata),file.names[i])
if(!all(rownames(Mydata)==rownames(out.file))) {print("error");print(file.names[i])}
out.file <- cbind(out.file, Mydata)
out.filename <- rbind(out.filename, Mydata2)
print(i)
}

Mydata <- t(out.file[,-1])
		
dim(Mydata)
dim(out.filename)
genotype11<-"AA"
genotype12<-"AB"
genotype22<-"BB"
for(i in 1:dim(Mydata)[2])
{
Mydata[c(which(Mydata[,i]==as.character(genotype11))),i]="0"
Mydata[c(which(Mydata[,i]==as.character(genotype12))),i]="1"
Mydata[c(which(Mydata[,i]==as.character(genotype22))),i]="2"
}
for(i in 1:dim(Mydata)[2])
{
Mydata[,i] <- as.numeric(Mydata[,i])
}
Sys.time()
d <- dist(as.matrix(Mydata)) 
Sys.time()

mode(d) <- 'numeric'  
cps.full <- cmdscale(d, eig = T, k = 10)
names(cps.full)
cps <- cps.full$points
pca=as.data.frame(cps)

cobb_lines=factor(substring(rownames(pca),1,3))
qplot(x=V1, y=V2, data=pca, colour=cobb_lines,xlab="PC1",ylab="PC2") 
ggsave(file = "pca_2019_chipdata.pdf",,width=10,height=10)


write.table(pca,file="pca_2019.csv",append=FALSE,quote=FALSE,sep=",",row.names=TRUE,col.names=NA)
write.table(out.filename,file="group_code_2019.csv",append=FALSE,quote=FALSE,sep=",",row.names=FALSE,col.names=FALSE)


# 3D

MAF0.2_DP10_Corrected_RD_171479_SNPs_COB_Samples.recode.renamed.annotate.vcf

