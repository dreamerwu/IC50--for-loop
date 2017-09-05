#1. calculate correlation between gene expression and IC50
#install packages
install.packages("Hmisc")

#load packages
library("pheatmap")
library("Hmisc")

#input data
data1=read.delim("D:/demo/xu_lab/drug/IC50_gene1.txt",head=T,sep="\t",row.names=1)
data2=read.delim("D:/demo/xu_lab/drug/IC50_gene2.txt",head=T,sep="\t",row.names=1)
data3=read.delim("D:/demo/xu_lab/drug/IC50_gene3.txt",head=T,sep="\t",row.names=1)

#calculate correlation and P value
cor1=rcorr(as.matrix(data1),type="pearson")
cor2=rcorr(as.matrix(data2),type="pearson")
cor3=rcorr(as.matrix(data3),type="pearson")

#output
corp1=cor1$P
corr1=cor1$r
corp1_p=corp1[,1:2]
corr1_r=corr1[,1:2]
write.csv(corp1_p,"D:/demo/corp1.csv")
write.csv(corr1_r,"D:/demo/corr1.csv")


corp2=cor2$P
corr2=cor2$r
corp2_p=corp2[,1:2]
corr2_r=corr2[,1:2]
write.csv(corp2_p,"D:/demo/corp2.csv")
write.csv(corr2_r,"D:/demo/corr2.csv")


corp3=cor3$P
corr3=cor3$r
corp3_p=corp3[,1:2]
corr3_r=corr3[,1:2]
write.csv(corp3_p,"D:/demo/corp3.csv")
write.csv(corr3_r,"D:/demo/corr3.csv")


#2. plot r and p
library("ggplot2")
library("ggthemes")
data=read.delim("D:/demo/xu_lab/drug/IC50_gene_r_p.txt",head=T,sep="\t")
p=ggplot(data,aes(x=Pearson_r,y=P_value))
pp=p+geom_point(size=2,aes(color=group,shape=group))+theme_bw()
pp+scale_color_manual(values=c("blue","purple","gray","red","seagreen2","lightcoral"))


#3. calculate IC50 at same time
library("drc")
data=read.delim("D:/demo/demo.txt",head=T,sep="\t")
Concentration=data[,1:1]    #defign the first column as Concentration
n=ncol(data)                #calculate column number
output=matrix(data=NA,5,383)  #defign a matrix


for(i in 2:n) {             #for loop
Response=data[,i:i]         #defign the Response group
df=data.frame(Concentration=Concentration,Response=Response)  #merge Concentration and Response into one group
fit=drm(Response~Concentration,data=df,fct=LL.4(),type="continuous",separate=TRUE)            #Built model
result=ED(fit,50)                                             #calculate IC50                                              #Extract IC50
output[1,(i-1):(i-1)]=result[1,1]    #Put the IC50 value into matrix

summary=summary(fit)     
coefficient=summary$coefficients
pvalue=coefficient[1:4,4:4]   #extract Pvalue
output[2:5,(i-1):(i-1)]=pvalue  #put the p value into matrix

}
write.csv(output,"D:/demo/result.csv")       #output    


#Screen cell lines with P value of three predictor less than 0.05





################################################################

plot(fit)

names(summary(fit))

summary=summary(fit)
names(summary)
coefficient=summary$coefficients
coefficient[1:4,4:4]


names(fit)

fit$indexMat2

print(fit)

































cor_matrix=cor(data)
output=cor_matrix[1:1,]
write.csv(output,"D:/demo/demo.csv")





data=read.delim("D:/demo/xu_lab/drug/IC50.txt",head=T,sep="\t")
ncol(data)
data2=data[,2:257]
pheatmap(data2,cluster_cols=FALSE,cluster_rows=FALSE,cellwidth=5,cellheight=10)




library("ggplot2")
library("ggthemes")
p=ggplot(data,aes(x=Order,y=IC50,color=Organ))
p+geom_point(size=5)+theme_bw()


IC50=data$IC50
hist(IC50,breaks=seq(0,800,5))


