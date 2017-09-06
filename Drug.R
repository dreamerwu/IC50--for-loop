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
data=read.delim("D:/demo/xu_lab/drug/IC50_drc.txt",head=T,sep="\t")

#Test the best model (fctlist:LL.2,LL.3,LL.3u,LL.4,LL.5,W1.2,W1.3,W1.4,W2.2,W2.3,W2.4,BC.4,BC.5,LL2.2,LL2.3,LL2.3u,LL2.4,LL2.5,AR.2,AR.3,MM.2,MM.3)
Concentration=data[,1:1]    #defign the first column as Concentration
Response=data[,83:83]         #defign the Response group
df=data.frame(Concentration=Concentration,Response=Response)  #merge Concentration and Response into one group
fit=drm(Response~Concentration,data=df,fct=LL.3(),type="continuous",separate=TRUE)            #Built model
result=ED(fit,50)                                             #calculate IC50                                              #Extract IC50
plot(fit)
summary(fit)
getMeanFunctions
#conclusion: LL.3 is the best model, 297 cell lins are successfully calculated IC50, but 85 cell lines failed
#LL.3 is Log-logistic(ED50 as parameter) with lower limit at 0 (3 parms)


#Calculate IC50 at the same time using the best model
Concentration=data[,1:1]    #defign the first column as Concentration
n=ncol(data)                #calculate column number
output=matrix(data=NA,5,383)  #defign a matrix
for(i in 2:n) {             #for loop
Response=data[,i:i]         #defign the Response group
df=data.frame(Concentration=Concentration,Response=Response)  #merge Concentration and Response into one group
fit=drm(Response~Concentration,data=df,fct=LL.3(),type="continuous",separate=TRUE)            #Built model
result=ED(fit,50)                                             #calculate IC50                                              #Extract IC50
output[2:2,(i-1):(i-1)]=result[1:1,1:1]    #Put the IC50 value into matrix

summary=summary(fit)     
coefficient=summary$coefficients
pvalue=coefficient[1:3,4:4]   #extract Pvalue
output[3:5,(i-1):(i-1)]=pvalue  #put the p value into matrix

}
colname=1:383
output[1:1,]=colname
toutput=t(output)
write.csv(toutput,"D:/demo/xu_lab/drug/output/result_LL3.csv")       #output   


#to plot IC50 curve with multiple cell lines (success vs failed)
#success cell lines
success=read.delim("D:/demo/xu_lab/drug/success_IC50_drc.txt",head=T,sep="\t")
Concentration=success[,1:1]
Response=success[,2:2]       
df=data.frame(Concentration=Concentration,Response=Response)
fit=drm(Response~Concentration,data=df,fct=LL.3(),type="continuous",separate=TRUE) 
plot(fit,ylim=c(0,1.5),xlim=c(0,30),col="red")
n=ncol(success)
for(i in 3:n) {             
  Response=success[,i:i]       
  df=data.frame(Concentration=Concentration,Response=Response)
  fit=drm(Response~Concentration,data=df,fct=LL.3(),type="continuous",separate=TRUE) 
  plot(fit,add=TRUE,ylim=c(0,1.5),xlim=c(0,30),col="red")
}
#failed cell lines
fail=read.delim("D:/demo/xu_lab/drug/fail_IC50_drc.txt",head=T,sep="\t")
Concentration=fail[,1:1]
Response=fail[,2:2]       
df=data.frame(Concentration=Concentration,Response=Response)
fit=drm(Response~Concentration,data=df,fct=LL.3(),type="continuous",separate=TRUE) 
plot(fit,ylim=c(0,1.5),xlim=c(0,30),col="blue")
n=ncol(fail)
for(i in 3:n) {             
  Response=fail[,i:i]       
  df=data.frame(Concentration=Concentration,Response=Response)
  fit=drm(Response~Concentration,data=df,fct=LL.3(),type="continuous",separate=TRUE) 
  plot(fit,add=TRUE,ylim=c(0,1.5),xlim=c(0,30),col="blue")
  
  }



#Screen cell lines with P value of three predictor less than 0.05 & evaluate the model
m=ncol(output)
new_matrix=matrix(data=NA,6,383)
t=1
for (i in 1:m) {
  p1=output[3:3,i:i]
  p2=output[4:4,i:i]
  p3=output[5:5,i:i]
  p4=output[6:6,i:i]
  
        if (p1 < 0.05) {
          if (p2 <0.05) {
             if (p3 <0.05) {
               if (p4 <0.05) {
              new_matrix[1:6,t:t]=output[1:6,i:i]
              t=t+1
              }
            }
          }
        }
}

write.csv(new_matrix,"D:/demo/good_cell_lines.csv")









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


