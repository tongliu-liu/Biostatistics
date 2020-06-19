library(glmnet)
setwd ( "E:/生物统计/结课论文数据/结课论文数据/Regression" )

###Load the data of input and output
x<-as.matrix(read.table(file="1130sgRNA_100bp_7mer.txt"))
sgRNA.1130<-read.csv(file="sgRNA.1130.csv",header=T)
sgRNA.1130.names<-as.character(sgRNA.1130[,1])
sgRNA.top20.id<-sgRNA.1130.names[1:20]
rownames(x)<-sgRNA.1130.names
y<-sgRNA.1130[,5]
names(y)<-sgRNA.1130.names
###Load the data of input and output

##Construct the training groups and testing groups
train.groups.id<-list()
test.groups.id<-list()
x.train<-list()
x.test<-list()
y.train<-list()
y.test<-list()
y.test.total<-c()
tmp<-c()
for(i in 1:10){
  test.groups.id[[i]]<-sgRNA.1130.names[which(sgRNA.1130[1:980,6]==i)]
  train.groups.id[[i]]<-sgRNA.1130.names[which(sgRNA.1130[,6]!=i)]
  x.train[[i]]<-x[train.groups.id[[i]],]
  x.test[[i]]<-x[test.groups.id[[i]],]
  y.train[[i]]<-y[train.groups.id[[i]]]
  y.test[[i]]<-y[test.groups.id[[i]]]
  tmp<-c(tmp,test.groups.id[[i]])
  y.test.total<-c(y.test.total,y.test[[i]])
}
names(y.test.total)<-tmp
##Construct the training groups and testing groups

##Perform Lasso
library(glmnet)
sum_y_top20<-rep(0,20)
sum_y_total<-rep(0,980)
for(m in 1:10){
  cat("m=",m,'\t')
  y.test.pred.total<-c()
  for(i in 1:10){
    Lasso.cv<-cv.glmnet(x.train[[i]],y.train[[i]])
    bestlambda<-Lasso.cv$lambda.min
    Lasso.model<-Lasso.cv$glmnet.fit
    y.test.pred.Lasso<-predict(Lasso.model,newx=x.test[[i]],s=bestlambda)
    y.test.pred.total<-c(y.test.pred.total,y.test.pred.Lasso)
  }
  names(y.test.pred.total)<-tmp
  sum_y_top20 <- y.test.pred.total[sgRNA.top20.id] + sum_y_top20
  sum_y_total <- y.test.pred.total+sum_y_total
}
y.test.pred.top20 <- sum_y_top20/10
y.test.pred.total.Lasso.mean <- sum_y_total/10
PCC.Lasso.7mer.top20<-cor(y.test.total[sgRNA.top20.id],y.test.pred.top20)
PCC.Lasso.7mer.all <- cor(y.test.total,y.test.pred.total.Lasso.mean)
save(y.test.total,y.test.pred.total.Lasso.mean,PCC.Lasso.7mer.top20,PCC.Lasso.7mer.all,file="PCC.Lasso.7mer.Rdata")

sse = sum((y.test.total-y.test.pred.total.Lasso.mean)^2)
sst = sum((y.test.total-mean(y.test.total))^2)
R2_lasso = (1 - sse/sst)^2

MSE_lasso = mean((y.test.pred.total.Lasso.mean-y.test.total)^2)

##Perfrom Ridge Regression
sum_y_top20<-rep(0,20)
sum_y_total<-rep(0,980)
for(m in 1:10){
  cat("m=",m,'\t')
  y.test.pred.total<-c()
  for(i in 1:10){
    Ridge.cv<-cv.glmnet(x.train[[i]],y.train[[i]],alpha=0)
    bestlambda<-Ridge.cv$lambda.min
    Ridge.model<-Ridge.cv$glmnet.fit
    y.test.pred.Ridge<-predict(Ridge.model,newx=x.test[[i]],s=bestlambda)
    y.test.pred.total<-c(y.test.pred.total,y.test.pred.Ridge)
  }
  names(y.test.pred.total)<-tmp
  sum_y_top20 <- y.test.pred.total[sgRNA.top20.id] + sum_y_top20
  sum_y_total <- y.test.pred.total+sum_y_total
}
y.test.pred.top20 <- sum_y_top20/10
y.test.pred.total.Ridge.mean <- sum_y_total/10
PCC.Ridge.7mer.top20<-cor(y.test.total[sgRNA.top20.id],y.test.pred.top20)
PCC.Ridge.7mer.all <- cor(y.test.total,y.test.pred.total.Ridge.mean)
save(y.test.total,y.test.pred.total.Ridge.mean,PCC.Ridge.7mer.top20,PCC.Ridge.7mer.all,file="PCC.Ridge.7mer.Rdata")

sse = sum((y.test.total-y.test.pred.total.Ridge.mean)^2)
sst = sum((y.test.total-mean(y.test.total))^2)
R2_RR = (1 - sse/sst)^2

MSE_RR = mean((y.test.pred.total.Ridge.mean-y.test.total)^2)