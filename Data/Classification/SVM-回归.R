install.packages("glmnet")
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

library(e1071)
rbf.par=c( c(0.00001,0.00002,0.00005,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.1,0.2,0.5,1), seq(5,100,20),seq(100,1000,200))
PCC.SVR.RBF.top20<-matrix(,2,length(rbf.par))
PCC.SVR.RBF.top20[1,]<-rbf.par
PCC.SVR.RBF.all<-matrix(,2,length(rbf.par))
PCC.SVR.RBF.all[1,]<-rbf.par
y.test.pred.SVR.total<-list()
for(m in 1:length(rbf.par)){
  cat("m=",m)
  y.test.pred.total<-c()
  SVR.kernel="radial"
  rbf.kpar=rbf.par[m]
  for(i in 1:10){
    SVR.model<-svm(x.train[[i]],y.train[[i]],kernel=SVR.kernel,gamma=rbf.kpar);
    y.test.pred.SVR<-predict(SVR.model,x.test[[i]]);
    y.test.pred.total<-c(y.test.pred.total,y.test.pred.SVR);
  }
  names(y.test.pred.total)<-tmp
  PCC.SVR.RBF.top20[2,m]<-cor(y.test.total[sgRNA.top20.id],y.test.pred.total[sgRNA.top20.id])
  PCC.SVR.RBF.all[2,m] <- cor(y.test.total,y.test.pred.total)
  y.test.pred.SVR.total[[m]] <- y.test.pred.total
}
save(PCC.SVR.RBF.top20,PCC.SVR.RBF.all,y.test.pred.SVR.total,y.test.total,file="SVR_RBF_PCC_7mer_100bp.Rdata")