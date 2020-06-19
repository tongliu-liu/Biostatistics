library(Biostrings)
library(kknn)
library(pROC)
library(gmodels)
library(ggplot2)
library(randomForest)
library(e1071)
library(caret)

# 计算程序运行时间
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
time1<-Sys.time()

setwd ( "E:/生物统计/结课论文数据/结课论文数据/Classification" )

# 载入数据
neg_data = readDNAStringSet("ABI5-neg.fasta")
pos_data = readDNAStringSet("ABI5-pos.fasta")
head(neg_data)
head(pos_data)

# 统计ATCG频率
A_neg_frq <- letterFrequency(neg_data,'A')/letterFrequency(neg_data,'AGCT')
T_neg_frq <- letterFrequency(neg_data,'T')/letterFrequency(neg_data,'AGCT')
C_neg_frq <- letterFrequency(neg_data,'C')/letterFrequency(neg_data,'AGCT')
G_neg_frq <- letterFrequency(neg_data,'G')/letterFrequency(neg_data,'AGCT')
A_pos_frq <- letterFrequency(pos_data,'A')/letterFrequency(pos_data,'AGCT')
T_pos_frq <- letterFrequency(pos_data,'T')/letterFrequency(pos_data,'AGCT')
C_pos_frq <- letterFrequency(pos_data,'C')/letterFrequency(pos_data,'AGCT')
G_pos_frq <- letterFrequency(pos_data,'G')/letterFrequency(pos_data,'AGCT')
A_neg_frq <- data.frame(A_neg_frq)
T_neg_frq <- data.frame(T_neg_frq)
C_neg_frq <- data.frame(C_neg_frq)
G_neg_frq <- data.frame(G_neg_frq)
A_pos_frq <- data.frame(A_pos_frq)
T_pos_frq <- data.frame(T_pos_frq)
C_pos_frq <- data.frame(C_pos_frq)
G_pos_frq <- data.frame(G_pos_frq)

neg_frq <- cbind(A_neg_frq,T_neg_frq,C_neg_frq,G_neg_frq)
pos_frq <- cbind(A_pos_frq,T_pos_frq,C_pos_frq,G_pos_frq)

category <- c(rep("neg",nrow(neg_frq)))
neg_frq <- cbind(neg_frq,category)
category <- c(rep("pos",nrow(pos_frq)))
pos_frq <- cbind(pos_frq,category)
head(neg_frq)
head(pos_frq)

# 规整全部数据
data_frq <- rbind(neg_frq,pos_frq)
data_frq[is.na(data_frq)]<-0


########KNN模型
# 模型训练，确定最佳k值、最佳kernel参数
data.kknn <- train.kknn(category~., data_frq, distance = 1, kernel = c("rectangular", "triangular", "epanechnikov", "optimal"))

# 不同k值对应预测错误率图
t = data.frame(data.kknn$MISCLASS[,1])
names(t) <- c("err.rate")
K <- c(1,2,3,4,5,6,7,8,9,10,11)
k <- data.frame(K)
err_rate <- cbind(K,t)
err_rate <- data.frame(err_rate)
class(err_rate)
ggplot(data=err_rate, aes(K, y=err.rate)) + geom_line(color = "red") +
  geom_point(color = "red") +
  scale_color_discrete(guide = guide_legend(title = NULL)) + theme_minimal() +
  ggtitle("KNN")

# 设置随机数
set.seed(7)
require(caret)
folds <- createFolds(data_frq$category,k=10)

# 十折交叉验证，训练模型
pre_knn_sum <- c()
fold_test_sum <- c()
err_rate_sum = 0
for(i in 1:10){
  fold_test <- data_frq[folds[[i]],]   #取folds[[i]]作为测试集
  fold_train <- data_frq[-folds[[i]],]
  data.kknn <- kknn(category~., fold_train,fold_test, distance = 1, k =  11, kernel = "rectangular")
  pre_knn <- fitted(data.kknn)
  pre_knn_sum <- c(pre_knn_sum,pre_knn)
  fold_test_sum <- c(fold_test_sum,fold_test$category)
  t = table(fold_test$category,pre_knn)
  err_rate_sum <- err_rate_sum+(sum(t)-sum(diag(t)))/sum(t)
  print((sum(t)-sum(diag(t)))/sum(t))
  
  # 绘制ROC曲线并计算AUC值
  knn_roc <- roc(fold_test$category,as.numeric(pre_knn))
  plot(knn_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='knn算法ROC曲线')
}
# 绘制平均AUC曲线并计算AUC值
knn_roc <- roc(fold_test_sum,as.numeric(pre_knn_sum))
plot(knn_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='knn算法ROC曲线')

# 求平均错误率
err_rate = err_rate_sum/10

##########逻辑回归算法
set.seed(7)
require(caret)
folds <- createFolds(data_frq$category,k=10)

err.rate.sum = 0
pre_logistic_sum <- c()
fold_test_sum <- c()
for(i in 1:10){
  
  fold_test <- data_frq[folds[[i]],]   #取folds[[i]]作为测试集
  fold_train <- data_frq[-folds[[i]],]   # 剩下的数据作为训练集
  
  # 逻辑回归
  logit.fit <- glm(category~.,
                   family = binomial(link = 'logit'),
                   data = fold_train)
  logit.predictions <- ifelse(predict(logit.fit) > 0,'neg', 'pos')
  pre_logistic<-as.numeric(predict(logit.fit,newdata=fold_test,type="response")>0.5)
  pre_logistic_sum <- c(pre_logistic_sum,pre_logistic)
  fold_test_sum <- c(fold_test_sum,fold_test$category)
  
  #输出混淆矩阵
  t = table(fold_test$category,pre_logistic) 
  err.rate.sum <- (sum(t)-sum(diag(t)))/sum(t) + err.rate.sum
  
  print((sum(t)-sum(diag(t)))/sum(t))
  pred <- roc(fold_test$category,pre_logistic)
  plot(pred, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='逻辑回归算法ROC曲线')
  
}
# 绘制平均AUC曲线并计算AUC值
logistic_roc <- roc(fold_test_sum,as.numeric(pre_logistic_sum))
plot(logistic_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='逻辑回归算法ROC曲线')

err.rate <- err.rate.sum/10

###########随机森林
set.seed(7)
require(caret)
folds <- createFolds(data_frq$category,k=10)

err.rate.sum = 0
pre_rf_sum <- c()
fold_test_sum <- c()
for(i in 1:10){
  
  fold_test <- data_frq[folds[[i]],]   #取folds[[i]]作为测试集
  fold_train <- data_frq[-folds[[i]],]   # 剩下的数据作为训练集
  
  # 随机森林
  rf_fit = randomForest(category~.,data=fold_train,mtry=1,importance=TRUE,ntree=50)
  
  pred = predict(rf_fit,fold_test)
  pre_rf_sum <- c(pre_rf_sum,pred)
  fold_test_sum <- c(fold_test_sum,fold_test$category)
  t = table(pred,fold_test$category)
  
  err.rate.sum <- (sum(t)-sum(diag(t)))/sum(t) + err.rate.sum
  
  print((sum(t)-sum(diag(t)))/sum(t))
  pred <- roc(fold_test$category,as.numeric(pred))
  plot(pred, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='随机森林算法ROC曲线')
  
} 
# 绘制平均AUC曲线并计算AUC值
rf_roc <- roc(fold_test_sum,as.numeric(pre_rf_sum))
plot(rf_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='随机森林算法ROC曲线')

err.rate <- err.rate.sum/10

###########SVM
set.seed(7)
require(caret)
folds <- createFolds(data_frq$category,k=10)

err.rate.sum = 0
pre_svm_sum <- c()
fold_test_sum <- c()
for(i in 1:10){
  
  fold_test <- data_frq[folds[[i]],]   #取folds[[i]]作为测试集
  fold_train <- data_frq[-folds[[i]],]   # 剩下的数据作为训练集
  
  # svm
  model = svm(category~.,data = fold_train, kernel="linear")
  
  pred = predict(model,fold_test)
  pre_svm_sum <- c(pre_svm_sum,pred)
  fold_test_sum <- c(fold_test_sum,fold_test$category)
  t = table(pred,fold_test$category)
  
  err.rate.sum <- (sum(t)-sum(diag(t)))/sum(t) + err.rate.sum
  
  print((sum(t)-sum(diag(t)))/sum(t))
  pred <- roc(fold_test$category,as.numeric(pred))
  plot(pred, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM算法ROC曲线')
  
} 
# 绘制平均AUC曲线并计算AUC值
svm_roc <- roc(fold_test_sum,as.numeric(pre_svm_sum))
plot(svm_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM算法ROC曲线')

err.rate <- err.rate.sum/10

#计算程序运行时间
f(time1)