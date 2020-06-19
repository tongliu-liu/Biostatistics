install.packages("randomForest")
library(Biostrings)
library(kknn)
library(pROC)
library(gmodels)
library(ggplot2)
library(randomForest)

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

data_frq <- rbind(neg_frq,pos_frq)
data_frq[is.na(data_frq)]<-0

val <- sample(1:nrow(data_frq), size = round(nrow(data_frq)*2/3), replace = FALSE)
data.learn <- data.frame(data_frq[val,])
data.test <- data.frame(data_frq[-val,])
tail(data.learn)
tail(data.test)
nrow(data.learn)
nrow(data.test)

# 随机森林
ntree <- c()
for (i in seq(1,100,5)){
rf_fit = randomForest(category~.,data=data.learn,mtry=1,importance=TRUE,ntree= i)
rf_fit$importance

pred = predict(rf_fit,data.test)
t = table(pred,data.test$category)
err.rate <- (sum(t)-sum(diag(t)))/sum(t)
ntree <- c(ntree,err.rate)
}
ntree_num <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
plot(ntree_num,ntree)

ntree_data <- cbind(ntree_num,ntree)
ntree_data <- data.frame(ntree_data)
names(ntree_data) <- c("ntree","err.rate")

ggplot(data=ntree_data, aes(ntree, y=err.rate)) + geom_line(color = "red") +
  geom_point(color = "red") +
  scale_color_discrete(guide = guide_legend(title = NULL)) + theme_minimal() +
  ggtitle("KNN")

set.seed(1234)
rf_fit = randomForest(category~.,data=data.learn,mtry=1,importance=TRUE,ntree= 400)
plot(rf_fit)