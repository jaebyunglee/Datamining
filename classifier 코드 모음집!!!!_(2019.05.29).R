rm(list=ls())

#install.packages("nnet")
#install.packages("ncpen")
#install.packages("glmnet")
#install.packages("robustreg")
#install.packages("spatstat")
#install.packages("class")
#install.packages("naivebayes")
#install.packages("rpart") 
#install.packages("rpart.plot")
#install.packages("e1071")
#install.packages("ipred")
#install.packages("randomForest") 
#install.packages("gbm")


### data analysis
work.dir = "D:\\S Kwon\\A Lecture\\A3 Data Mining\\A Sunghoon Kwon\\2019 01 Under Graduate\\"
data = read.csv("C:\\Users\\kis91\\Desktop\\new.sam.for.log.fit.csv")
dim(data)
# loading samples
xy.df = data[1:100,1:30]
nxy.df = data[101:145,1:30]
xy.df[1:5,1:5] 

# level check 
levels(xy.df[,1]) # the first level, case, becomes the baseline or y=0
xy.df[,1] = relevel(xy.df[,1],ref="non-case")
nxy.df[,1] = relevel(nxy.df[,1],ref="non-case")
levels(xy.df[,1])

# matrix form
library(spatstat)
xy.mat = dummify(xy.df); nxy.mat = dummify(nxy.df) 
xy.mat[1:5,1:5] 

xy.mat = xy.mat[,-1]; nxy.mat = nxy.mat[,-1]
xy.mat[1:5,1:5] 

# matrix for comparison 
c.mat = matrix(0,2,11)
colnames(c.mat) = c("logistic","lasso","scad"
                    ,"knn","naive","tree","svm","lda"
                    ,"bagging","r-forest","boosting")
rownames(c.mat) = c("training error", "test error") 

########################
### logistic regression
########################

# logistic 
log.fit = glm(disease~.,data=xy.df,family="binomial")
tab = table(xy.df[,1],predict(log.fit,newdata= xy.df[,-1])>0)
tr.err = 1-sum(diag(tab))/sum(tab)
tab = table(nxy.df[,1],predict(log.fit,newdata= nxy.df[,-1])>0)
ts.err = 1-sum(diag(tab))/sum(tab)
c.mat[,1] = c(tr.err,ts.err)

################################################
### lasso + logistic regression  
################################################

library(glmnet)

set.seed(1)
lasso.fit = cv.glmnet(x=xy.mat[,-1],y=xy.mat[,1],family="binomial")
tab = table(xy.df[,1],predict(lasso.fit,newx=xy.mat[,-1])>0)
tr.err = 1-sum(diag(tab))/sum(tab)
tab = table(nxy.df[,1],predict(lasso.fit,newx=nxy.mat[,-1])>0)
ts.err = 1-sum(diag(tab))/sum(tab)
c.mat[,2] = c(tr.err,ts.err)

################################################
### scad + logistic regression 
################################################

library(ncpen)

set.seed(1)
scad.fit = cv.ncpen(x.mat=xy.mat[,-1],y.vec=xy.mat[,1],family="binomial")
opt = which.min(scad.fit$like)
tab = table(xy.mat[,1],predict(scad.fit$ncpen.fit,new.x.mat=xy.mat[,-1])[,opt])
tr.err = 1-sum(diag(tab))/sum(tab)
tab = table(nxy.mat[,1],predict(scad.fit$ncpen.fit,new.x.mat=nxy.mat[,-1])[,opt])
ts.err = 1-sum(diag(tab))/sum(tab)
c.mat[,3] = c(tr.err,ts.err)

################################################
### knn 
################################################

library(class)

# choice of k???
set.seed(1)
k.grid = 50
cv.err.vec = rep(0,k.grid)
for(k in 1:k.grid){
  knn.fit = knn.cv(train=xy.df[,-1],cl=xy.df[,1],k=k)
  tab = table(xy.df[,1],knn.fit)
  cv.err.vec[k] = 1-sum(diag(tab))/sum(tab)
}
plot(cv.err.vec)
opt = which.min(cv.err.vec)
print(opt)

# fit with the optimal k 
knn.fit = knn(train=xy.df[,-1],test=xy.df[,-1],cl=xy.df[,1],k=opt)
tab = table(xy.df[,1],knn.fit)
tr.err = 1-sum(diag(tab))/sum(tab)
knn.fit = knn(train=xy.df[,-1],test=nxy.df[,-1],cl=xy.df[,1],k=opt)
tab = table(nxy.df[,1],knn.fit)
ts.err = 1-sum(diag(tab))/sum(tab)
c.mat[,4] = c(tr.err,ts.err)

################################################
### naivebayes
################################################

library(naivebayes)

na.fit = naive_bayes(x=xy.mat[,-1],y=xy.mat[,1])
tab = table(xy.mat[,1],predict(na.fit,newdata=xy.mat[,-1]))
tr.err = 1-sum(diag(tab))/sum(tab)
tab = table(nxy.mat[,1],predict(na.fit,newdata=nxy.mat[,-1]))
ts.err = 1-sum(diag(tab))/sum(tab)
c.mat[,5] = c(tr.err,ts.err)

################################################
### tree: CART
################################################

library(rpart)
library(rpart.plot)

cart.fit = rpart(disease~.,data=xy.df,method="class")
rpart.plot(cart.fit)
cart.fit$variable.importance
tab = table(xy.df[,1],predict(cart.fit,newdata=xy.df[,-1],type="class"))
tr.err = 1-sum(diag(tab))/sum(tab)
tab = table(nxy.df[,1],predict(cart.fit,newdata=nxy.df[,-1],type="class"))
ts.err = 1-sum(diag(tab))/sum(tab)
c.mat[,6] = c(tr.err,ts.err)

################################################
### SVM
################################################

library(e1071)

svm.fit = tune.svm(disease~.,data=xy.df,gamma=2^(-5:5),cost=2^(-5:5),type="C")
names(svm.fit)
opt = svm.fit$best.parameters
svm.fit = svm(disease~.,data=xy.df,type="C",gamma=opt[1],cost=opt[2]) # radial basis kernel
tab = table(predict(svm.fit,newdata=xy.df[,-1]),xy.df[,1])
tr.err = 1-sum(diag(tab))/sum(tab)
tab = table(predict(svm.fit,newdata=nxy.df[,-1]),nxy.df[,1])
ts.err = 1-sum(diag(tab))/sum(tab)
c.mat[,7] = c(tr.err,ts.err)


################################################
### LDA linear dicriminant analysis
################################################

library(MASS)
lda.fit = lda(disease~.,data=xy.df)
names(lda.fit)
tab = table(xy.df[,1],predict(lda.fit,newdata=xy.df[,-1])$class)
tr.err = 1-sum(diag(tab))/sum(tab)
tab = table(nxy.df[,1],predict(lda.fit,newdata=nxy.df[,-1])$class)
ts.err = 1-sum(diag(tab))/sum(tab)
c.mat[,8] = c(tr.err,ts.err)

################################################
### bagging with trees
################################################

library(ipred)

# choice of nbagg 
set.seed(1)
n.grid = 50
boot.err.vec = rep(0.5,n.grid)
for(k in 10:n.grid){
  bag.fit = bagging(disease~.,data=xy.df,nbagg=k,coob=T)
  boot.err.vec[k] = bag.fit$err
}
plot(boot.err.vec)
opt = which.min(boot.err.vec)

# fit bagging 
bag.fit = bagging(disease~.,data=xy.df,nbagg=opt)
bag.fit$mtrees[[1]]
bag.fit$mtrees[[25]]
length(bag.fit$mtrees)
tab = table(xy.df[,1],predict(bag.fit,newdata=xy.df[,-1]))
tab; levels(predict(bag.fit,newdata=xy.df[,-1]))
tr.err = sum(diag(tab))/sum(tab)
tab = table(nxy.df[,1],predict(bag.fit,newdata=nxy.df[,-1]))
ts.err = sum(diag(tab))/sum(tab)
c.mat[,9] = c(tr.err,ts.err) 
c.mat

################################################
### random Froest
################################################
library(randomForest)

set.seed(1)
randomForest(disease~.,data=xy.df)
set.seed(1)
randomForest(disease~.,data=xy.df,mtry=3) #default = 3
set.seed(1)
randomForest(disease~.,data=xy.df,mtry=4) 
set.seed(1)
randomForest(disease~.,data=xy.df,mtry=4,nodesize=1) #default = 1
set.seed(1)
randomForest(disease~.,data=xy.df,mtry=4,nodesize=2) #default = 1

set.seed(1)
rf.fit = rfcv(trainx=xy.df[,-1],trainy=xy.df[,1],cv.fold=5,nodesize=2)
names(rf.fit)
rf.fit$n.var
rf.fit$error.cv
opt = rf.fit$n.var[which.min(rf.fit$error.cv)]

rf.fit = randomForest(disease~.,data=xy.df,mtry=opt,nodesize=2) #default = 1
tab = table(xy.df[,1],predict(rf.fit,newdata=xy.df[,-1]))
tr.err = 1-sum(diag(tab))/sum(tab)
tab = table(nxy.df[,1],predict(rf.fit,newdata=nxy.df[,-1]))
ts.err = 1-sum(diag(tab))/sum(tab)

c.mat[,10] = c(tr.err,ts.err)



################################################
### boosting 
################################################

library(gbm)

# data processing 
set.seed(1)
txy.df = xy.df; tnxy.df = nxy.df  
txy.df[,1] = txy.df[,1]=="case"
tnxy.df[,1] = tnxy.df[,1]=="case"

# example 
gb.fit = gbm(disease~.,data=txy.df,distribution="adaboost",n.trees = 10)
gb.fit = gbm(disease~.,data=txy.df,distribution="adaboost",n.trees = 100)
#set.seed(1)
#gb.fit = gbm.fit(x=xy.mat[,-1],y=xy.mat[,1],distribution="bernoulli",n.trees = 10)

set.seed(1)
gb.fit = gbm(disease~.,data=txy.df,distribution="adaboost",n.trees=100,cv.folds=5)
opt = gbm.perf(gb.fit,method="cv",plot.it=TRUE)
tab = table(txy.df[,1],predict(gb.fit,newdata=txy.df[,-1],n.trees=opt)>0)
tr.err = 1-sum(diag(tab))/sum(tab)
tab = table(tnxy.df[,1],predict(gb.fit,newdata=tnxy.df[,-1],n.trees=opt)>0)
ts.err = 1-sum(diag(tab))/sum(tab)
c.mat[,11] = c(tr.err,ts.err)
c.mat
