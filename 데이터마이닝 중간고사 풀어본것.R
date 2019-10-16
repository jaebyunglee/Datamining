rm(list=ls())
xy.df = read.csv("C:\\Users\\kis91\\Desktop\\old.sam.for.reg.fit.csv")
head(xy.df)

#A
full.fit = lm(sensitivity~.,data=xy.df)
null.fit = lm(sensitivity~type+pressure,data=xy.df)
A=step(null.fit, direction="forward",scope=list(lower=null.fit,upper=full.fit),k=2)

#B - 1
fit1 = lm(sensitivity~.,data=xy.df[,c(1,2,3,4)])
abs(coef(fit1)[5])

#B - 2
fit2 = lm(sensitivity~.,data=xy.df[,c(1,2,3,5)])
abs(coef(fit2)[5])

#B - 3
b.vec = c()
for(i in 4:21){
  fit = lm(sensitivity~.,data=xy.df[,c(1,2,3,i)])
  b.vec = c(b.vec,abs(coef(fit)[5]))
}
rank.vec = data.frame(b.vec,rank(-b.vec))

#B - 4
fit3 = lm(sensitivity~.,data=xy.df[,c(1,2,3,4)])
AIC(fit3)

#B - 5
fit4 = lm(sensitivity~.,data=xy.df[,c(1,2,3,4,10)])
AIC(fit4)

#B - 6
aic.vec = rep(NA,18)
c.set = 4:21
for(i in 1:18){
  pos = rank.vec[,2]<=i
  fit5 = lm(sensitivity~.,data=xy.df[,c(1,2,3,c.set[pos])])
  aic.vec[i]=AIC(fit5)
}
#B - 7
b.opt=which.min(aic.vec)

#C - 1
m = dim(xy.df)[1]
c.mat = matrix(NA,18,50)
rownames(c.mat) = names(b.vec)
id = sample(1:m)[1:(m*0.7)]
cb.vec = c()
for(i in 4:21){
  fit = lm(sensitivity~.,data=xy.df[id,c(1,2,3,i)])
  cb.vec = c(cb.vec,abs(coef(fit)[5]))
}
rank(-cb.vec)

#C - 2
for(j in 1:50){
  id = sample(1:dim(xy.df)[1],dim(xy.df)[1],replace = T)
  cb.vec = c()
  for(i in 4:21){
    fit = lm(sensitivity~.,data=xy.df[id,c(1,2,3,i)])
    cb.vec = c(cb.vec,abs(coef(fit)[5]))
  }
  c.mat[,j] = rank(-cb.vec)
}
mrank.vec = rank(rowMeans(c.mat))

#C - 3
fit6 = lm(sensitivity~.,data=xy.df[,c(1,2,3,4)])
AIC(fit6)

#C - 4
fit7 = lm(sensitivity~.,data=xy.df[,c(1,2,3,4,10)])
AIC(fit7)

#C - 5
c.aic.vec = rep(NA,18)
c.set = 4:21
for(i in 1:18){
  pos = mrank.vec<=i
  fit5 = lm(sensitivity~.,data=xy.df[,c(1,2,3,c.set[pos])])
  c.aic.vec[i]=AIC(fit5)
}
c.aic.vec

#C - 6
c.opt=which.min(c.aic.vec)

#D - 1
trid = sample(1:dim(xy.df)[1])[1:800]
A.fit = lm(sensitivity~.,data = A$model[trid,])
mean((xy.df[-trid,1]-predict(A.fit,newdata = A$model[-trid,]))^2)
B.fit = lm(sensitivity~.,data = xy.df[trid,c(1,2,3,4:(b.opt+3))])
mean((xy.df[-trid,1]-predict(B.fit,newdata = xy.df[-trid,]))^2)
C.fit = lm(sensitivity~.,data = xy.df[trid,c(1,2,3,4:(c.opt+3))])
mean((xy.df[-trid,1]-predict(C.fit,newdata = xy.df[-trid,]))^2)



#D - 2
D.mat = matrix(NA,50,3)
for(i in 1:50){
  trid = sample(1:dim(xy.df)[1])[1:800]
  A.fit = lm(sensitivity~.,data = A$model[trid,])
  D.mat[i,1]=mean((xy.df[-trid,1]-predict(A.fit,newdata = A$model[-trid,]))^2)
  B.fit = lm(sensitivity~.,data = xy.df[trid,c(1,2,3,4:(b.opt+3))])
  D.mat[i,2]=mean((xy.df[-trid,1]-predict(B.fit,newdata = xy.df[-trid,]))^2)
  C.fit = lm(sensitivity~.,data = xy.df[trid,c(1,2,3,4:(c.opt+3))])
  D.mat[i,3]=mean((xy.df[-trid,1]-predict(C.fit,newdata = xy.df[-trid,]))^2)
}
boxplot(D.mat)
