rm(list=ls())
library(glmnet)
library(spatstat) #dummify
library(ncpen)
library(ncvreg)
###logistic
xy.df = read.csv("C:\\Users\\kis91\\Desktop\\new.sam.for.log.fit.csv")
levels(xy.df[,1]) = c("non-case","case")
tem = paste(names(xy.df)[-1],collapse = "+")
f.mod = paste(names(xy.df)[1],"~",tem)
f.fit = glm(f.mod,data=xy.df,family = "binomial")
r.mod = paste(names(xy.df)[1],"~1")
r.fit = glm(r.mod,data=xy.df,family = "binomial")
fit = step(r.fit,f.mod , dirdction="forward",k=2)
aic.fit = fit
tab = table(xy.df[,1],predict(fit,xy.df)>0)
sen = tab[2,2]/sum(tab[2,])
spc = tab[1,1]/sum(tab[1,])
acc = sum(diag(tab))/sum(tab)

x.mat = as.matrix(xy.df[,-1])
y.vec = dummify(as.vector(xy.df[,1]))[,2] 
fit.glmnet = glmnet(x.mat,y.vec,family = "binomial")
colSums(coef(fit.glmnet)!=0)
fit.cv = cv.glmnet(x.mat,y.vec,family = "binomial")
b.vec = coef(fit.cv)
b.vec[b.vec!=0]

fit.ncvreg = ncvreg(x.mat,y.vec,family = "binomial")
fit.cv.ncv = cv.ncvreg(x.mat,y.vec,family = "binomial")
b.vec = (coef(fit.cv.ncv))
b.vec[b.vec!=0]

