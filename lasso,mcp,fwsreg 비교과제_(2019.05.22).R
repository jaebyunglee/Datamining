rm(list=ls())
library(glmnet)
library(spatstat) #dummify
library(ncpen)
library(ncvreg)
###logistic
xy.df = read.csv("C:\\Users\\kis91\\Desktop\\new.sam.for.log.fit.csv")
xy.df = xy.df[,1:10]

final.mat = matrix(NA,5,3)
rownames(final.mat) = c("fwreg","lasso","lasso-aic","mcp","mcp-aic")
colnames(final.mat) = c("acc","sen","spc")

fwreg.mat = NULL
lasso.mat = NULL
lasso.aic.mat = NULL
mcp.mat = NULL
mcp.aic.mat = NULL
iter = 100
n.vec = c(9,1) #split rate
b = cumsum(n.vec)/sum(n.vec)
pos = xy.df[,1]=="case"
ntb = ceiling(sum(pos)*b); nfb = ceiling(sum(!pos)*b)
for(i in 1:iter){
  print(i)
  nta = sample(1:sum(pos)); nfa = sample(1:sum(!pos)) 
  t.trid = nta[1:ntb[1]] ; t.teid = nta[(ntb[1]+1):ntb[2]]
  f.trid = nfa[1:nfb[1]]; f.teid = nfa[(nfb[1]+1):nfb[2]]
  test = rbind(xy.df[pos,][t.teid,],xy.df[!pos,][f.teid,])
  train = rbind(xy.df[pos,][t.trid,],xy.df[!pos,][f.trid,])
  train[,1] = train[,1]=="case"
  test[,1] = test[,1]=="case"
  
  ### fwreg
  tem = paste(names(xy.df)[-1],collapse = "+")
  r.mod = paste(names(xy.df)[1],"~1")
  f.mod = paste(names(xy.df)[1],"~",tem)
  r.fit = glm(r.mod,data=train,family = "binomial")
  fw.fit = step(r.fit,f.mod,direction="forward",k=2,trace = "F")
  tab = table(test[,1],predict(fw.fit,test)>0)
  acc = sum(diag(tab))/sum(tab)
  sen = tab[2,2]/sum(tab[2,])
  spc = tab[1,1]/sum(tab[1,])
  fwreg.mat = rbind(fwreg.mat,cbind(acc,sen,spc))
  
  ### lasso
  try.vec = train[,1]; trx.mat = as.matrix(train[,-1])
  tey.vec = test[,1]; tex.mat = as.matrix(test[,-1])
  cv.glm = cv.glmnet(trx.mat,try.vec,family = "binomial")
  tab = table(test[,1],predict(cv.glm,tex.mat)>0)
  acc = sum(diag(tab))/sum(tab)
  sen = tab[2,2]/sum(tab[2,])
  spc = tab[1,1]/sum(tab[1,])
  lasso.mat = rbind(lasso.mat,cbind(acc,sen,spc))
  
  ### lasso - aic
  tl = cv.glm$glmnet.fit$nulldev-deviance(cv.glm$glmnet)
  k = cv.glm$glmnet.fit$df
  n = cv.glm$glmnet.fit$nobs
  aic.vec = -tl+2*k+2*k*(k+1)/(n-k-1)
  # aic.vec = NULL
  # for(i in 1:length(cv.glm$glmnet.fit$lambda)){
  #   aic.vec = c(aic.vec,aic.fun(coef(cv.glm$glmnet)[,i],cbind(1,trx.mat),try.vec))
  # }
  aic.pos = which.min(aic.vec)
  tab = table(test[,1],predict(cv.glm$glmnet.fit,tex.mat)[,aic.pos]>0)
  acc = sum(diag(tab))/sum(tab)
  sen = tab[2,2]/sum(tab[2,])
  spc = tab[1,1]/sum(tab[1,])
  lasso.aic.mat = rbind(lasso.aic.mat,cbind(acc,sen,spc))
  
  ### mcp
  cv.ncv = cv.ncvreg(trx.mat,try.vec,family = "binomial",penalty = "MCP")
  tab = table(test[,1],predict(cv.ncv,tex.mat)>0)
  acc = sum(diag(tab))/sum(tab)
  sen = tab[2,2]/sum(tab[2,])
  spc = tab[1,1]/sum(tab[1,])
  mcp.mat = rbind(mcp.mat,cbind(acc,sen,spc))
  
  ### mcp - aic
  tl = -2*cv.ncv$fit$loss
  k = colSums(coef(cv.ncv$fit)!=0)
  n = length(try.vec)
  aic.vec = -tl+2*k+2*k*(k+1)/(n-k-1)
  aic.pos = which.min(aic.vec)
  tab = table(test[,1],predict(cv.ncv$fit,tex.mat)[,aic.pos]>0)
  acc = sum(diag(tab))/sum(tab)
  sen = tab[2,2]/sum(tab[2,])
  spc = tab[1,1]/sum(tab[1,])
  mcp.aic.mat = rbind(mcp.aic.mat,cbind(acc,sen,spc))
}

final.mat[1,] = colMeans(fwreg.mat)
final.mat[2,] = colMeans(lasso.mat)
final.mat[3,] = colMeans(lasso.aic.mat)
final.mat[4,] = colMeans(mcp.mat)
final.mat[5,] = colMeans(mcp.aic.mat)
final.mat
