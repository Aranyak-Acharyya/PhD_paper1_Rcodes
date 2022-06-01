library(Matrix)
library(tidyverse)
library(MASS)
library(MVA)
library(maps)
library(MCMCpack)
library(doBy)
library(utilities)
library(irlba)
library(scatterplot3d)
library(ggplot2)
library(reshape2)
library(latex2exp)
library(vegan3d)
library(car)




s<-3
n_vec<-seq(500,750,50)
T<-100
beta<-15.0

risk_vec<-vector()


for(n in n_vec)
{

m<-(n-s)
ts<-c(0.786,0.452,0.813)
loss_vec<-vector()



for(trial in 1:T)
{
  
  tm<-runif(m,min=0,max=1)
  t<-c(ts,tm)
  rvec<-rank(t)[1:s]
  
  t_ord<-sort(t,decreasing = FALSE)
  
  X<-cbind(t_ord^2,2*t_ord*(1-t_ord),(1-t_ord)^2)
  
  
  e<-rnorm(s,mean=0,sd=0.1)
  y<-beta*ts+e
  P<-X%*%t(X)
  
  
  A<-matrix(nrow=n,ncol=n)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      if(i<=j)
      {
        v<-P[i,j]
        A[i,j]<-rbinom(1,1,v)
      }
      if(i>j)
      {
        A[i,j]<-A[j,i]
      }
    }
  }
  
  
  A_irlba<-irlba(A,3)
  X_hat_raw<-A_irlba$u%*%diag(A_irlba$d)^0.5
  
  
  vals<-svd(t(X_hat_raw)%*%X)
  W<-vals$u%*%t(vals$v)
  X_hat<-X_hat_raw%*%W
  
  
  
  
  t_ord_hat_raw<-isomap(vegdist(X_hat, method="euclidean"),
                        k=20,ndim=1,
                        path="shortest")$points
  t_ord_hatt<-t_ord_hat_raw+0.5  
  t_ord_hat<-round(c(scales::rescale(t_ord_hatt,c(0,1))),4)  
  t_ord_Hat<-sort(t_ord_hat,decreasing = FALSE)
  
  ts_Hat<-t_ord_Hat[rvec]
  
 
  
  beta_true<-(t(ts)%*%y)/(t(ts)%*%ts)
  #loss1<-(beta_true-beta)^2
  #loss1_vec<-c(loss1_vec,loss1)
  
  beta_naive<-(t(ts_Hat)%*%y)/(t(ts_Hat)%*%ts_Hat)
  #loss2<-(beta_naive-beta)^2
  #loss2_vec<-c(loss2_vec,loss2)
  
  loss<-(beta_naive-beta_true)^2
  loss_vec<-c(loss_vec,loss)
}



risk<-mean(loss_vec)
risk_vec<-c(risk_vec,risk)



dec<-c(n,risk)
print(dec)
}


plot(n_vec,risk_vec,type="l",col="green")



x<-seq(0,1,0.1)
y<-5*x+rnorm(length(x),mean=0,sd=0.1)
plot(x,y,type="l")