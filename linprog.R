library(Matrix)
#library(tidyverse)
#library(MASS)
#library(MVA)
#library(maps)
#library(MCMCpack)
#library(doBy)
#library(utilities)
library(irlba)
#library(scatterplot3d)
#library(ggplot2)
#library(reshape2)
#library(latex2exp)
#library(vegan3d)
#library(car)
#library(ggtext)


sigma=0.1
T=20
TT<-20
d<-4
beta=15.00
m_vec<-c(0.5,0.5,0.5,0.5)

n_vec<-seq(100,300,100)


tr_vec<-1:T


risk1_vec<-vector()
risk2_vec<-vector()
risk3_vec<-vector()
risk4_vec<-vector()

for(n in n_vec)
{
  
  loss1_vec<-vector()
  loss2_vec<-vector()
  loss3_vec<-vector()
  loss4_vec<-vector()
  
  for(trial in 1:T)
  {
    t<-runif(n,min=0,max=1)
    
    
    e<-rnorm(n,mean=0,sd=sigma)
    y<-vector()
    y<-beta*t+e
    
    X<-matrix(nrow=n,ncol=d)
    
    X<-cbind(t/2,t/2,t/2,t/2)
    
    P<-X%*%t(X)
    
    
    A<-matrix(nrow=n, ncol=n)
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        if(i==j)
        {
          u<-P[i,j]
          A[i,j]<-rbinom(1,1,u)
          
        }
        if(i<j)
        {
          u<-P[i,j]
          A[i,j]<-rbinom(1,1,u)
        }
        if(i>j)
        {
          A[i,j]<-A[j,i]
        }
      }
    }
    
    
    A_irlba<-irlba(A,4)
    X_hat_raw<-A_irlba$u%*%diag(A_irlba$d)^0.5
    
    vals<-svd(t(X_hat_raw)%*%X)
    W<-vals$u%*%t(vals$v)
    X_hat<-X_hat_raw%*%W
    
    t_hat<-X_hat%*%m_vec
    
    
    
    
    
    
    beta_true=(t(t) %*% y)/(t(t) %*% t)
    #beta1_vec=c(beta1_vec,beta_true)
    loss1=(beta-beta_true)^2
    
    beta_naive=(t(t_hat) %*% y)/(t(t_hat) %*% (t_hat))
    #beta2_vec=c(beta2_vec,beta_naive)
    loss2=(beta-beta_naive)^2
    
    delta<-matrix(0,d,d)
    
    
    for(j in 1:n)
    {
      delta<-delta+(X_hat[j,]%*%t(X_hat[j,]))
    }
    Delta<-(1.0/n)*delta
    
    
    incr_vec<-vector()
    for(i in 1:n)
    {
      sum_mat<-matrix(0,d,d)
      
      for(j in 1:n)
      {
        
        F1<-(t(X_hat[i,])%*%X_hat[j,])*(1-t(X_hat[i,])%*%X_hat[j,])
        f1<-as.double(F1)
        num<-f1*(X_hat[j,]%*%t(X_hat[j,]))
        sum_mat<-sum_mat+num
      }
      
      Sum_mat<-(1.0/n)*sum_mat
      
      incr<-t(m_vec)%*%solve(Delta)%*%Sum_mat%*%solve(Delta)%*%m_vec
      
      incr_vec<-c(incr_vec,incr)
      
    }
    
    sum_gamma<-mean(incr_vec)
    
    
    beta_adj_par_hat=(t(t_hat) %*% y)/((t(t_hat)%*%t_hat)-sum_gamma)
    #beta4_vec<-c(beta4_vec,beta_adj_par_hat)
    loss4=(beta-beta_adj_par_hat)^2
    
    loss1_vec<-c(loss1_vec,loss1)
    loss2_vec<-c(loss2_vec,loss2)
    #loss3_vec<-c(loss3_vec,loss3)
    loss4_vec<-c(loss4_vec,loss4)
    
    
  }
  
  risk1<-mean(loss1_vec)
  risk1_vec<-c(risk1_vec,risk1)
  
  risk2<-mean(loss2_vec)
  risk2_vec<-c(risk2_vec,risk2)
  
  risk4<-mean(loss4_vec)
  risk4_vec<-c(risk4_vec,risk4)
  
  dec<-c(n,risk1,risk2,risk4)
  print(dec)
}

save(risk4_vec,file="mse.RData")
