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
library(vegan3d)
#library(car)




set.seed(12345)


s<-5
l_vec<-seq(100,500,100)
T<-100
d<-4

a<-2.0
b<-1.0



#probably need larger n for isomap on X_hat
risk_vec<-vector()
gl_min_vec<-vector()




for(l in l_vec)
{
  
 
  m<-l-s
  n<-(s+m)
  
  KK<-l^(0.5)
  K<-as.integer(KK)
  
  
  
  loss_vec<-vector()
  
  U<-vector()
  U_alt<-vector()
  
  qt_vec<-vector()
  
  for(trial in 1:T)
  {
    trial<-1
    ts<-runif(s,min=0,max=1)
    #pre-images of auxiliary latent positions
    tm<-runif(m,min=0,max=1)
    
    #forming entire set of pre-images by concatenation
    t<-c(ts,tm)
    
    
    
    
    
    #forming the matrix whose rows are true latent positions
    X<-cbind(t/2,t/2,t/2,t/2)
    
    
    
    #forming probability matrix
    P<-X%*%t(X)
    
    #forming adjacency matrix
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
    
    
    #finding adjacency spectral embedding of desired dimension
    A_irlba<-irlba(A,d)
    X_hat_raw<-A_irlba$u%*%diag(A_irlba$d)^0.5
    
    #finding consistent adjacency spectral estimates of latent positions
    vals<-svd(t(X_hat_raw)%*%X)
    W<-vals$u%*%t(vals$v)
    X_hat<-X_hat_raw%*%W
    
    
    
    XX_hat<-X_hat[1:l,]
    
    
    
    #1st step of dimension reduction of X_hat by isomap
    z_hat<-isomap(vegdist(XX_hat, method="euclidean"),
                      k=K,ndim=1,
                      path="shortest")$points
    
    #scaling the dimension-reduced 1d embeddings to estimate original scalar pre-images
    #t_ord_hatt<-t_ord_hat_raw+0.5  
    #t_hat<-round(c(scales::rescale(t_hat_raw,c(0,1))),4) 
    
    
    
    
    qt<-0
    for(i in 1:s)
    {
      for(j in 1:s)
      {
        qt<-qt+(abs(z_hat[i]-z_hat[j])-abs(t[i]-t[j]))^2
      }
    }
    
    
    qt_vec<-c(qt_vec,qt)
    
    #ordering the scaled embeddings once again
    #ts_hat<-t_hat[1:s]
    #t_bar<-mean(t_hat_raw)
    
    #loss<-max(abs(ts_hat-ts))
    #loss_vec<-c(loss_vec,loss)
    
    #estimating CDF with known scalar pre-images
    
    
  }
  
  gl_min<-mean(qt_vec)
  gl_min_vec<-c(gl_min_vec,gl_min)
  
  dec<-c(l,gl_min)
  print(dec)
  
  
  #dec<-c(n,risk)
  #print(dec)
  
  #u1<-min(U)-0.5
  #u2<-max(U)+0.5
  
  #u1_alt<-min(U_alt)-0.5
  #u2_alt<-max(U_alt)+0.5
  
  #hist(U,breaks=50,xlim=c(u1,u2))
  
  #hist(U_alt,breaks=50,xlim=c(u1_alt,u2_alt))
    
}

df<-data.frame(l_vec,gl_min_vec)
save(df,"unkn.RData")