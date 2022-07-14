library(Matrix)
#library(tidyverse)
library(MASS)
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


Rcpp::cppFunction("
                  NumericMatrix generateAdjacencyMatrix(NumericMatrix pMatrix) {
                  
                  int n = pMatrix.cols();
                  NumericMatrix A(n,n);
                  for(int i = 0; i < n; i ++) {
                  for (int j = i + 1; j < n; j++) {
                  A(i,j) = (int)(rand()%100 < (pMatrix(i,j)* 100));
                  A(j,i) = A(i,j);
                  }
                  }
                  return A;
                  }
                  ")






set.seed(12345)


s<-20
l_vec<-seq(1000,2000,100)
T<-500
d<-4

alpha<-2.0
beta<-5.0
sig_ep<-0.01



#probably need larger n for isomap on X_hat

n_vec<-vector()
pred_diff_vec<-vector()



for(l in l_vec)
{
  
  n<-5*l
  m<-(n-s)
  
  
  KK<-l/5
  K<-as.integer(KK)
  
  
  
  loss_vec<-vector()
  
  qt_vec<-vector()
  
  for(trial in 1:T)
  {
    
    #creating pre-images of regression-related latent positions
    ts<-runif(s,min=0,max=1)
    
    #pre-images of auxiliary latent positions
    tm<-runif(m,min=0,max=1)
    
    #forming entire set of pre-images by concatenation
    t<-c(ts,tm)
    
    
    
    
    
    #forming the matrix whose rows are true latent positions
    X<-cbind(t/2,t/2,t/2,t/2)
    
    
    
    #forming probability matrix
    P<-X%*%t(X)
    
    #generating adjacency matrix
    A<-generateAdjacencyMatrix(P)
    
    
    
    
    #finding adjacency spectral embedding of desired dimension
    A_irlba<-irlba(A,d)
    X_hat_raw<-A_irlba$u%*%diag(A_irlba$d)^0.5
    
    #finding consistent adjacency spectral estimates of latent positions
    vals<-svd(t(X_hat_raw)%*%X)
    W<-vals$u%*%t(vals$v)
    X_hat<-X_hat_raw%*%W
    
    
    #extracting rows on which isomap will run
    XX_hat<-X_hat[1:l,]
    
    
    
    #1st step of dimension reduction of X_hat by isomap
    z_hat<-isomap(vegdist(XX_hat, method="euclidean"),
                  k=K,ndim=1,
                  path="shortest")$points
    
    zs<-z_hat[1:s]
    
    
    #generating auxiliary responses
    ys<-alpha+beta*ts+rnorm(s,mean=0,sd=sig_ep)
    
    #true regression parameter estimates
    beta_true<-cov(ys,ts)/var(ts)
    alpha_true<-mean(ys)-beta_true*mean(ts)
    
    #naive regression parameter estimates
    beta_naive<-cov(zs,ys)/var(zs)
    alpha_naive<-mean(ys)-beta_naive*mean(zs)
    
    #Goal: to check if naive and true estimators predict equally well
    
    tnew<-t[s+1]
    znew<-z_hat[s+1]
    
    #generating auxiliary response variables
    ynew<-alpha+beta*tnew+rnorm(1,mean=0,sd=sig_ep)
    
    #predicting by true regressors
    ynew_true<-alpha_true+beta_true*tnew
    
    #predicting by isomap embeddings
    ynew_naive<-alpha_naive+beta_naive*znew
    
    #proximity of predicted responses
    qt<-abs(ynew_naive-ynew_true)
    
    
    
    
    
    
    
    
    
    qt_vec<-c(qt_vec,qt)
    
    #dec1<-c(trial,qt)
    #print(dec1)
    
  }
  
  
  pred_diff<-mean(qt_vec)
  
  pred_diff_vec<-c(pred_diff_vec,pred_diff)
  
  n_vec<-c(n_vec,n)
  #dec<-c(l,gl_min)
  #print(dec)
  
  
}


#print(gl_min_vec)
#l_new_vec<-l_vec[1:length(gl_min_vec)]
#plot(l_new_vec,gl_min_vec,type="l")


df<-data.frame(l_vec,n_vec,pred_diff_vec)
save(df,file = "isopred.RData")



