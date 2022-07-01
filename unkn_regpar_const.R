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


s<-5
l_vec<-seq(100,800,100)
T<-100
d<-4

a<-2.0
b<-1.0

alpha<-10.0
beta<-25.0



#probably need larger n for isomap on X_hat
risk_vec<-vector()
pred_diff_vec<-vector()




for(l in l_vec)
{
  
  #total no. of nodes of the graph
  n<-5*l
  
  #no. of auxiliary latent positions to be generated
  m<-(n-s)
  
  #determining nearest neighbour rule on constructed graph for isomap
  KK<-l/5
  K<-as.integer(KK)
  
  qt_vec<-vector()
  
  for(trial in 1:T)
  {
   
    #pre-images of regression-related latent positions
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
    
    
    #creating the regression model
    y<-alpha+beta*ts+rnorm(s,mean=0,sd=0.01)
    
    #hypothetical estimates of regression parameters
    # if true values of regressors were known    
    beta_true<-(cov(ts,y))/var(ts)
    alpha_true<-mean(y)-beta_true*mean(ts)
    
    
    
    #extracting rows on which isomap will run
    XX_hat<-X_hat[1:l,]
    
    
    
    #1st step of dimension reduction of X_hat by isomap
    z_hat<-isomap(vegdist(XX_hat, method="euclidean"),
                  k=K,ndim=1,
                  path="shortest")$points
    
    
    
    zs_hat<-z_hat[1:s]
    zm_hat<-z_hat[(s+1):l]
    
    #finding naive estimates of regression parameters based on
    # isomap embeddings
    beta_naive<-(cov(zs_hat,y))/var(zs_hat)
    alpha_naive<-mean(y)-beta_naive*mean(zs_hat)
    
    
    
    #Goal:to check if prediction errors of the two estimators match
    
    #generating new response variable for auxiliary latent positions
    ynew<-alpha+beta*t[(s+1):l]+rnorm(l-s,mean=0,sd=0.01)
    
    #estimating response variables by true estimates
    ynew_true<-alpha_true+beta_true*t[(s+1):l]
    
    
    #estimates of new response variable by isomap embedding estimates
    ynew_naive<-alpha_naive+beta_naive*zm_hat
    
    
   
    
    yy<-(ynew_naive-ynew_true)/ynew
    
    qt<-mean(yy^2)
    
    
    
    
    qt_vec<-c(qt_vec,qt)
    
    dec1<-c(trial,qt)
    print(dec1)
    
  }
  
  
  pred_diff<-mean(qt_vec)
  
  pred_diff_vec<-c(pred_diff_vec,pred_diff)
  
  dec<-c(l,pred_diff)
  print(dec)
  
}

plot(l_vec,pred_diff_vec,type="l")