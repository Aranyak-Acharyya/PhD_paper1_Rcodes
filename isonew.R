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


s<-5
l_vec<-seq(300,500,100)
T<-100
d<-4

a<-2.0
b<-1.0



#probably need larger n for isomap on X_hat
risk_vec<-vector()
gl_min_vec<-vector()




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
    
    #scaling the dimension-reduced 1d embeddings to estimate original scalar pre-images
    #t_ord_hatt<-t_ord_hat_raw+0.5  
    #t_hat<-round(c(scales::rescale(t_hat_raw,c(0,1))),4) 
    
    #pairwise difference matrix of regressors
    M_t<-matrix(t[1:s],nrow=s,ncol=s,byrow=FALSE)
    
    #pairwise difference matrix of isomap embeddings
    M_z_hat<-matrix(z_hat[1:s],nrow=s,ncol=s,byrow=FALSE)
    
    #maximum distance between pairwise differences 
    D<-abs(M_t-t(M_t))-abs(M_z_hat-t(M_z_hat))
    qt<-norm(D,type="M")
    
    
    
    
    
    #qt_vec<-c(qt_vec,qt)
    
    #dec1<-c(trial,qt)
    #print(dec1)
    
  }
  
  
  gl_min<-mean(qt_vec)
  
  gl_min_vec<-c(gl_min_vec,gl_min)
  
  dec<-c(l,gl_min)
  print(dec)
  
  
}


#print(gl_min_vec)
#l_new_vec<-l_vec[1:length(gl_min_vec)]
plot(l_new_vec,gl_min_vec,type="l")


df<-data.frame(l_new_vec,gl_min_vec)
save(df,file = "unkntest2.RData")