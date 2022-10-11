library(parallel)
library(iterators)
no_of_cores<-detectCores()
clust<-parallel::makeCluster(no_of_cores)
library(foreach)
library(doParallel)
registerDoParallel()



library(Matrix)
library(irlba)
library(princurve)
library(vegan3d)


e <- new.env()
e$libs <- c("irlba","Matrix","princurve","vegan3d",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))


s<-5
l_vec<-seq(100,1000,100)
T<-100
d<-4

clusterExport(clust,list("s","d"))



gl_min_vec<-vector()



for(l in l_vec)
{
  n<-5*l
  m<-(n-s)
  
  
  KK<-l/5
  K<-as.integer(KK)
  
  
  clusterExport(clust,list("l","n","m","K"))
  
  L<-foreach(i=1:100,.combine = 'c') %dopar%
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
      pvec<-c(P[lower.tri(P,diag=TRUE)])
      avec<-rbinom(length(pvec),1,pvec)
      A<-matrix(nrow=nrow(P),ncol=ncol(P))
      A[lower.tri(A,diag=TRUE)]<-avec
      A[upper.tri(A)]<-t(A)[upper.tri(A)]
      
      
      
      
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
      
      qt
      
    
      
    }
  
  gl_min<-mean(L)
  
  gl_min_vec<-c(gl_min_vec,gl_min)
  
  dec<-c(n,l,gl_min,length(L))
  print(dec)
  
  
  
}

n_vec<-5*l_vec
plot(n_vec,gl_min_vec,type="l")

