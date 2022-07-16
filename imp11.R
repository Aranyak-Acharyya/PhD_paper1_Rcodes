library(parallel)
library(iterators)
no_of_cores<-detectCores()
clust<-parallel::makeCluster(no_of_cores)
library(foreach)
library(doParallel)
registerDoParallel()

library(Matrix)
library(irlba)
library(vegan3d)


e <- new.env()
e$libs <- c("irlba","Matrix","vegan3d",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))









d<-4

s<-50

n_vec<-seq(1000,50000,1000)



alpha<-2.0
beta<-5.0
sig_ep<-0.01

clusterExport(clust,list("d","s","alpha","beta","sig_ep"))

RP<-matrix(,ncol=4)


for(n in n_vec)
{
  
  ll<-n/5
  l<-as.integer(ll)
  
  
  KK<-l/5
  K<-as.integer(KK)
  
  clusterExport(clust,list("n","l","K"))
  
  
  opts<-list(preschedule=FALSE)
  
  B<-foreach(trial=1:250,.combine='rbind',.options.multicore=opts) %dopar%
    {
      #generating regressors
      ts<-runif(s,min=0,max=1)
      
      #generating pre-images for auxiliary latent positions
      tm<-runif(n-s,min=0,max=1)
      
      #combining to form set of all pre-images
      t<-c(ts,tm)
      
      
      #forming the matrix whose rows are latent positions
      X<-cbind(t/2,t/2,t/2,t/2)
      
      #forming the probability matrix
      P<-X%*%t(X)
      
      #generating the adjacency matrix
      pvec<-c(P[lower.tri(P,diag=TRUE)])
      avec<-rbinom(length(pvec),1,pvec)
      A<-matrix(nrow=nrow(P),ncol=ncol(P))
      A[lower.tri(A,diag=TRUE)]<-avec
      A[upper.tri(A)]<-t(A)[upper.tri(A)]
      
      #finding ASE estimates
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
      
      #scaling the dimension-reduced 1d embeddings to estimate original scalar pre-images
      #t_ord_hatt<-t_ord_hat_raw+0.5  
      #t_hat<-round(c(scales::rescale(t_hat_raw,c(0,1))),4) 
      
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
      qt1<-(ynew_naive-ynew_true)^2
      
      qt2<-(ynew_true-ynew)^2
      
      qt3<-(ynew_naive-ynew)^2
      
      dec<-c(qt1,qt2,qt3)
      dec
      
    }
  
  new_row<-c(n,apply(B,2,mean))
  print(new_row)
  RP<-rbind(RP,new_row)
}

stopCluster(clust)

RP<-RP[-1,]

df<-data.frame(RP)
save(df,file="new11.RData")





