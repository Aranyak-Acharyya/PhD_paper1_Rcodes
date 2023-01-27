library(parallel)
library(iterators)
no_of_cores<-detectCores()
clust<-parallel::makeCluster(no_of_cores)
library(foreach)
library(doParallel)
registerDoParallel()


library(Matrix)
library(MASS)
library(irlba)
library(vegan3d)

e <- new.env()
e$libs <- c("irlba","Matrix","vegan3d",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))




M<-matrix(,ncol=3)

x_vec<-seq(-10,10,by=0.001)

#n<-1000
n_vec<-seq(5000,100000,5000)

d<-4

s<-1000



clusterExport(clust,list("d","s","x_vec"))

diff_vec<-vector()


for(n in n_vec)
{
  
  ll<-n/5
  l<-as.integer(ll)
  
  
  KK<-l/5
  K<-as.integer(KK)
  
  clusterExport(clust,list("n","l","K"))
  
  B<-foreach(trial=1:100,.combine='c') %dopar%
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
      
      zs<-z_hat[1:s]+0.5
      
      
      
      F_hat_org<-density(ts,bw="nrd0",adjust=1,
                         kernel = "epanechnikov",
                         n=512,from=-0.05,to=1.05)$y
      
      F_hat_iso<-density(zs,bw="nrd0",adjust=1,
                         kernel = "epanechnikov",
                         n=512,from=-0.05,to=1.05)$y
      
      a<-max(abs(F_hat_org-F_hat_iso))
      a
      
      
    }
  
  
  
  new_row<-c(length(B),n,mean(B))
  print(new_row)
  
  M<-rbind(M,new_row)
  
  diff_vec<-c(diff_vec,mean(B))
  
}

stopCluster(clust)

M<-M[-1,]
print(M)

df<-data.frame(M[,1],M[,2],M[,3])
print(df)

save(df,file="P3new2.RData")
