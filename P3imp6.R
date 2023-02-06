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

w_vec<-seq(-10,10,by=0.001)

#n<-1000
n_vec<-seq(1000,15000,1000)

d<-4

s1<-50
s2<-50



clusterExport(clust,list("d","s1","s2","w_vec"))

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
      ts1<-rbeta(s1,shape1 = 1.2,shape2 = 2.8)
      ts2<-rbeta(s2,shape=2.0,shape2=5.0)
      
      #generating pre-images for auxiliary latent positions
      tm<-runif(n-s1-s2,min=0,max=1)
      
      #combining to form set of all pre-images
      t<-c(ts1,ts2,tm)
      
      
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
      
      zs1<-unique(z_hat[1:s1])
      zs2<-unique(z_hat[(s1+1):(s1+s2)])
      
      Ft1<-ecdf(ts1)
      Ft2<-ecdf(ts2)
      Fz1<-ecdf(zs1)
      Fz2<-ecdf(zs2)
      
      m1<-max(abs(Ft1(w_vec)-Ft2(w_vec)))
      m2<-max(abs(Fz1(w_vec)-Fz2(w_vec)))
      
      q<-abs(m1-m2)
      q
      
    }
  
  
  
  
  
  new_row<-c(length(B),n,mean(B))
  print(new_row)
  
  M<-rbind(M,new_row)
  
  diff_vec<-c(diff_vec,mean(B))
  
}

stopCluster(clust)


M<-M[-1,]
df<-as.data.frame(M)
save(df,file="P3new6.RData")