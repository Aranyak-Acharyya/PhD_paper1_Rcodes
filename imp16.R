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

s<-15

n_vec<-seq(200,800,25)



alpha<-2.0
beta<-5.0
sig_ep<-0.01

u_vec<-seq(0,100,0.00001)

clusterExport(clust,list("d","s","alpha","beta","sig_ep"))

RP<-matrix(,ncol=4)

max_diff_vec<-vector()

for(n in n_vec)
{
  
  ll<-n/5
  l<-as.integer(ll)
  
  
  KK<-l/5
  K<-as.integer(KK)
  
  clusterExport(clust,list("n","l","K"))
  
  
  opts<-list(preschedule=FALSE)
  
  B<-foreach(trial=1:100,.combine='rbind',.options.multicore=opts) %dopar%
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
      
      #true predicted response
      ys_pred_true<-alpha_true+beta_true*ts
      
      #naive predicted responses
      ys_pred_naive<-alpha_naive+beta_naive*zs
      
      #true F-statistic
      num_true<-sum((ys_pred_true-mean(ys))^2)
      den_true<-sum((ys-ys_pred_true)^2)
      F_true<-(s-2)*(num_true/den_true)
      
      #naive F-statistic
      num_naive<-sum((ys_pred_naive-mean(ys))^2)
      den_naive<-sum((ys-ys_pred_naive)^2)
      F_naive<-(s-2)*(num_naive/den_naive)
      
      
      
      dec<-c(F_true,F_naive)
      dec
      
      
    }
  
   
   F1<-function(u)
   {
     p<-mean(ifelse(B[,1]<u,1,0))
     return(p)
   }
   
   F2<-function(u)
   {
     p<-mean(ifelse(B[,2]<u,1,0))
     return(p)
   }
   
   max_diff<-max(abs(F1(u_vec)-F2(u_vec)))
   
  
  
    
  new_row<-c(n,l,max_diff,nrow(B))
  print(new_row)
  RP<-rbind(RP,new_row)
}

stopCluster(clust)

RP<-RP[-1,]
df<-data.frame(RP)
save(df,file="isoFstat.RData")
max_diff_vec<-RP[,3]

ggplot(data = df, aes(x=RP[,1],y=RP[,3],
                      colour="red")) +
  geom_point()+
  geom_line() +
  xlab(TeX("no. of nodes (n) \\rightarrow"))+
  ylab(TeX(" $sup_{x \\in R}|\\hat{F}_K(x)- \\hat{F}(x)|$
           "))




