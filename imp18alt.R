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
library(MASS)
library(igraph)
library(smacof)


e <- new.env()
e$libs <- c("irlba","Matrix","vegan3d","smacof",
            "MASS","igraph",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))




d<-4

s<-5

n_vec<-seq(1000,10000,1000)



alpha<-2.0
beta<-5.0
sig_ep<-0.01

lsg<-0.01
thres<-qf(lsg,1,s-2,lower.tail = FALSE)

u_vec<-seq(0,100,0.00001)

clusterExport(clust,list("d","s","alpha","beta","sig_ep"))

RP<-matrix(,ncol=3)

power_diff_vec<-vector()

for(n in n_vec)
{
  
  lambda<-0.6
  l<-as.integer(n/10)
  
  
  clusterExport(clust,list("n","lambda"))
  
  
  #opts<-list(preschedule=FALSE)
  
  B<-foreach(trial=1:50,.combine='c') %dopar%
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
      X_hat<-A_irlba$u%*%diag(A_irlba$d)^0.5
      
      #finding consistent adjacency spectral estimates of latent positions
      #vals<-svd(t(X_hat_raw)%*%X)
      #W<-vals$u%*%t(vals$v)
      #X_hat<-X_hat_raw%*%W
      
      
      B0<-as.matrix(dist(X_hat[1:l,],method = "euclidean",
              diag=TRUE,upper=TRUE))
      
      BB<-ifelse(B0<lambda,B0,0)
      
      
      colnames(BB)<-as.character(seq(1,nrow(B0),1))
      
      
      
      g<-graph_from_adjacency_matrix(BB,
                                      mode="undirected",
                                      weighted=TRUE,
                                      diag=FALSE,
                                      add.colnames = NULL,
                                      add.rownames = NA)
      
     
      
      
      
      #matrix of shortest path distances
      D<-shortest.paths(g, v=V(g),to=V(g))
      D<-as.matrix(D)
      Ds<-D[1:s,1:s]
      
      
      
      MM<-mds(Ds,ndim = 1,type = "interval",
          weightmat = NULL,
          init = "torgerson")
      
      
      zs<-as.vector(MM$conf)
      
      Dts<-dist(ts,method = "euclidean",
                diag = TRUE,
                upper = TRUE)
      
      Dzs<-dist(zs,method = "euclidean",
                diag = TRUE,
                upper = TRUE)
      
      
      dec<-max(abs(Dts-Dzs))
      #print(dec)
      dec
      
      
    }
  
  
  
  
  
  #power_true<-mean(ifelse(B[,1]<thres,1,0))
  #power_naive<-mean(ifelse(B[,2]<thres,1,0))
  
  new_row<-c(n,length(B),mean(B))
  print(new_row)
  
  

  RP<-rbind(RP,new_row)
}

stopCluster(clust)

RP<-RP[-1,]
df<-data.frame(RP)
save(df,file="P1new18alt.RData")
