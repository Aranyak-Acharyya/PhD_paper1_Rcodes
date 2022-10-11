library(parallel)
library(iterators)
no_of_cores<-detectCores()
clust<-parallel::makeCluster(no_of_cores)
library(foreach)
library(doParallel)
registerDoParallel()




e <- new.env()
e$libs <- c("irlba","Matrix","vegan3d",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))




load("parl_reg.RData")
n<-df$n_vec[length(df$n_vec)]+250



#n<-1000

d<-4

s<-5


ll<-n/5
l<-as.integer(ll)


KK<-l/5
K<-as.integer(KK)


alpha<-2.0
beta<-5.0
sig_ep<-0.01

clusterExport(clust,list("d","n","s","l","K",
                         "alpha","beta","sig_ep"))


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

stopCluster(clust)

print(nrow(B))


print(apply(B,2,mean))

#only from 2nd n onwards, meant to append RData file
load("parl_reg.RData")
n_vec<-c(df$n_vec,n)
pred_val_diff_mat<-rbind(df[,-1],apply(B,2,mean))
df<-data.frame(cbind(n_vec,pred_val_diff_mat))
save(df,file = "parl_reg.RData")
load("parl_reg.RData")
print(df$n_vec)
print(df[,-1])





load("parl_reg.RData")
plot(df$n_vec,df$X1,type="l",col="red")
plot(df$n_vec,df[,2],type="l",col="green")





#only for first n, meant to initiate RData file
n_vec<-n
pred_val_diff_row<-apply(B,2,mean)
df<-data.frame(n_vec,t(pred_val_diff_row))
save(df,file="parl_reg.RData")



library(tidyverse)
library(ggplot2)
library(reshape2)
library(latex2exp)
load("parl_reg.RData")
print(df)
#df<-data.frame(df$n_vec,df$X1)
#print(df)
#dfm<-melt(df, id.var='df.n_vec')
#print(dfm)
ggplot(data = df, aes(x=n_vec, y=X1, colour = "red")) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = "red",
                    labels=unname(TeX(c("$E(\\hat{y}_{true}-\\hat{y}_{naive})^2$"
                                      )))) +
  labs(x="no. of nodes (n)",y="mean sq predicted difference") +
  theme(legend.box.background = element_rect(color="red")) 
