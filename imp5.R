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
library(princurve)

e <- new.env()
e$libs <- c("Matrix","MASS","irlba","princurve",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))



sig_ep<-0.1
d<-4
m_vec<-c(0.5,0.5,0.5,0.5)
n_vec<-seq(600,2500,100)


beta<-5.0
alpha<-2.0


clusterExport(clust,list("d","m_vec","alpha","beta"))

RS<-matrix(,ncol=8)


for(n in n_vec)
{
  
  clusterExport(clust,"n")
  
  L<-foreach(i=1:100,.combine = 'rbind') %dopar%
    {
      #generating regressors
      t<-runif(n,min=0,max=1)
      
      #generating response variable values
      e<-rnorm(n,mean=0,sd=sig_ep)
      y<-alpha+beta*t+e
      
      
      #forming the matrix of latent positions
      X<-cbind(t/2,t/2,t/2,t/2)
      
      #forming probability matrix
      P<-X%*%t(X)
      
      #generating adjacency matrix
      pvec<-c(P[lower.tri(P,diag=TRUE)])
      avec<-rbinom(length(pvec),1,pvec)
      A<-matrix(nrow=nrow(P),ncol=ncol(P))
      A[lower.tri(A,diag=TRUE)]<-avec
      A[upper.tri(A)]<-t(A)[upper.tri(A)]
      
      
      
      #finding adjacency spectral embedding
      A_irlba<-irlba(A,d)
      X_hat_raw<-A_irlba$u%*%diag(A_irlba$d)^0.5
      
      #rotating ASE to find consistent estimates of latent positions
      vals<-svd(t(X_hat_raw)%*%X)
      W<-vals$u%*%t(vals$v)
      X_hat<-X_hat_raw%*%W
      
      
      #estimating regressors from projections
      t_hat<-X_hat%*%m_vec
      
      t_ct<-t-mean(t)
      t_hat_ct<-t_hat-mean(t_hat)
      
      reg_par<-c(alpha,beta)
      
      
      #true estimator
      beta_true<-(t(t_ct)%*%y)/(t(t_ct)%*%t_ct)
      alpha_true<-mean(y)-beta_true*mean(t)
      reg_par_true<-c(alpha_true,beta_true)
      loss1<-(norm(reg_par_true-reg_par,type = "2"))^2
      loss1a<-(alpha_true-alpha)^2
      loss1b<-(beta_true-beta)^2
      
      
      
      #naive estimator
      beta_naive<-(t(t_hat_ct)%*%y)/(t(t_hat_ct)%*%t_hat_ct)
      alpha_naive<-mean(y)-beta_naive*mean(t_hat)
      reg_par_naive<-c(alpha_naive,beta_naive)
      loss2<-(norm(reg_par_naive-reg_par,type = "2"))^2
      loss2a<-(alpha_naive-alpha)^2
      loss2b<-(beta_naive-beta)^2
      
      
      
      
      
      dec<-c(loss1a,loss1b,loss2a,loss2b,loss1,loss2)
      dec
      
    }
  
  risk_all<-apply(L,2,mean)
  new_row<-c(n,risk_all,nrow(L))
  RS<-rbind(RS,new_row)
  
  
  print(new_row)
  
  
  
}

stopCluster(clust)


RS<-RS[-1,]

df<-data.frame(RS)
save(df,file="P1new5.RData")
load("P1new5.RData")

risk1a_vec<-RS[,2]
risk1b_vec<-RS[,3]
risk2a_vec<-RS[,4]
risk2b_vec<-RS[,5]
risk1_vec<-RS[,6]
risk2_vec<-RS[,7]


library(ggplot2)
library(reshape2)
library(latex2exp)

df<-data.frame(n_vec,risk1_vec,risk2_vec)
dfm<-melt(df, id.vars = 'n_vec')
print(dfm)
ggplot(dfm, aes(x=n_vec, y=value, 
                colour = variable)) +
  geom_point() +
  geom_line() +
  ylab(TeX("sample MSEs $\\rightarrow$")) +
  xlab(TeX("number of nodes(n) $\\rightarrow$")) +
  #ggtitle("Consistency of regression parameter estimates on
   #       known linear manifold") +
  scale_colour_manual(values = c("red","orange"),
                      labels=unname(TeX(c(
                        "sample MSE of  $\\hat{\\theta}_{true}$",
                        "sample MSE of $\\hat{\\theta}_{naive}$"))))
