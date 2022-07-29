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
e$libs <- c("irlba","Matrix","princurve",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))



sig_ep<-0.1
d<-4
m_vec<-c(0.5,0.5,0.5,0.5)
n_vec<-seq(500,1500,50)


beta<-5.0


clusterExport(clust,list("d","m_vec","beta","sig_ep"))

RS<-matrix(,ncol=5)


for(n in n_vec)
{
  
  clusterExport(clust,"n")
  
  L<-foreach(i=1:100,.combine = 'rbind') %dopar%
    {
      #generating regressors
      t<-runif(n,min=0,max=1)
      
      #generating response variable values
      e<-rnorm(n,mean=0,sd=sig_ep)
      y<-beta*t+e
      
      
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
      
      
      
      delta<-matrix(0,d,d)
      
      
      for(j in 1:n)
      {
        delta<-delta+(X_hat[j,]%*%t(X_hat[j,]))
      }
      Delta<-(1.0/n)*delta
      
      Delinv<-ginv(Delta,tol = sqrt(.Machine$double.eps))
      
      incr_vec<-vector()
      for(i in 1:n)
      {
        sum_mat<-matrix(0,d,d)
        
        for(j in 1:n)
        {
          
          F1<-(t(X_hat[i,])%*%X_hat[j,])*(1-t(X_hat[i,])%*%X_hat[j,])
          f1<-as.double(F1)
          num<-f1*(X_hat[j,]%*%t(X_hat[j,]))
          sum_mat<-sum_mat+num
        }
        
        Sum_mat<-(1.0/n)*sum_mat
        
        incr<-t(m_vec)%*%Delinv%*%Sum_mat%*%Delinv%*%m_vec
        
        incr_vec<-c(incr_vec,incr)
        
      }
      
      sum_Gamma_hat<-mean(incr_vec)
      
      
      
      #true estimator
      beta_true<-(t(t)%*%y)/(t(t)%*%t)
      loss1<-(beta_true-beta)^2
      
      
      
      #naive estimator
      beta_naive<-(t(t_hat)%*%y)/(t(t_hat)%*%t_hat)
      loss2<-(beta_naive-beta)^2
      
      
      #ME-adjusted estimator with estimated unknown variance
      beta_adj_par_hat<-(t(t_hat)%*%y)/(t(t_hat)%*%t_hat-sum_Gamma_hat)
      loss4<-(beta_adj_par_hat-beta)^2
      
      
      dec<-c(loss1,loss2,loss4)
      dec
      
    }
  
  risk_all<-apply(L,2,mean)
  new_row<-c(n,risk_all,nrow(L))
  RS<-rbind(RS,new_row)
  
  
  print(new_row)
  
  
  
}

stopCluster(clust)

df<-data.frame(RS)
save(df,file = "P1new2.RData")
load("P1new2.RData")


RS<-RS[-1,]
risk1_vec<-RS[,2]
risk2_vec<-RS[,3]
risk4_vec<-RS[,4]


library(ggplot2)
library(reshape2)
library(latex2exp)

df<-data.frame(n_vec,risk1_vec,risk2_vec,risk4_vec)
dfm<-melt(df, id.vars = 'n_vec')
print(dfm)
ggplot(dfm, aes(x=n_vec, y=value, 
                colour = variable)) +
  geom_point() +
  geom_line() +
  ylab(TeX("sample MSEs $\\rightarrow$")) +
  xlab(TeX("number of nodes(n) $\\rightarrow$")) +
  #ggtitle("Consistency of regression parameter estimates on
    #      known linear manifold") +
  scale_colour_manual(values = c("red","orange","blue"),
                      labels=unname(TeX(c(
                        "sample MSE of  $\\hat{\\beta}_{true}$",
                        "sample MSE of $\\hat{\\beta}_{naive}$",
                        "sample MSE of $\\hat{\\beta}_{adj,\\hat{\\sigma}}$"
                      )))) +
  theme(axis.title.x = element_text(size = 15, vjust = -.2))

