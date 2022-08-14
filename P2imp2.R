library(Matrix)
library(MASS)
library(irlba)
library(vegan3d)

S<-5
sigma_ep<-0.01
D<-4
d<-1
alpha<-2.0
beta<-5.0

T<-200




N<-15
n<-2000
ind_vec<-vector()
qt_vec<-vector()
pred_diff_vec<-vector()

MM<-matrix(,ncol=2)

for(trial in 1:T) 
{
  #generating regressors
  tS<-runif(S,min=0.5,max=1)
  
  #generating auxiliary indices
  tM<-runif(N-S,min=0.5,max=1)
  
  #concatenating all indices
  t<-c(tS,tM)
  
  #randomly selecting set of iid outcomes
  omega<-runif(n,min=0.5,max=1)
  
  
  
  #matrix whose t-th column is ASE estimate for index=t
  X_hat_mat<-matrix(nrow=n,ncol=N)
  
  #creating X_hat_mat by concatenating ASE estimates by columns
  for(r in 1:length(t))
  {
    
    #latent position matrix
    X<-omega*t[r]
    
    
    #probability matrix
    P<-X%*%t(X)
    
    
    
    
    #generating adjacency matrix
    pvec<-c(P[lower.tri(P,diag=TRUE)])
    avec<-rbinom(length(pvec),1,pvec)
    A<-matrix(nrow=nrow(P),ncol=ncol(P))
    A[lower.tri(A,diag=TRUE)]<-avec
    A[upper.tri(A)]<-t(A)[upper.tri(A)]
    
    
    
    
    #finding adjacency spectral embedding
    A_irlba<-irlba(A,d)
    X_hat_raw<-A_irlba$u*(as.double(A_irlba$d)^0.5)
    
    
    #rotating to find consistent ASE estimates
    vals<-svd(t(X_hat_raw)%*%X)
    W<-vals$u%*%t(vals$v)
    X_hat<-X_hat_raw%*%W
    
    
    
    X_hat_mat[,r]<-X_hat
  }
  
  
  #approximate dissimilarity matrix upon which CMDS will run 
  D_2_phi_hat<-matrix(nrow=N,ncol=N)
  
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      mm<-(n^(-0.5))*(min(
        norm(X_hat_mat[,i]-X_hat_mat[,j],type="2"),
        norm(X_hat_mat[,i]+X_hat_mat[,j],type="2")
      ))
      D_2_phi_hat[i,j]<-mm^2
    }
  }
  
  
  #adjustment matrix for CMDS
  H<-matrix(nrow=N,ncol=N)
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      H[i,j]<-ifelse(i==j,1-1.0/N,-1.0/N)
    }
  }
  
  
  #Applying CMDS to find the embeddings
  F<-(-0.5)*(H%*%D_2_phi_hat%*%H)
  F_irlba<-irlba(F,D)
  psi_hat<-F_irlba$u%*%diag(F_irlba$d)^0.5
  
  
  
  z_hat<-isomap(vegdist(psi_hat, method="euclidean"),
                k=7,ndim=1,
                path="shortest")$points
  
  
  
  

  
  #generating responses
  yS<-alpha+beta*tS+rnorm(S,mean=0,sd=sigma_ep)
  
  #estimates based on true regressors
  beta_true<-cov(yS,tS)/var(tS)
  alpha_true<-mean(yS)-beta_true*mean(tS)
  
  
  #naive estimates
  beta_naive<-cov(yS,z_hat[1:S])/var(z_hat[1:S])
  alpha_naive<-mean(yS)-beta_naive*mean(z_hat[1:S])
  
  
  #predicted response based on true regressors
  t0<-t[S+1]
  y0<-alpha+beta*t0+rnorm(1,mean=0,sd=sigma_ep)
  y0_true<-alpha_true+beta_true*t0
  y0_naive<-alpha_naive+beta_naive*z_hat[S+1]
  
  pred_diff<-abs(y0_true-y0_naive)
  pred_diff_vec<-c(pred_diff_vec,pred_diff)
  print(c(trial,pred_diff))
  
  
  
}

print(c(n,mean(pred_diff_vec)))




