rm(list=ls(all=TRUE))

library(kernlab)
library(mvtnorm)
calRand=function(A)
{
  n=sum(A)
  nTrue=apply(A,1,sum)
  nEst=apply(A,2,sum)
  pairA=choose(A,2)
  pairTrue=choose(nTrue,2)
  pairEst=choose(nEst,2)
  adjRand=(sum(pairA)-(sum(pairTrue)*sum(pairEst))/choose(n,2))/(1/2*(sum(pairTrue)+sum(pairEst))-(sum(pairTrue)*sum(pairEst))/choose(n,2))
  a=sum(pairA)
  b=sum(pairTrue)-sum(pairA)
  c=sum(pairEst)-sum(pairA)
  d=choose(n,2)-a-b-c
  Rand=(a+d)/choose(n,2)
  return(list(Rand=Rand, adjRand=adjRand))
}

matchLabel=function(trueLabel,myLabel)
{ 
  mTrue=max(trueLabel)
  mEst=max(myLabel)
  A=matrix(0,mTrue,mEst)
  for (i in 1:mTrue)
    for (j in 1:mEst)
    {
      A[i,j]=sum((trueLabel==i)&(myLabel==j))
      
    }
  return(calRand(A))
  
}


EM <- function(sample,p_inits,para_inits,maxit=5000,tol=10^-4)
{
  # Estimation of parameter(Initial)
  flag <- 0
  p_cur <- p_inits; para_cur <- para_inits
  # Iterate between expectation and maximization parts
  
  for(i in 1:maxit){
    new <- maximization(sample,expectation(sample,p_cur,para_cur),para_cur)
    p_new <- new[[1]]; para_new=new[[2]]
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    diff=max(max(abs(p_cur-p_new)),max(abs(unlist(para_new)-unlist(para_cur))))

    if( diff < tol ){ flag <- 1; break}
    
    # Otherwise continue iteration
    p_cur <- p_new; para_cur <- para_new
    
  }
  if(!flag) warning("Didn't converge\n")
  e=apply(expectation(sample,p_new,para_new),1,which.max)
  return(list(e=e, Pi=p_new,para=para_new))
}


expectation <- function(sample,p,para)
{ 
  K=length(p)
  p_e=matrix(0,length(sample[,1]),K)
  for (i in 1:length(sample[,1])){
    t=sapply(para,NN,x=sample[i,],w=w.2[i,])
    p_e[i,]=t*p/sum(t*p)
  }
  return(p_e)
}

maximization <- function(sample,epart,para){
  # estimate p

  p_temp <- apply(epart,2,sum)/sum(epart)
  para_temp=list()
  w=w.2
  mu=apply(w*sample,1,sum)/apply(w,1,sum)
  gm=(apply(w,1,sum)-1)/apply(w*(sample-mu)^2,1,sum)
  for(i in 1:length(epart[1,])){
    a=para[[i]][[1]];b=para[[i]][[2]]
    roots=optim(c(a,b),model.IG,df=sample,epart=epart[,i],w=w.2,mu=mu,gm=gm,gr=NULL,
                method = "BFGS")$par
    para_temp=append(para_temp,list(roots))
  }
  return(list(p_temp,para_temp))   
}



BIC_JCGS=function(A,e)
  
{
  K=max(e)
  n=dim(A)[1];
  L_A=0;
  M=matrix(0,K,K);
  N=matrix(0,K,K);
  
  for (k in 1:K)
    for (l in 1:K)
    {  M[k,l]=sum(A[e==k,e==l])
    if (k==l) N[k,l]=sum(e==k)*(sum(e==k)-1)
    if (k!=l) N[k,l]=sum(e==k)*sum(e==l)
    if ((M[k,l]>0) & (M[k,l]<N[k,l]))
      L_A=L_A+M[k,l]*log(M[k,l]/N[k,l])+(N[k,l]-M[k,l])*log(1-M[k,l]/N[k,l])
    
    } 
  logLik1=L_A/2;
  
  bic=logLik1-K*(K+1)/2*log(n)
  return(bic)
  
}

CBIC=function(A,e,lambda)
  
{
  K=max(e)
  n=dim(A)[1];
  
  L_A=0;
  M=matrix(0,K,K);
  N=matrix(0,K,K);
  
  for (k in 1:K)
    for (l in 1:K)
    {  M[k,l]=sum(A[e==k,e==l])
    if (k==l) N[k,l]=sum(e==k)*(sum(e==k)-1)
    if (k!=l) N[k,l]=sum(e==k)*sum(e==l)
    if ((M[k,l]>0) & (M[k,l]<N[k,l]))
      L_A=L_A+M[k,l]*log(M[k,l]/N[k,l])+(N[k,l]-M[k,l])*log(1-M[k,l]/N[k,l])
    
    } 
  logLik1=L_A/2;
  
  bic=logLik1-(lambda*n*log(K)+K*(K+1)/2*log(n))
  return(bic)
}

BIC_X=function(X,pi,para){
  p=dim(X)[1]; K=length(pi);p_e=matrix(0,p,K)
  for (i in 1:p){
    t=sapply(para,NN,x=X[i,],w=w.2[i,])
    p_e[i,]=t*pi
  }
  L=sum(log(apply(p_e,1,sum)))
  return(2*L-(3*K-1)*log(p))
}


label=function(x){
  return(which(x==1))
}


int.x=function(df,K){
  df.input=t(apply(cbind(df,w.2),1,function(x){y=x[1:n];z=x[(n+1):(2*n)];
  y[which(y==0)]=sum(y*z)/sum(z);return(y)}))
  cl=kmeans(df.input,K)$cluster; p_int=rep(0,K); para_int=list()
  for (i in 1:K){p_int[i]=sum(cl==i)/length(cl)}
  for (i in 1:K){
    w=w.2[cl==i,]
    df1=df[cl==i,]
    mu=apply(df1*w,1,sum)/apply(w,1,sum)
    gm=(apply(w,1,sum)-1)/apply(w*(df1-mu)^2,1,sum)
    a=mean(mu); b=1/var(mu)
    para=optim(c(a,b),model,df=df1,gm=gm,
               mu=mu,w=w,gr=NULL,method = "BFGS")$par
    para_int=append(para_int,list(para))
  }
  return(list(cl=cl,para_int=para_int, p_int=p_int))
}


model=function(x,df,gm,mu,w){
  x1=x[1]; x2=x[2]
  -sum(log(x2/(gm*apply(w,1,sum)+x2)))+
    sum(gm/(gm*apply(w,1,sum)+x2)*(x2*apply(w*(df-x1)^2,1,sum)+gm*apply(w,1,sum)*
                                     apply(w*(df-mu)^2,1,sum)))
}


pseudo=function(A,e)
{
  n=length(e);
  K=max(e);
  T=10
  for (t in 1:T)
  {
    
    Pi=rep(0,K);
    for (k in 1:K)
      Pi[k]=sum(e==k)/n;
    
    lambda=matrix(0,K,K);   
    P=matrix(0,K,K);
    for (l in 1:(K))
      for (k in 1:(K))
      {    
        ind1=which(e==l)
        ind2=which(e==k)
        if ((length(ind1)>=2) && (length(ind2)>=2))
        {
          if (l!=k)
          {P[l,k]=sum(A[ind1,ind2])/(length(ind1)*length(ind2));}
          if (l==k)       
          {P[l,k]=sum(A[ind1,ind1])/(length(ind1)*(length(ind1)-1)); }
        }
      }
    for (l in 1:(K))  
      for (k in 1:(K))
      {
        lambda[l,k]=sum(e==k)*P[l,k]
      }
    
    b=matrix(0,n,K);
    for (i in 1:n)
      for (k in 1:(K))
        b[i,k]=sum(A[i,e==k]);
    diff=1;
    
    while (diff>1e-3)
    {
      for (l in 1:K)
        for (k in 1:K)
          if (lambda[l,k]==0)
            lambda[l,k]=1e-6
          
          z=matrix(0,n,K);
          for (i in 1:n)
          {
            z[i,l]=log(Pi[l])
            for (l in 1:K)
            {
              for (k in 1:(K))
              { z[i,l]=z[i,l]+(b[i,k]*log(lambda[l,k])-lambda[l,k]); }
            }

            
            z[i,]=z[i,]-max(z[i,])
            z[i,]=exp(z[i,]);
            z[i,]=z[i,]/sum(z[i,]);   
            
          }
          lambdaNew=matrix(0,K,K);
          for (l in 1:(K))
            for (k in 1:(K))
            {
              lambdaNew[l,k]=sum(b[,k]*z[,l])/sum(z[,l]);
              if (is.nan(lambdaNew)[l,k]) lambdaNew[l,k]=0
            }
          
          PiNew=apply(z,2,sum);
          PiNew=PiNew[1:K]/sum(PiNew[1:K]);
          diff=max(max(abs(Pi-PiNew)),max(abs(lambda-lambdaNew)))
          
          if (is.na(diff))
            return(list(error=1));
          if (abs(diff)>200)
          {return(list(error=1)) }
          
          Pi=PiNew;
          lambda=lambdaNew;
          
    }
    e=apply(z,1,which.max);

    
  }
  
  return(list(error=0,e=e, Pi=Pi,lambda=lambda,b=b,z=z));
}

IntegrativeEM=function(A,e,X,para_int)
{
  p=length(e);
  K=max(e);
  T=12
  for (t in 1:T)
  {
    
    Pi=rep(0,K);
    for (k in 1:K)
      Pi[k]=sum(e==k)/p;
    
    lambda=matrix(0,K,K);   
    P=matrix(0,K,K);
    for (l in 1:(K))
      for (k in 1:(K))
      {    
        ind1=which(e==l)
        ind2=which(e==k)
        
        if ((length(ind1)>=2) && (length(ind2)>=2))
        {
          if (l!=k)
          {P[l,k]=sum(A[ind1,ind2])/(length(ind1)*length(ind2));}
          if (l==k)       
          {P[l,k]=sum(A[ind1,ind1])/(length(ind1)*(length(ind1)-1)); }
        }
      }
    for (l in 1:(K))  
      for (k in 1:(K))
      {
        lambda[l,k]=sum(e==k)*P[l,k]
      }
    
    for (l in 1:K)
      for (k in 1:K)
        if (lambda[l,k]==0)
          lambda[l,k]=1e-6
        
        
        b=matrix(0,p,K);
        for (i in 1:p)
          for (k in 1:(K))
            b[i,k]=sum(A[i,e==k]);
        
        diff=1
        para_cur=para.mg(para_int,lambda);p_cur=Pi
        
        t=1
        while ((diff>1e-3)&&(t<1000))
        { 
          t=t+1
          
          z=expectation.IG(sample1=X,sample2=b,pi=p_cur,para=para_cur)
          
          ##################update Lambda###################
          lambdaNew=matrix(0,K,K);
          for (l in 1:(K))
            for (k in 1:(K))
            {
              lambdaNew[l,k]=sum(b[,k]*z[,l])/sum(z[,l]);
              if (is.nan(lambdaNew)[l,k]) lambdaNew[l,k]=0
              
            }
          
          for (l in 1:K)
            for (k in 1:K)
              if (lambdaNew[l,k]==0)
                lambdaNew[l,k]=1e-6
          
          ################update \mu_k, \gamma_k##############
          para_x_cur=para.cl(para_cur)[[1]]
          para_x_new=maximization.IG(sample=X,epart=z,para=para_x_cur)
          
          ###############merge to para_new###############
          para_new=para.mg(para_x_new,lambdaNew)
          
          ###############update \pi_k################
          p_new=apply(z,2,sum);
          p_new=p_new[1:K]/sum(p_new[1:K]);
          diff=max(max(abs(p_cur-p_new)),max(abs(unlist(para_new)-unlist(para_cur))))
          
          if (is.na(diff))
            return(list(error=1));
          if (abs(diff)>200)
          {return(list(error=1)) }
          
          para_cur=para_new;p_cur=p_new
          
          Pi=p_new;
          lambda=lambdaNew;
          para_x=para_x_new
          
        }
        e=apply(z,1,which.max);

        
  }
  
  return(list(error=0,para_x=para_x,e=e, Pi=Pi,P=P,lambda=lambda,b=b,z=z));
}

para.mg=function(para1,para2){
  K=length(para1);para=list()
  for (i in 1:K){
    para=append(para,list(list(para1[[i]],para2[i,])))
  }
  return(para)
}

expectation.IG <- function(sample1,sample2,pi,para)
{ 
  K=dim(sample2)[2]
  p_e=matrix(0,p,K)
  for (i in 1:p){
    t=sapply(para,IG,x=sample1[i,],w=w.2[i,],B=sample2[i,])
    t=t-max(t)
    t=exp(t)
    p_e[i,]=t*pi/sum(t*pi)
  }
  return(p_e)
}


IG=function(para,x,w,B){
  n=length(x)
  mu=sum(w*x)/sum(w); gm=(sum(w)-1)/sum(w*(x-mu)^2)

  
  log.pdf=-sum(para[[2]])+sum(B*log(para[[2]]))+log((2*pi)^(-sum(w)/2)*gm^(sum(w)/2)*sqrt(para[[1]][2]/(para[[1]][2]+gm*sum(w))))-0.5*(gm/(gm*sum(w)+para[[1]][2]))*(
    para[[1]][2]*sum(w*(x-para[[1]][1])^2)+sum(w)*gm*sum(w*(x-mu)^2))
  
  
  
  return(log.pdf)
}

para.cl=function(para){
  K=length(para);para1=list();para2=list()
  for (i in 1:K){
    para1=append(para1, list(para[[i]][[1]]))
    para2=append(para2,list(para[[i]][[2]]))
  }
  return(list(para1,para2))
}



maximization.IG <- function(sample,epart,para){
  para_temp=list()
  w=w.2
  mu=apply(w*sample,1,sum)/apply(w,1,sum)
  gm=(apply(w,1,sum)-1)/apply(w*(sample-mu)^2,1,sum)
  K=length(para)
  
  for(i in 1:K){
    a=para[[i]][[1]];b=para[[i]][[2]]
    roots=rep(0,2)
    roots=optim(c(a,b),model.IG,df=sample,epart=epart[,i],w=w.2,mu=mu,gm=gm,gr=NULL,
                method = "BFGS")$par
    para_temp=append(para_temp,list(roots))
  }
  return(para_temp)   
}

model.IG=function(x,df,epart,w,mu,gm){
  x1=x[1]; x2=x[2]
  -sum(epart*log(x2/(gm*apply(w,1,sum)+x2)))+
    sum(epart*gm/(gm*apply(w,1,sum)+x2)*(x2*apply(w*(df-x1)^2,1,sum)+gm*apply(w,1,sum)*
                                           apply(w*(df-mu)^2,1,sum)))
}

NN=function(para,x,w){
  n=length(x)
  mu=sum(w*x)/sum(w); gm=(sum(w)-1)/sum(w*(x-mu)^2)
  pdf=(2*pi)^(-sum(w)/2)*gm^(sum(w)/2)*sqrt(para[[2]]/(para[[2]]+
                                                         gm*sum(w)))*exp(-0.5*(gm/(gm*sum(w)+para[[2]]))*(
                                                           para[[2]]*sum(w*(x-para[[1]])^2)+sum(w)*gm*sum(w*(x-mu)^2)))
  return(pdf)
}





###Generate simulated data


p.noise.x=0.2
p.noise.A=0.2
K=4

p=1000; n=40; 
m=rmultinom(p, size = 1, prob = rep(1/K, K))

l=apply(m,2,label)
para1=list(-0.5,0.6)
para2=list(2,0.6)
para3=list(4,0.8)
para4=list(1, 0.6)

para=list(para1, para2, para3, para4)

P=matrix(rep(0.06, K^2), K, K); diag(P)=0.12
p_0=0.0001
u=rep(0,p)

for (i in 1:p){
  u[i]=rnorm(1, para[[l[i]]][[1]], para[[l[i]]][[2]])
}
se=runif(p,0.2,1.5)
df=matrix(0,p,n)

for(i in 1:p){
  df[i,]=rnorm(n,u[i],se[i])
}

df.noise=df
for (i in 1:p){
  if(rbinom(1,1,p.noise.x)==1){df.noise[i,]=rnorm(n,2.5,1)}
}

A.noise=matrix(0,p,p); 
for (i in 1:(p-1)){
  for (j in (i+1):p){
    if (rbinom(1,1,1-p.noise.A)==1){
      if (rbinom(1,1,P[l[i],l[j]])==1){A.noise[i,j]=1}
    }
    else {
      if (rbinom(1,1,p_0)==1) {A.noise[i,j]=1}
    }
  }
}
A.noise=A.noise+t(A.noise)


### Main part

K_hat=4


e=tryCatch(specc(as.kernelMatrix(A.noise),K_hat),error=function(e){return(NaN)})
if (is.nan(e[1]))
  e=specc(A.noise,K_hat)

r.A=pseudo(A.noise,e)

w.2=matrix(1,p,n)
x_initials = int.x(df.noise,K_hat)
r.X=EM(df.noise,x_initials$p_int,x_initials$para_int,maxit=5000,tol=10^-3)


r.AX=IntegrativeEM(A.noise,e,df.noise,x_initials$para_int)

BIC_value=BIC_JCGS(A.noise,r.AX$e)
CBIC_value=CBIC(A.noise, r.AX$e, lambda=0.5)



ARI.X=matchLabel(l,r.X$e)$adjRand
ARI.A=matchLabel(l,r.A$e)$adjRand
ARI.IG=matchLabel(l,r.AX$e)$adjRand
ARI.AX=matchLabel(r.X$e, r.A$e)$adjRand

### r.X result of NHM
### r.A result of SBM
### r.AX result of integrative method

### ARI.X ARI for NHM
### ARI.A ARI for SBM
### ARI.IG ARI for integrative method
### ARI.AX ARI for comparing NHM and SBM 


