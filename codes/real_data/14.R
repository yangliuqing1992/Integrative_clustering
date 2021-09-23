library(stats)
library(fcd)
label=function(x){
  return(which(x==1))
}
zerout_mean=function(x){
  if(length(which(x==0))==0)return(mean(x))
  else return(mean(x[-(which(x==0))]))
}
zerout_var=function(x){
  if(length(which(x==0))==0)return(var(x))
  else return(var(x[-(which(x==0))]))
}

zero.omit=function(x){return(x[x!=0])}

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
        if (l!=k)
        {P[l,k]=sum(A[ind1,ind2])/(length(ind1)*length(ind2));}
        if (l==k)       
        {P[l,k]=sum(A[ind1,ind1])/(length(ind1)*(length(ind1)-1)); }
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
    
    while (diff>1e-4)
    {
      z=matrix(0,n,K);
      for (i in 1:n)
      {
        
        # for (l in 1:K)
        # {
        #   z[i,l]=Pi[l]
        #   for (k in 1:(K))
        #     #z[i,l]=z[i,l]*exp(-lambda[l,k]+log(lambda[l,k])*b[i,k]);
        #     z[i,l]=z[i,l]*exp(b[i,k]*log(lambda[l,k])-lambda[l,k]-log.fact(b[i,k]));
        #   if(is.na(z[i,l])) z[i,l]=0;
        # }
        # 
        # z[i,]=z[i,]/sum(z[i,]);   
        temp=rep(0,K)
        for (l in 1:K){
          temp[l]=log(Pi[l])
          for (k in 1:K){
            if(lambda[l,k]==0){lambda[l,k]=0.00001}
            temp[l]=temp[l]+b[i,k]*log(lambda[l,k])-lambda[l,k]
          }
        }
        max1=max(temp)
        z[i,]=exp(temp-max1)
        z[i,]=z[i,]/sum(z[i,])
        
      }
      lambdaNew=matrix(0,K,K);
      for (l in 1:(K))
        for (k in 1:(K))
        {
          lambdaNew[l,k]=sum(b[,k]*z[,l])/sum(z[,l]);
          
        }
      
      PiNew=apply(z,2,sum);
      PiNew=PiNew[1:K]/sum(PiNew[1:K]);
      diff=max(max(abs(Pi-PiNew)),max(abs(lambda-lambdaNew)))
      
      #cat("diff=",diff,"\n")
      if (is.na(diff))
        return(list(error=1));
      if (abs(diff)>100)
      {return(list(error=1)) }
      
      #print(cat("diff=",diff,"\n"));
      Pi=PiNew;
      lambda=lambdaNew;
      
    }
    e=apply(z,1,which.max);
    #cat(e,"\n")
    
    
  }
  
  return(list(error=0,e=e, Pi=Pi,lambda=lambda,b=b,z=z));
}

NN=function(para,x,w){
  n=length(x)
  #w=weights1(x); 
  mu=sum(w*x)/sum(w); gm=(sum(w)-1)/sum(w*(x-mu)^2)
  pdf=(2*pi)^(-sum(w)/2)*gm^(sum(w)/2)*sqrt(para[[2]]/(para[[2]]+
                                                         gm*sum(w)))*exp(-0.5*(gm/(gm*sum(w)+para[[2]]))*(
                                                           para[[2]]*sum(w*(x-para[[1]])^2)+sum(w)*gm*sum(w*(x-mu)^2)))
  return(pdf)
}


#Expectation
expectation <- function(sample,p,para)
{ 
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
  #w=1/(1-probs)
  #w=t(apply(sample,1,weights))
  w=w.2
  mu=apply(w*sample,1,sum)/apply(w,1,sum)
  gm=(apply(w,1,sum)-1)/apply(w*(sample-mu)^2,1,sum)
  #gm=1/(se)^2
  for(i in 1:length(epart[1,])){
    a=para[[i]][[1]];b=para[[i]][[2]]
    # roots=optim(c(a,b),model,df=sample,epart=epart[,i],w=w,mu=mu,gm=gm,gr=NULL,
    #             method = "L-BFGS-B",lower=c(-100,0))$par
    roots=optim(c(a,b),model.IG,df=sample,epart=epart[,i],w=w.2,mu=mu,gm=gm,gr=NULL,
                method = "BFGS")$par
    para_temp=append(para_temp,list(roots))
  }
  return(list(p_temp,para_temp))   
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
    #cat("diff=",diff,"\n")
    # if( all(abs(p_cur - p_new) < tol)& all(abs(para_cur[[1]] - para_new[[1]]) < tol)
    #     &all(abs(para_cur[[2]] - para_new[[2]]) < tol)&
    #     all(abs(para_cur[[3]] - para_new[[3]]) < tol)){ flag <- 1; break}
    if( diff < tol ){ flag <- 1; break}
    
    # Otherwise continue iteration
    p_cur <- p_new; para_cur <- para_new
    
  }
  if(!flag) warning("Didn't converge\n")
  e=apply(expectation(sample,p_new,para_new),1,which.max)
  return(list(e=e, Pi=p_new,para=para_new))
  #return(list(p_cur,para_cur))
}


weights2=function(x,data){
  w=rep(0,length(x))
  y=x[x!=0]; min=min(y); max=max(y); r=max-min
  if(sum(x==0)==0) w=rep(1,length(x))
  else{
    for (j in 1:15){
      data1=matrix(0,1,length(x))
      for (i in 1:length(data[,1])){
        t=zero.omit(data[i,])
        if (min(t)>=(min-j/10*r)&max(t)<=(max+j/10*r)){
          data1=rbind(data1,data[i,])
        }
      }
      if(length(data1[,1])>=11) break
    }
    data1=data1[-1,]
    props=apply(data1,1,function(x){sum(x==0)/length(x)})
    means=apply(data1,1,zerout_mean)
    lm=lm(props~means); w[x!=0]=1/(1-predict(lm,data.frame(means=x[x!=0])))
  }
  return(w)
}

log.fact=function(x){
  sapply(x,function(x){
    if(x==0) {return(0)} 
    else{n=x;r=0
    for(i in 1:n){
      r=r+log(i)
    }
    }
    return(r)})
}

# IG=function(para,x,w,B){
#   n=length(x)
#   mu=sum(w*x)/sum(w); gm=(sum(w)-1)/sum(w*(x-mu)^2)
#   # pdf=exp(-sum(para[[2]])+sum(B*log(para[[2]]))-sum(log.fact(B)))*
#   #   (2*pi)^(-sum(w)/2)*gm^(sum(w)/2)*sqrt(para[[1]][2]/(para[[1]][2]+
#   #            gm*sum(w)))*exp(-0.5*(gm/(gm*sum(w)+para[[1]][2]))*(
#   #               para[[1]][2]*sum(w*(x-para[[1]][1])^2)+sum(w)*gm*sum(w*(x-mu)^2)))
#   pdf1=exp(-sum(para[[2]])+sum(B*log(para[[2]]))-sum(log.fact(B)))
#   if (is.na(pdf1)){pdf1=0}
#   pdf2=(2*pi)^(-sum(w)/2)*gm^(sum(w)/2)*sqrt(para[[1]][2]/(para[[1]][2]+
#            gm*sum(w)))*exp(-0.5*(gm/(gm*sum(w)+para[[1]][2]))*(
#               para[[1]][2]*sum(w*(x-para[[1]][1])^2)+sum(w)*gm*sum(w*(x-mu)^2)))
#   #print(pdf1); print(pdf2)
#   pdf=pdf1*pdf2
#   return(pdf)
# }
# 
# 
# expectation.IG <- function(sample1,sample2,pi,para)
# { 
#   p_e=matrix(0,p,K)
#   for (i in 1:p){
#     t=sapply(para,IG,x=sample1[i,],w=w.2[i,],B=sample2[i,])
#     #print(t)
#     p_e[i,]=t*pi/sum(t*pi)
#   }
#   return(p_e)
# }

expectation.IG <- function(sample1,sample2,pi,para){ 
  K=length(pi)
  p_e=matrix(0,p,K)
  for (i in 1:p){
    temp=rep(0,K); pdf2=rep(0,K)
    w=w.2[i,]; x=sample1[i,]
    mu=sum(w*x)/sum(w); gm=(sum(w)-1)/sum(w*(x-mu)^2)
    for (l in 1:K){
      pdf2[l]=(2*pi)^(-sum(w)/2)*gm^(sum(w)/2)*sqrt(para[[l]][[1]][2]/(para[[l]][[1]][2]+
                                                                         gm*sum(w)))*exp(-0.5*(gm/(gm*sum(w)+para[[l]][[1]][2]))*(
                                                                           para[[l]][[1]][2]*sum(w*(x-para[[l]][[1]][1])^2)+sum(w)*gm*sum(w*(x-mu)^2)))
      temp[l]=log(pi[l])
      for (k in 1:K){
        if (para[[l]][[2]][k]==0){para[[l]][[2]][k]=0.00001}
        #if(i<3){print(c(para[[l]][[2]][k],sample2[i,k]))}
        temp[l]=temp[l]+sample2[i,k]*log(para[[l]][[2]][k])-para[[l]][[2]][k]
        #if(i<3){print(temp[l])}
      }
    }
    
    temp=temp-max(temp); pdf1=exp(temp)
    #print(pdf1); #print(pdf2)
    t=pdf1*pdf2
    p_e[i,]=t*pi/sum(t*pi)
  }
  return(p_e)
}

model.IG=function(x,df,epart,w,mu,gm){
  x1=x[1]; x2=x[2]
  -sum(epart*log(x2/(gm*apply(w,1,sum)+x2)))+
    sum(epart*gm/(gm*apply(w,1,sum)+x2)*(x2*apply(w*(df-x1)^2,1,sum)+gm*apply(w,1,sum)*
                                           apply(w*(df-mu)^2,1,sum)))
}



maximization.IG <- function(sample,epart,para){
  para_temp=list()
  w=w.2
  mu=apply(w*sample,1,sum)/apply(w,1,sum)
  gm=(apply(w,1,sum)-1)/apply(w*(sample-mu)^2,1,sum)
  for(i in 1:K){
    a=para[[i]][[1]];b=para[[i]][[2]]
    roots=rep(0,2)
    roots=optim(c(a,b),model.IG,df=sample,epart=epart[,i],w=w.2,mu=mu,gm=gm,gr=NULL,
                method = "BFGS")$par
    para_temp=append(para_temp,list(roots))
  }
  return(para_temp)   
}


para.mg=function(para1,para2){
  K=length(para1);para=list()
  for (i in 1:K){
    para=append(para,list(list(para1[[i]],para2[i,])))
  }
  return(para)
}

para.cl=function(para){
  K=length(para);para1=list();para2=list()
  for (i in 1:K){
    para1=append(para1, list(para[[i]][[1]]))
    para2=append(para2,list(para[[i]][[2]]))
  }
  return(list(para1,para2))
}


IntegrativeEM=function(A,e,X,para_int)
{
  p=length(e);
  K=max(e)  
  T=12
  for (t in 1:T)
  {
    if(length(unique(e))!=K){print("Hi");e=clean_label(e); K=max(e);para_int=int.IG.x(X,e)$para_int}
    #print(K)
    Pi=rep(0,K);
    for (k in 1:K)
      Pi[k]=sum(e==k)/p;
    #print(Pi)
    lambda=matrix(0,K,K);   
    P=matrix(0,K,K);
    for (l in 1:(K))
      for (k in 1:(K))
      {    
        ind1=which(e==l)
        ind2=which(e==k)
        if (l!=k)
        {P[l,k]=sum(A[ind1,ind2])/(length(ind1)*length(ind2));}
        if (l==k)       
        {
          if (length(ind1)==1){P[l,k]=0}
          else {P[l,k]=sum(A[ind1,ind1])/(length(ind1)*(length(ind1)-1))} }
      }
    for (l in 1:(K))  
      for (k in 1:(K))
      {
        lambda[l,k]=sum(e==k)*P[l,k]
      }
    #print(lambda)
    b=matrix(0,p,K);
    for (i in 1:p)
      for (k in 1:(K))
        b[i,k]=sum(A[i,e==k]);
    
    diff=1
    para_cur=para.mg(para_int,lambda);p_cur=Pi
    #print(para_cur)
    while (diff>1e-4)
    { 
      
      z=expectation.IG(sample1=X,sample2=b,pi=p_cur,para=para_cur)
      #print(z)
      ##################update Lambda###################
      lambdaNew=matrix(0,K,K);
      for (l in 1:(K))
        for (k in 1:(K))
        {
          lambdaNew[l,k]=sum(b[,k]*z[,l])/sum(z[,l]);
        }
      #print(lambdaNew)
      ################update \mu_k, \gamma_k##############
      para_x_cur=para.cl(para_cur)[[1]]
      para_x_new=maximization.IG(sample=X,epart=z,para=para_x_cur)
      #print(para_x_cur)
      ###############merge to para_new###############
      para_new=para.mg(para_x_new,lambdaNew)
      
      ###############update \pi_k################
      p_new=apply(z,2,sum);
      p_new=p_new[1:K]/sum(p_new[1:K]);
      diff=max(max(abs(p_cur-p_new)),max(abs(unlist(para_new)-unlist(para_cur))))
      
      #cat("diff=",diff,"\n"); print(abs(p_cur-p_new))
      # print(abs(unlist(para_new)-unlist(para_cur)))
      #print(para_cur);print(para_new)
      if (is.na(diff))
        return(list(error=1));
      if (abs(diff)>1000)
      {return(list(error=1)) }
      
      para_cur=para_new;p_cur=p_new
      
      #print(cat("diff=",diff,"\n"));
      Pi=p_new;
      lambda=lambdaNew;
      para_x=para_x_new
      
    }
    e=apply(z,1,which.max); #print(pk(e))
    #cat(e,"\n")
  }
  
  return(list(error=0,para_x=para_x,e=e, Pi=Pi,P=P,lambda=lambda,b=b,z=z));
}

model=function(x,df,gm,mu,w){
  x1=x[1]; x2=x[2]
  -sum(log(x2/(gm*apply(w,1,sum)+x2)))+
    sum(gm/(gm*apply(w,1,sum)+x2)*(x2*apply(w*(df-x1)^2,1,sum)+gm*apply(w,1,sum)*
                                     apply(w*(df-mu)^2,1,sum)))
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

int.IG.x=function(df,e){
  df.input=t(apply(cbind(df,w.2),1,function(x){y=x[1:n];z=x[(n+1):(2*n)];
  y[which(y==0)]=sum(y*z)/sum(z);return(y)}))
  K=max(e)
  p_int=rep(0,K); para_int=list()
  for (i in 1:K){p_int[i]=sum(e==i)/length(e)}
  for (i in 1:K){
    w=w.2[e==i,]
    df1=df[e==i,]
    mu=apply(df1*w,1,sum)/apply(w,1,sum)
    gm=(apply(w,1,sum)-1)/apply(w*(df1-mu)^2,1,sum)
    a=mean(mu); b=1/var(mu)
    para=optim(c(a,b),model,df=df1,gm=gm,
               mu=mu,w=w,gr=NULL,method = "BFGS")$par
    para_int=append(para_int,list(para))
  }
  return(list(para_int=para_int, p_int=p_int))
}


pk=function(e){
  K=max(e)
  return(sapply(1:K, function(x){return(sum(e==x))}))
}

BIC_FengYang=function(A,e)
  
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

clean_label=function(e){
  u=unique(e)
  return(sapply(e,function(x){return(which(u==x))}))
}



### Read real data sets

set.seed(57397)
df=read.csv("X.csv",row.names = 1)
n=dim(df)[2]; p=dim(df)[1]; K=14
w.2=t(apply(df,1,weights2,data=df))
initial=int.x(df,K)

net=read.csv("A.csv",header=F)
A=as.matrix(net)

Xaff=as.matrix(read.csv("Xaff.csv",header=F))


### Main part

mi=0
while(mi<2){
  e=spectral.clustering(A+Xaff,K=K)
  mi=min(pk(e))
}
IG=IntegrativeEM(A,e,df,initial$para_int)
l.pred=IG$e

mi=0
while(mi<2){
  e.A=spectral.clustering(A,K=K)
  mi=min(pk(e.A))
}
PS=pseudo(A,e.A)
l.pred.A=PS$e

HC=EM(df,initial$p_int,initial$para_int,maxit=5000,tol=10^-4)
l.pred.X=HC$e

## write labels for each method out

result=data.frame(rownames(df),e,l.pred,initial$cl,l.pred.X,e.A,l.pred.A)
names(result)=c("gene","ig_int","ig","km_x","hc","sc_a","pl_sc")
write.csv(result,paste("labels_",as.character(K),".csv",sep=""),row.names = F)


## Print out BIC and CBIC(lambda = 0.5, 2) values
cat("BIC_A_IGlabels, K=",K,"\n","BIC=",BIC_FengYang(A,l.pred),"\n","CBIC0.5=",CBIC(A, l.pred, lambda=0.5),
    "\n","CBIC2=",CBIC(A,l.pred,lambda=2),"\n")