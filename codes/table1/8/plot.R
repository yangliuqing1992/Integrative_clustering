naout_mean=function(x){
  x=x[x!=0]
  if(sum(is.na(x))==0)return(mean(x))
  else return(mean(x[(!is.na(x))]))
}

a=read.csv('1.csv');a=a[,-1]
b=read.csv('2.csv');b=b[,-1]
c=read.csv('3.csv');c=c[,-1]
d=read.csv('4.csv');d=d[,-1]
#e=read.csv('5.csv');e=e[,-1]
# f=read.csv('6.csv');f=f[,-1]
# g=read.csv('7.csv');g=g[,-1]
# h=read.csv('8.csv');h=h[,-1]
# i=read.csv('9.csv');i=i[,-1]

par(mfrow=c(2,2),mar=c(4, 2.6, 1.8, 1.5))
#t=0.05

avg=apply(a,1,naout_mean)
ARI.A=avg[1];ARI.IG=avg[2:10]; ARI.X=avg[11:19]; ARI.AX=avg[20:28]
p.noise.A=seq(0.1,0.9,length=9)
plot(p.noise.A,ARI.IG,ylim=c(0,0.5),type='b',col='red',
     ylab='Adjusted Rand', xlab='Proportion of noise in A, p.noise.X=0.1')
lines(p.noise.A,ARI.X,type='b')
lines(p.noise.A, ARI.AX, type='b',col='green')
abline(h=ARI.A,col='blue')
#abline(h=t,col='gold')


avg=apply(b,1,naout_mean)
ARI.A=avg[1];ARI.IG=avg[2:10]; ARI.X=avg[11:19]; ARI.AX=avg[20:28]
p.noise.A=seq(0.1,0.9,length=9)
plot(p.noise.A,ARI.IG,ylim=c(0,0.5),type='b',col='red',
     ylab='Adjusted Rand', xlab='Proportion of noise in A, p.noise.X=0.2')
lines(p.noise.A,ARI.X,type='b')
lines(p.noise.A, ARI.AX, type='b',col='green')
abline(h=ARI.A,col='blue')
#abline(h=t,col='gold')


avg=apply(c,1,naout_mean)
ARI.A=avg[1];ARI.IG=avg[2:10]; ARI.X=avg[11:19]; ARI.AX=avg[20:28]
p.noise.A=seq(0.1,0.9,length=9)
plot(p.noise.A,ARI.IG,ylim=c(0,0.5),type='b',col='red',
     ylab='Adjusted Rand', xlab='Proportion of noise in A, p.noise.X=0.3')
lines(p.noise.A,ARI.X,type='b')
lines(p.noise.A, ARI.AX, type='b',col='green')
abline(h=ARI.A,col='blue')
#abline(h=t,col='gold')

avg=apply(d,1,naout_mean)
ARI.A=avg[1];ARI.IG=avg[2:10]; ARI.X=avg[11:19]; ARI.AX=avg[20:28]
p.noise.A=seq(0.1,0.9,length=9)
plot(p.noise.A,ARI.IG,ylim=c(0,0.5),type='b',col='red',
     ylab='Adjusted Rand', xlab='Proportion of noise in A, p.noise.X=0.4')
lines(p.noise.A,ARI.X,type='b')
lines(p.noise.A, ARI.AX, type='b',col='green')
abline(h=ARI.A,col='blue')
#abline(h=t,col='gold')


legend('topright',cex=0.6,c('IG','SBM','NHM','NHM vs SBM'),lty=c(1,1,1,1),
       col=c('red','black','blue','green'),pch=c(1,1,0,1))

avg=apply(e,1,naout_mean)
ARI.A=avg[1];ARI.IG=avg[2:10]; ARI.X=avg[11:19]; ARI.AX=avg[20:28]
p.noise.A=seq(0.1,0.9,length=9)
plot(p.noise.A,ARI.IG,ylim=c(0,1),type='b',col='red',
     ylab='Adjusted Rand', xlab='Proportion of noise in A, p.noise.X=0.5')
lines(p.noise.A,ARI.X,type='b')
lines(p.noise.A, ARI.AX, type='b',col='green')
abline(h=ARI.A,col='blue')

legend('topright',cex=0.6,c('IG','SBM','NHM','NHM vs SBM'),lty=c(1,1,1,1),
       col=c('red','black','blue','green'),pch=c(1,1,0,1))