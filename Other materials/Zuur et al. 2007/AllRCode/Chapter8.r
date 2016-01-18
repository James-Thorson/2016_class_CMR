#    Mixed Effects Models and Extensions in Ecology with R (2009)
#    Zuur, Ieno, Walker, Saveliev and Smith.    Springer
#    This file was produced by Alain Zuur (highstat@highstat.com)
#    www.highstat.com

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.



library(AED); data(Sparrows)

op <- par(mfrow=c(2,2))
hist(Sparrows$wt,nclass=15,main="Observed data",xlab="Weight")
hist(Sparrows$wt,nclass=15,main="Observed
     data",xlab="Weight",freq=F)
Y<-rnorm(1281,mean=mean(Sparrows$wt),sd=sd(Sparrows$wt))
hist(Y,nclass=15,main="Simulated data",xlab="Weight")
X<-seq(0,30,length=200)
Y<-dnorm(X,mean=mean(Sparrows$wt),sd=sd(Sparrows$wt))
plot(X,Y,type="l",xlab="Weight",ylab="Probablities",
     ylim=c(0,0.25), xlim=c(0,30),
     main="Normal density curve")
par(op)



x1 <-0:10; Y1<-dpois(x1,lambda=3)
x2<-0:10; Y2<-dpois(x2,lambda=5)
x3<-0:40; Y3<-dpois(x3,lambda=10)
x4<-50:150; Y4<-dpois(x4,lambda=100)
XLab<-"Y values";YLab<-"Probabilities"
par(mfrow=c(2,2))
plot(x1,Y1,type="h",xlab=XLab,ylab=YLab,
     main="Poisson with mean 3")
plot(x2,Y2,type="h", xlab=XLab,ylab=YLab,
     main="Poisson with mean 5")
plot(x3,Y3,type="h", xlab=XLab,ylab=YLab,
     main="Poisson with mean 10")
plot(x4,Y4,type="h", xlab=XLab,ylab=YLab,
     main="Poisson with mean 100")
     
     
     
     
     
     
#Figure 8.3    Negative binomial
library(stats)
## Alternative parametrization

 mu1B=1 ; k1B=0.1
 mu1C=1 ; k1C=1
 mu1D=1 ; k1D=100000

 mu2B=10 ; k2B=0.1
 mu2C=10 ; k2C=1
 mu2D=10 ; k2D=100000

 mu3B=100 ; k3B=0.1
 mu3C=100 ; k3C=1
 mu3D=100 ; k3D=100000


x1B<-0:10; Y12<-dnbinom(x1B,mu=mu1B,size=k1B)
x1C<-0:10; Y13<-dnbinom(x1C,mu=mu1C,size=k1C)
x1D<-0:10; Y14<-dnbinom(x1D,mu=mu1D,size=k1D)

x2B<-0:20; Y22<-dnbinom(x2B,mu=mu2B,size=k2B)
x2C<-0:20; Y23<-dnbinom(x2C,mu=mu2C,size=k2C)
x2D<-0:20; Y24<-dnbinom(x2D,mu=mu2D,size=k2D)


x3B<-0:200; Y32<-dnbinom(x3B,mu=mu3B,size=k3B)
x3C<-0:200; Y33<-dnbinom(x3C,mu=mu3C,size=k3C)
x3D<-0:200; Y34<-dnbinom(x3D,mu=mu3D,size=k3D)


par(mfrow=c(3,3))
Xlab="Y values"
Ylab="Probabilities"
plot(x1B,Y12,type="h",xlab=Xlab,ylab=Ylab,main=paste("NB(",mu1B,",",k1B,")"))
plot(x1C,Y13,type="h",xlab=Xlab,ylab=Ylab,main=paste("NB(",mu1C,",",k1C,")"))
plot(x1D,Y14,type="h",xlab=Xlab,ylab=Ylab,main=paste("NB(",mu1D,",",k1D,")"))


plot(x2B,Y22,type="h",xlab=Xlab,ylab=Ylab,main=paste("NB(",mu2B,",",k2B,")"))
plot(x2C,Y23,type="h",xlab=Xlab,ylab=Ylab,main=paste("NB(",mu2C,",",k2C,")"))
plot(x2D,Y24,type="h",xlab=Xlab,ylab=Ylab,main=paste("NB(",mu2D,",",k2D,")"))


plot(x3B,Y32,type="h",xlab=Xlab,ylab=Ylab,main=paste("NB(",mu3B,",",k3B,")"))
plot(x3C,Y33,type="h",xlab=Xlab,ylab=Ylab,main=paste("NB(",mu3C,",",k3C,")"))
plot(x3D,Y34,type="h",xlab=Xlab,ylab=Ylab,main=paste("NB(",mu3D,",",k3D,")"))




#########################################################
#Binomial
#Figure 8.5

Xlab="Y values"
Ylab="Probabilities"

n11<-20; x11<-1:n11; p11<-0.2
n12<-20; x12<-1:n12; p12<-0.5
n13<-20; x13<-1:n13; p13<-0.7

n21<-10; x21<-1:n21; p21<-0.2
n22<-10; x22<-1:n22; p22<-0.5
n23<-10; x23<-1:n23; p23<-0.7


n31<-100; x31<-1:n31; p31<-0.2
n32<-100; x32<-1:n32; p32<-0.5
n33<-100; x33<-1:n33; p33<-0.7


prop11<-dbinom(x11, size=n11, prob=p11)
prop12<-dbinom(x12, size=n12, prob=p12)
prop13<-dbinom(x13, size=n13, prob=p13)


prop21<-dbinom(x21, size=n21, prob=p21)
prop22<-dbinom(x22, size=n22, prob=p22)
prop23<-dbinom(x23, size=n23, prob=p23)


prop31<-dbinom(x31, size=n31, prob=p31)
prop32<-dbinom(x32, size=n32, prob=p32)
prop33<-dbinom(x33, size=n33, prob=p33)


par(mfrow=c(3,3))


plot(x21,prop21,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(",p21,",",n21,")"))
plot(x22,prop22,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(",p22,",",n22,")"))
plot(x23,prop23,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(",p23,",",n23,")"))

plot(x11,prop11,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(",p11,",",n11,")"))
plot(x12,prop12,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(",p12,",",n12,")"))
plot(x13,prop13,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(",p13,",",n13,")"))


plot(x31,prop31,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(",p31,",",n31,")"))
plot(x32,prop32,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(",p32,",",n32,")"))
plot(x33,prop33,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(",p33,",",n33,")"))




#Figure 8.6
#Bernouli

par(mfrow=c(2,2))
prop11<-dbinom(0:1, size=1, prob=0.2)
plot(0:1,prop11,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(0.2,1)"),ylim=c(0,1))


prop11<-dbinom(0:1, size=1, prob=0.5)
plot(0:1,prop11,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(0.5,1)"),ylim=c(0,1))


prop11<-dbinom(0:1, size=1, prob=0.7)
plot(0:1,prop11,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(0.7,1)"),ylim=c(0,1))

prop11<-dbinom(0:1, size=1, prob=1)
plot(0:1,prop11,type="h",xlab=Xlab,ylab=Ylab,main=paste("B(1,1)"),ylim=c(0,1))




