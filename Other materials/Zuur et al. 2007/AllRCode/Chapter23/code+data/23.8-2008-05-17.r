library(coda)
library(BRugs)

setwd("D:\\SAA\\Alain\\Seal-houl-out\\R\\Alain\\")
setwd("E:\\Alain\\Consultancy\\_2008\\Seal-houl-out\\R\\Alain2\\")
setwd("D:/SAA/Alain/R/Alain")
Seals<-read.table("countsBK1.txt",header=T)

library(AED); data(Seals)
Seals$fSite <- factor(Seals$Site)

Seals$Time=Seals$Year+(Seals$Week-1)/52

library(lattice)
xyplot(Abun~Time|fSite, data=Seals,
 ylab="Abundance",xlab="Time (years)",
 panel=function(x,y){
   panel.loess(x,y,span=0.3,col=1)
   panel.xyplot(x,y,col=1)
 }
)

Seals$fSeason<-Seals$Month
I1 <- Seals$Month==1 | Seals$Month==2 | Seals$Month==12
I2 <- Seals$Month==3 | Seals$Month==4 | Seals$Month==5
I3 <- Seals$Month==6 | Seals$Month==7 | Seals$Month==8
I4 <- Seals$Month==9 | Seals$Month==10| Seals$Month==11

Seals$fSeason[I1] <-"a"
Seals$fSeason[I2]  <-"b"
Seals$fSeason[I3]  <-"c"
Seals$fSeason[I4] <-"d"
Seals$fSeason<-as.factor(Seals$fSeason)

Seals$fWind2<-Seals$Winddir
Seals$fWind2[Seals$Winddir==1 | Seals$Winddir==2]<-1
Seals$fWind2[Seals$Winddir==3 | Seals$Winddir==4]<-2
Seals$fWind2[Seals$Winddir==5 | Seals$Winddir==6]<-3
Seals$fWind2[Seals$Winddir==7 | Seals$Winddir==8]<-4
Seals$fWind2<-factor(Seals$fWind2)

fSeason2<-Seals$fSeason
library(mgcv)
Seals$TDay <- Seals$Timeofday
M1<-gamm(Abun ~ s(Month,TDay) + fWind2+
        fSite,
        weights=varIdent(form=~1|fSeason2),data=Seals)
        
plot(M1$gam,pers=TRUE)
anova(M1$gam)

Seals$Month1=as.double(scale(Seals$Month))
Seals$Month2=Seals$Month1^2
Seals$TDay1<-as.double(scale(Seals$TDay))
Seals$TDay2<-Seals$TDay1^2
Seals$TDay3<-Seals$TDay1^3

M2 <- glm(Abun~Month1+Month2+TDay1+TDay2+Month1:TDay1+
              fWind2+fSite, data = Seals,
              family=poisson)
summary(M2)
drop1(M2,test="Chi")

Seals$E2 <- resid(M2,type="pearson")
xyplot(E2~Time|fSite,data=Seals,
  ylab="Pearson residuals",xlab="Time (years)",
  panel=function(x,y){
    panel.loess(x,y,span=0.5,col=1)
    panel.xyplot(x,y,col=1)})

########  Did not make it into the book
#Show 3-dim picture of Month and Tday1   TDay3+
M2 <- glm(Abun~Month1+Month2+TDay1+TDay2+
               Month1:TDay1+
              fWind2+fSite, data = Seals,
              family=poisson)

a<-seq(min(Seals$Month1),max(Seals$Month1),length.out=100)
b<-seq(min(Seals$TDay1),max(Seals$TDay1),length.out=100)

P2<-matrix(nrow=100,ncol=100)
for (i in 1:100){
    MyData<-data.frame(Month1=a,
                   Month2=a^2,
                   TDay1=b[i],
                   TDay2=b[i]^2,
                   TDay3=b[i]^3,
                   fWind2="1",
                   fSite="1")
    P2[,i]<-predict(M2,newdata=MyData,type="response")
    }

wireframe(P2)

###################################
#Section 23.5


modelCheck("Modelglm-25.8.txt")
modelData("Sealmatrix3.txt")
modelData("data_param.txt")

modelCompile(numChains = 3)
modelInits("InitializeParam1.txt"); 
modelInits("InitializeParam2.txt"); 
modelInits("InitializeParam3.txt"); 
for(i in 1:3) modelInits("Sealnitmatrix.txt"); 
for(i in 1:3) modelInits("Sealnitmatrix4.txt"); 
for(i in 1:3) modelInits("InitializeScale.txt")


#Burn in
modelUpdate(200,thin=50)

dicSet()
samplesSet("alpha")
samplesSet("b")
samplesSet("W")
samplesSet("S")
samplesSet("SS")
samplesSet("SS.prd")
samplesSet("BP.s1")
samplesSet("BP.s2")
samplesSet("rho1")
samplesSet("rho2")
samplesSet("u")
samplesSet("u.tau")
samplesSet("u.sigma")
samplesSet("size")

#Monitor model parameters
modelUpdate(10000,thin=10)

tmp.mcmc<-buildMCMC("alpha"); a.diag1<-gelman.diag(tmp.mcmc); a.diag2<-geweke.diag(tmp.mcmc)
a.diag1; rbind(a.diag2[[1]]$z,a.diag2[[2]]$z,a.diag2[[3]]$z)

tmp.mcmc<-buildMCMC("b");  b.diag1<-gelman.diag(tmp.mcmc); b.diag2<-geweke.diag(tmp.mcmc)
b.diag1; rbind(b.diag2[[1]]$z,b.diag2[[2]]$z,b.diag2[[3]]$z)

tmp.mcmc<-buildMCMC("W");  w.diag1<-gelman.diag(tmp.mcmc); w.diag2<-geweke.diag(tmp.mcmc)
w.diag1; rbind(w.diag2[[1]]$z,w.diag2[[2]]$z,w.diag2[[3]]$z)

tmp.mcmc<-buildMCMC("S");  s.diag1<-gelman.diag(tmp.mcmc); s.diag2<-geweke.diag(tmp.mcmc)
s.diag1; rbind(s.diag2[[1]]$z,s.diag2[[2]]$z,s.diag2[[3]]$z)

tmp.mcmc<-buildMCMC("rho1"); r1.diag1<-gelman.diag(tmp.mcmc); r1.diag2<-geweke.diag(tmp.mcmc)
r1.diag1; rbind(r1.diag2[[1]]$z,r1.diag2[[2]]$z,r1.diag2[[3]]$z)

tmp.mcmc<-buildMCMC("rho2"); r2.diag1<-gelman.diag(tmp.mcmc); r2.diag2<-geweke.diag(tmp.mcmc)
r2.diag1; rbind(r2.diag2[[1]]$z,r2.diag2[[2]]$z,r2.diag2[[3]]$z)

tmp.mcmc<-buildMCMC("u.tau"); t.diag1<-gelman.diag(tmp.mcmc); t.diag2<-geweke.diag(tmp.mcmc)
t.diag1; rbind(t.diag2[[1]]$z,t.diag2[[2]]$z,t.diag2[[3]]$z)

tmp.mcmc<-buildMCMC("size"); z.diag1<-gelman.diag(tmp.mcmc); z.diag2<-geweke.diag(tmp.mcmc)
z.diag1; rbind(z.diag2[[1]]$z,z.diag2[[2]]$z,z.diag2[[3]]$z)

samplesStats("alpha")
samplesStats("b")
samplesStats("W")
samplesStats("S")
samplesStats("SS")
samplesStats("SS.prd")
samplesStats("BP.s1")
samplesStats("BP.s2")
samplesStats("u.sigma")
samplesStats("rho1")
samplesStats("rho2")
samplesStats("size")
dicStats()
 
p.param<-data.frame(alpha=as.double(plotHistory("alpha",plot=F)))
for(i in 1:5) {
  tmp.names<-names(p.param)
  sbuf<-sprintf("b[%d]",i)
  p.param<-data.frame(p.param,as.double(plotHistory(sbuf,plot=F)))
  names(p.param)<-c(tmp.names,sprintf("b%d",i))
}
tmp.names<-names(p.param)
p.param<-data.frame(p.param,as.double(plotHistory("S[2]",plot=F)))
names(p.param)<-c(tmp.names,"S2")
for(i in 2:4) {
  tmp.names<-names(p.param)
  sbuf<-sprintf("W[%d]",i)
  p.param<-data.frame(p.param,as.double(plotHistory(sbuf,plot=F)))
  names(p.param)<-c(tmp.names,sprintf("W%d",i))
}
tmp.cor<-cor(p.param)
round(tmp.cor,2)
write.table(round(tmp.cor,2),"tmp_cor3.dat",sep="\t")

p.eps<-data.frame(b1=as.double(plotHistory("u[1]",plot=F)))
for(i in 2:98) {
  tmp.names<-names(p.eps)
  sbuf<-sprintf("u[%d]",i)
  p.eps<-data.frame(p.eps,as.double(plotHistory(sbuf,plot=F)))
  names(p.eps)<-c(tmp.names,sprintf("eps%d",i))
}
u.mean<-apply(p.eps,2,mean)
u.sd<-apply(p.eps,2,sd)
hist(u.mean)
hist(u.sd)
mean(u.sd^2)

p.dist<-matrix(ncol=98,nrow=98)
for(i in 1:98) {
  for(j in 1:98) {
    p.dist[i,j]<-abs(Seals$WeekTime[i]-Seals$WeekTime[j])
  }
}
p.u.cor<-cor(p.eps)
is.na(diag(p.u.cor))<-T
hist(p.u.cor)

plot(as.double(p.dist),as.double(p.u.cor),main="Correlation against WeekTime difference",
     xlab="WeekTime difference",ylab="Correlation",pch=20,cex=0.3) 

plot(as.double(p.dist[1:51,1:51]),as.double(p.u.cor[1:51,1:51]),main="Correlation against WeekTime difference",
     xlab="WeekTime difference",ylab="Correlation",pch=20,cex=0.3,ylim=c(-0.7,0.7)) 

plot(as.double(p.dist[52:98,52:98]),as.double(p.u.cor[52:98,52:98]),main="Correlation against WeekTime difference",
     xlab="WeekTime difference",ylab="Correlation",pch=20,cex=0.3,ylim=c(-0.7,0.7)) 

--------------------------------------------------------------------
var.data<-data.frame(mu=seq(50,300,by=1))
var.data$Vnb<-var.data$mu+var.data$mu^2/13.77
plot(var.data$mu,var.data$Vnb,type="l",xlab="mu",ylab="variance",main="variance model")
lines(var.data$mu,var.data$mu,lty=2)
