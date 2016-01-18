library(coda)
library(BRugs)

setwd("D:\\SAA\\Alain\\Seal-houl-out\\R\\Alain\\")
setwd("D:\\SAA\\Alain\\R\\Alain\\")
setwd("E:\\Alain\\Consultancy\\_2008\\Seal-houl-out\\R\\Alain\\")
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
#Section 23.6

#Prepare the data file
# ------------- data matrix ------------------
fn.data<-"Sealmatrix.txt"
p.names<-c("Abun","Month1","Month2","TDay1","TDay2","Wind2","Site")
Datamatrix<-cbind(Seals$Abun,Seals$Month1,Seals$Month2,Seals$TDay1,Seals$TDay2,Seals$fWind2,Seals$fSite)
pncol<-ncol(Datamatrix)
sbuf<-""
for(i in 1:pncol) {
  sbuf<-sprintf("%s%s[] ",sbuf,p.names[i])
}
write(sbuf,fn.data,append=F)
write.table(Datamatrix,fn.data,row.names=F,col.names=F,append=T)
write("END",fn.data,append=T)

# ------------- prediction init matrix ------------------
fn.data<-"Sealnitmatrix.txt"
p.names<-c("Aprd")
Datamatrix<-cbind(Seals$Abun)
pncol<-ncol(Datamatrix)
sbuf<-""
for(i in 1:pncol) {
  sbuf<-sprintf("%s%s[] ",sbuf,p.names[i])
}
write(sbuf,fn.data,append=F)
write.table(Datamatrix,fn.data,row.names=F,col.names=F,append=T)
write("END",fn.data,append=T)

fn.data<-"Sealnitmatrix2.txt"
p.names<-c("Aprd","eps")
Datamatrix<-cbind(Seals$Abun,rep(0,98))
pncol<-ncol(Datamatrix)
sbuf<-""
for(i in 1:pncol) {
  sbuf<-sprintf("%s%s[] ",sbuf,p.names[i])
}
write(sbuf,fn.data,append=F)
write.table(Datamatrix,fn.data,row.names=F,col.names=F,append=T)
write("END",fn.data,append=T)

modelCheck("Modelglm-25.6.txt")
modelData("Sealmatrix.txt")
modelData("data_param.txt")

modelCompile(numChains = 3)
modelInits("InitializeParam1.txt"); 
modelInits("InitializeParam2.txt"); 
modelInits("InitializeParam3.txt"); 
modelInits("Sealnitmatrix2.txt"); 
modelInits("Sealnitmatrix2.txt"); 
modelInits("Sealnitmatrix2.txt"); 
modelInits("InitializeTau.txt")
modelInits("InitializeTau.txt")
modelInits("InitializeTau.txt")


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
samplesSet("eps")
samplesSet("eps.tau")
samplesSet("eps.sigma")

#Monitor model parameters
modelUpdate(10000,thin=10)
samplesStats("alpha")
samplesStats("b")
samplesStats("W")
samplesStats("S")
samplesStats("SS")
samplesStats("SS.prd")
samplesStats("BP.s1")
samplesStats("BP.s2")
samplesStats("eps.tau")
samplesStats("eps.sigma")
dicStats()

eps.tau  <-as.double(plotHistory("eps.tau",plot=F))
eps.sigma<-sqrt(1/eps.tau)
sigma.stat<-c(mean(eps.sigma),quantile(eps.sigma,c(0.05,0.5,0.975)))
round(sigma.stat,3)

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
write.table(round(tmp.cor,2),"tmp_cor2.dat",sep="\t")

p.eps<-data.frame(b1=as.double(plotHistory("eps[1]",plot=F)))
for(i in 2:98) {
  tmp.names<-names(p.eps)
  sbuf<-sprintf("eps[%d]",i)
  p.eps<-data.frame(p.eps,as.double(plotHistory(sbuf,plot=F)))
  names(p.eps)<-c(tmp.names,sprintf("eps%d",i))
}
eps.mean<-apply(p.eps,2,mean)
eps.sd<-apply(p.eps,2,sd)
hist(eps.mean)
hist(eps.sd)
mean(eps.sd^2)

p.dist<-matrix(ncol=98,nrow=98)
for(i in 1:98) {
  for(j in 1:98) {
    p.dist[i,j]<-abs(Seals$WeekTime[i]-Seals$WeekTime[j])
  }
}
p.eps.cor<-cor(p.eps)
is.na(diag(p.eps.cor))<-T
hist(p.eps.cor)

plot(as.double(p.dist),as.double(p.eps.cor),main="Correlation against WeekTime difference",
     xlab="WeekTime difference",ylab="Correlation",pch=20,cex=0.3) 

plot(as.double(p.dist[1:51,1:51]),as.double(p.eps.cor[1:51,1:51]),main="Correlation against WeekTime difference",
     xlab="WeekTime difference",ylab="Correlation",pch=20,cex=0.3) 

plot(as.double(p.dist[52:98,52:98]),as.double(p.eps.cor[52:98,52:98]),main="Correlation against WeekTime difference",
     xlab="WeekTime difference",ylab="Correlation",pch=20,cex=0.3) 
