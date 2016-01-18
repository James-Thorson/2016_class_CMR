library(coda)
library(BRugs)

setwd("E:\\Alain\\Consultancy\\_2008\\Seal-houl-out\\R\\Alain2\\")

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

#----------- Burn in demonstration --------------------

modelCheck("Modelglm-25.5.txt")
modelData("Sealmatrix.txt")
modelData("data_param.txt")

modelCompile(numChains = 3)
modelInits("InitializeParam1.txt");
modelInits("InitializeParam2.txt"); 
modelInits("InitializeParam3.txt"); 
modelInits("Sealnitmatrix.txt")
modelInits("Sealnitmatrix.txt")
modelInits("Sealnitmatrix.txt")

samplesSet("alpha")
samplesSet("b")
samplesSet("W")
samplesSet("S")

# samle each 50-th simulation; total is the same 200x50=10,000, but picture is better
#modelUpdate(200,thin=50)
modelUpdate(200,thin=1)

# plot it for in B/W for the book, mention that mixing can be better seen in color
plotHistory("alpha",colour=c("black","black","black"))

# plot it for in colour, to see mixing
plotHistory("alpha")

#we can plot all to see (select one for the book)
plotHistory("b[1]")
plotHistory("b[2]")
plotHistory("b[3]")
plotHistory("b[4]")
plotHistory("b[5]")
plotHistory("S[2]")
plotHistory("W[2]")
plotHistory("W[3]")
plotHistory("W[4]")
 
#---------------- start from scratch -------------------------

modelCheck("Modelglm-25.5.txt")
modelData("Sealmatrix.txt")
modelData("data_param.txt")

modelCompile(numChains = 3)
modelInits("InitializeParam1.txt");
modelInits("InitializeParam2.txt"); 
modelInits("InitializeParam3.txt"); 
modelInits("Sealnitmatrix.txt")
modelInits("Sealnitmatrix.txt")
modelInits("Sealnitmatrix.txt")

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

#Monitor model parameters
modelUpdate(10000,thin=10)

#------------- Convergence checks ----------------------

tmp.mcmc<-buildMCMC("alpha");  a.diag1<-gelman.diag(tmp.mcmc); a.diag2<-geweke.diag(tmp.mcmc)
a.diag1; rbind(a.diag2[[1]]$z,a.diag2[[2]]$z,a.diag2[[3]]$z)

tmp.mcmc<-buildMCMC("b");  b.diag1<-gelman.diag(tmp.mcmc); b.diag2<-geweke.diag(tmp.mcmc)
b.diag1; rbind(b.diag2[[1]]$z,b.diag2[[2]]$z,b.diag2[[3]]$z)

tmp.mcmc<-buildMCMC("W");  w.diag1<-gelman.diag(tmp.mcmc); w.diag2<-geweke.diag(tmp.mcmc)
w.diag1; rbind(w.diag2[[1]]$z,w.diag2[[2]]$z,w.diag2[[3]]$z)

tmp.mcmc<-buildMCMC("S");  s.diag1<-gelman.diag(tmp.mcmc); s.diag2<-geweke.diag(tmp.mcmc)
s.diag1; rbind(s.diag2[[1]]$z,s.diag2[[2]]$z,s.diag2[[3]]$z)


#------------- Sampled Data Statistics -----------------

samplesStats("alpha")
samplesStats("b")
samplesStats("W")
samplesStats("S")
samplesStats("SS")
samplesStats("SS.prd")
samplesStats("BP.s1")
samplesStats("BP.s2")
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
write.table(round(tmp.cor,2),"tmp_cor.dat",sep="\t")



