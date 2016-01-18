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



library(AED); data(RiceFieldBirds)
RFBirds<-RiceFieldBirds
RFBirds$Richness<-rowSums(RFBirds[,8:56] > 0)
RFBirds$fField <-factor(RFBirds$FIELD)
library(lattice)
xyplot(Richness~Time|fField, data= RFBirds,
    panel=function(x,y){
      panel.grid(h=-1, v= 2)
      panel.points(x,y,col=1)
      panel.loess(x,y,col=1,lwd=2)})
#Or
xyplot(Richness ~ Time | fField, data= RFBirds,
       type = c("p", "smooth", "grid"),
       span=0.5,lwd=2, col=1, col.line=1)


RFBirds$LA<-log(RFBirds$AREA)
RFBirds$fSptreat <- factor(RFBirds$SPTREAT)
RFBirds$DEPTH2 <- RFBirds$DEPTH^2
M0 <- glm(Richness~offset(LA)+fSptreat+DEPTH+
             DEPTH2,family=quasipoisson,
             data=RFBirds)
summary(M0)



library(geepack)
M.gee1<-geeglm(Richness~offset(LA)+fSptreat+DEPTH+DEPTH2,
            data = RFBirds,
            family=poisson,id=fField,corstr = "ar1")
summary(M.gee1)


M.gee2<-geeglm(Richness~offset(LA)+ fSptreat,
         data = RFBirds,
         family=poisson,id=FIELD,corstr = "ar1")
anova(M.gee1,M.gee2)





#Owls
library(AED) ; data(Owls)
library(nlme)
Owls$NCalls<-Owls$SiblingNegotiation
Owls$LBroodSize<-log(Owls$BroodSize)

Form<-formula(NCalls~offset(LBroodSize)+
              SexParent*FoodTreatment+
              SexParent*ArrivalTime)
O1<-glm(Form,family=quasipoisson,data=Owls)
drop1(O1,test="F")
O2<-update(O1,.~. -SexParent*ArrivalTime)
drop1(O2,test="F")

Form<-formula(NCalls~offset(LBroodSize)+
              FoodTreatment+
              ArrivalTime)
O3<-glm(Form,family=quasipoisson,data=Owls)
drop1(O3,test="F")

#GEE:
library(geepack)
Form<-formula(NCalls~offset(LBroodSize)+
              SexParent*FoodTreatment+
              SexParent*ArrivalTime)
O2<-geeglm(Form,data=Owls,
            family=poisson,id=Nest,corstr = "exchangeable")


N<-length(Owls$Nest)
NewLevels<-c(paste(unique(Owls$Nest),".Dep",sep=""),
             paste(unique(Owls$Nest),".Sat",sep=""))
Owls$NestNight<-factor(levels=NewLevels)

for (i in 1:N){
  if (Owls$FoodTreatment[i]=="Deprived"){Owls$NestNight[i] <- paste(Owls$Nest[i],".Dep",sep="")}
  if (Owls$FoodTreatment[i]=="Satiated"){Owls$NestNight[i] <- paste(Owls$Nest[i],".Sat",sep="")}
}
cbind(Owls$Nest,Owls$FoodTreatment,Owls$NestNight)[1:5,]

cbind(Owls$Nest,Owls$FoodTreatment,Owls$NestNight)[1:5,]

#Alternative to the for loop with two if-statements:
Owls$NestNight<-factor(ifelse(Owls$FoodTreatment == "Deprived",
       paste(Owls$Nest, ".Dep",sep=""),
       paste(Owls$Nest, ".Sat",sep="")))



O3<-geeglm(Form,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")
summary(O3)



O3.A<-geeglm(NCalls~offset(LBroodSize)+
              SexParent+FoodTreatment+
              SexParent*ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")
            

O3.B<-geeglm(NCalls~offset(LBroodSize)+
              SexParent*FoodTreatment+
              SexParent+ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")

anova(O3,O3.A)
anova(O3,O3.B)


O4<-geeglm(NCalls~offset(LBroodSize)+
              SexParent+FoodTreatment+
              SexParent*ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")

O4.A<-geeglm(NCalls~offset(LBroodSize)+
              SexParent+FoodTreatment+
              SexParent+ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")


O4.B<-geeglm(NCalls~offset(LBroodSize)+
              SexParent+
              SexParent*ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")


anova(O4,O4.A)
anova(O4,O4.B)


            
O5<-geeglm(NCalls~offset(LBroodSize)+
              SexParent+FoodTreatment+ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")

O5.A<-geeglm(NCalls~offset(LBroodSize)+
              FoodTreatment+ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")

O5.B<-geeglm(NCalls~offset(LBroodSize)+
              SexParent+ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")

O5.C<-geeglm(NCalls~offset(LBroodSize)+
              SexParent+ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")

anova(O5,O5.A)
anova(O5,O5.B)
anova(O5,O5.C)

O6<-geeglm(NCalls~offset(LBroodSize)+
              FoodTreatment+ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")


O6.A<-geeglm(NCalls~offset(LBroodSize)+
              ArrivalTime,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")


O6.B<-geeglm(NCalls~offset(LBroodSize)+
              FoodTreatment,data=Owls,
            family=poisson,id=NestNight,corstr = "ar1")


anova(O6,O6.A)
anova(O6,O6.B)
summary(O6)






  #unstructured = alphs_ij


####################################


           
geese2<-function (formula = formula(data), sformula = ~1, id, waves = NULL,
    data = parent.frame(), subset = NULL, na.action = na.omit,
    contrasts = NULL, weights = NULL, zcor = NULL, corp = NULL,
    control = geese.control(...), b = NULL, alpha = NULL, gm = NULL,
    family = gaussian(), mean.link = NULL, variance = NULL, cor.link = "identity",
    sca.link = "identity", link.same = TRUE, scale.fix = TRUE,
    scale.value = 1, corstr = "independence", ...)
{
    scall <- match.call()
    mnames <- c("", "formula", "data", "offset", "weights", "subset",
        "na.action", "id", "waves", "corp")
    cnames <- names(scall)
    cnames <- cnames[match(mnames, cnames, 0)]
    mcall <- scall[cnames]
    if (is.null(mcall$id))
        mcall$id <- as.name("id")
    mcall[[1]] <- as.name("model.frame")
    m <- eval(mcall, parent.frame())
    y <- model.extract(m, response)
    if (is.null(dim(y)))
        N <- length(y)
    else N <- dim(y)[1]
    mterms <- attr(m, "terms")
    x <- model.matrix(mterms, m, contrasts)
    offset <- model.extract(m, offset)
    if (is.null(offset))
        offset <- rep(0, N)
    w <- model.extract(m, weights)
    if (is.null(w))
        w <- rep(1, N)
    id <- model.extract(m, id)
    waves <- model.extract(m, waves)
    corp <- model.extract(m, corp)
    if (is.null(id))
        stop("id variable not found.")
    mcall$formula <- formula
    mcall$formula[3] <- switch(match(length(sformula), c(0, 2,
        3)), 1, sformula[2], sformula[3])
    m <- eval(mcall, parent.frame())
    terms <- attr(m, "terms")
    zsca <- model.matrix(terms, m, contrasts)
    soffset <- model.extract(m, offset)
    if (is.null(soffset))
        soffset <- rep(0, N)
    if (is.character(family))
        family <- get(family)
    if (is.function(family))
        family <- family()
    ans <- geese.fit(x, y, id, offset, soffset, w, waves, zsca,
        zcor, corp, control, b, alpha, gm, family, mean.link,
        variance, cor.link, sca.link, link.same, scale.fix, scale.value,
        corstr, ...)
    ans <- c(ans, list(call = scall, formula = formula))
    class(ans) <- "geese"
    ans
}









#########################

#Farm	Sex	Length	Ecervi
library(AED); data(DeerEcervi)

DeerEcervi$Ecervi.01 <- DeerEcervi$Ecervi
DeerEcervi$Ecervi.01[DeerEcervi$Ecervi>0]<-1
DeerEcervi$fSex <- factor(DeerEcervi$Sex)
DeerEcervi$CLength <- DeerEcervi$Length - mean(DeerEcervi$Length)
#DeerEcervi$CLength <- scale(DeerEcervi$Length,scale=FALSE)

DE.glm<-glm(Ecervi.01 ~ CLength * fSex, data = DeerEcervi,
         family = binomial)
drop1(DE.glm, test = "Chi")
summary(DE.glm)


library(geepack)
DE.gee<-geeglm(Ecervi.01 ~ CLength * fSex, data = DeerEcervi,
           family=binomial,
           id=Farm,corstr = "exchangeable")
summary(DE.gee)


A3<-geese2(Vas1~LCT1*factor(Sex1),
           family=binomial,
           id=Farm1,corstr = "exchangeable")
summary(A3)

