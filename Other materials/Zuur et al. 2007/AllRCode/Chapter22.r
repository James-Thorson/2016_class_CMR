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




library(AED); data(BadgersFarmSurveys)
Badgers <- BadgersFarmSurveys

colSums(sapply(Badgers,FUN=is.na))




BadgersFarmSurveysNoNA=read.table(
    file="C:/Bookdata/BadgersFarmSurveysNoNA.txt",
    header=T)

library(AED); data(BadgersFarmSurveysNoNA)
Badgers <- BadgersFarmSurveysNoNA

Badgers$fSeason <- factor(Badgers$Season)
Badgers$fFeed.store  <- factor(Badgers$Accessible_feed_store_present)
Badgers$fCattle.house  <- factor(Badgers$Accessible_cattle_house_present)
Badgers$fFeed.present <- factor(Badgers$Accessible_feed_present)
Badgers$fGrass.silage <- factor(Badgers$Grass_silage)
Badgers$fCereal.silage <- factor(Badgers$Cereal_silage)
Badgers$fHayStraw     <- factor(Badgers$HayStraw)
Badgers$fCereal.grains <- factor(Badgers$Cereal_grains)
Badgers$fConcentrates  <- factor(Badgers$Concentrates)
Badgers$fSugarbeet     <- factor(Badgers$Sugarbeet)
Badgers$fMolasses      <- factor(Badgers$Molasses)

B.glm <- glm(Signs_in_yard~fSeason+
       No_setts_in_fields+No_buildings +
       fFeed.store +
       fCattle.house +
       fFeed.present +
       fGrass.silage + fCereal.silage + fHayStraw +
       fCereal.grains + fConcentrates + fSugarbeet +
       fMolasses,
       family = binomial,
       data = Badgers)
drop1(B.glm, test = "Chi")

B2.glm <- glm(Signs_in_yard~No_setts_in_fields + fFeed.store,
       family = binomial,
       data = Badgers)
drop1(B2.glm, test = "Chi")

summary(B2.glm)



library(geepack)
B.gee<-geese2(Signs_in_yard~No_setts_in_fields+
             fFeed.store,
             family=binomial,
             id=farm_code_numeric,corstr = "ar1",
             waves=Survey,
             data = Badgers)
summary(B.gee)




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




library(lme4)
B.glmm<-lmer(Signs_in_yard~No_setts_in_fields+
             fFeed.store+
            (1|farm_code_numeric),family=binomial,
            data=Badgers)

summary(B.glmm)






#FULL SELECTION
#
B.gee <- geese2(Signs_in_yard~fSeason+
       No_setts_in_fields+No_buildings +
       fFeed.store +
       fCattle.house +
       fFeed.present +
       fGrass.silage + fCereal.silage + fHayStraw +
       fCereal.grains + fConcentrates + fSugarbeet +
       fMolasses,
       family=binomial,
       id=farm_code_numeric,corstr = "ar1",
       waves=Survey,
       data = Badgers)
 summary(B.gee)
 
 #Drop least significant term
 B.gee1 <- geese2(Signs_in_yard~
       No_setts_in_fields+
       fSugarbeet,
       family=binomial,
       id=farm_code_numeric,corstr = "ar1",
       waves=Survey,
       data = Badgers)
 summary(B.gee1)


library(lme4)
B.glmm<-lmer(Signs_in_yard~fSeason++
       No_setts_in_fields+No_buildings +
       fFeed.store +
       fCattle.house +
       fFeed.present +
       fGrass.silage + fCereal.silage + fHayStraw +
       fCereal.grains + fConcentrates + fSugarbeet +
       fMolasses+
       (1|farm_code_numeric),family=binomial,
       data=Badgers)

summary(B.glmm)


B.glmm<-lmer(Signs_in_yard~
       No_setts_in_fields+  fSugarbeet+
       (1|farm_code_numeric),family=binomial,
       data=Badgers)

summary(B.glmm)
