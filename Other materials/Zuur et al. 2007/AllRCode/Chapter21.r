#    Mixed Effects Models and Extensions in Ecology with R (2009)
#    Zuur, Ieno, Walker, Saveliev and Smith.    Springer

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.



#This code was written by: Jonathan Rhodes (j.rhodes@uq.edu.au)
##############################################################

library(AED); data(Koalas)

cor(Koalas[, 6:22], method = "spearman")


#New variables: phss
Koalas$phss_2.5km_new <- Koalas[,"phss_2.5km"] -
              Koalas[,"phss_5km"]
Koalas$phss_1km_new <- Koalas[,"phss_1km"] -
              Koalas[,"phss_2.5km"]


#pm
Koalas$pm_2.5km_new <- Koalas[,"pm_2.5km"] -
              Koalas[,"pm_5km"]
Koalas$pm_1km_new <- Koalas[,"pm_1km"] -
              Koalas[,"pm_2.5km"]


#pdens
Koalas$pdens_2.5km_new <- Koalas[,"pdens_2.5km"] -
              Koalas[,"pdens_5km"]
Koalas$pdens_1km_new <- Koalas[,"pdens_1km"] -
              Koalas[,"pdens_2.5km"]



#edens
Koalas$edens_2.5km_new <- Koalas[,"edens_2.5km"] -
              Koalas[,"edens_5km"]
Koalas$edens_1km_new <- Koalas[,"edens_1km"] -
              Koalas[,"edens_2.5km"]


#rdens
Koalas$rdens_2.5km_new <- Koalas[,"rdens_2.5km"] -
              Koalas[,"rdens_5km"]
Koalas$rdens_1km_new <- Koalas[,"rdens_1km"] -
              Koalas[,"rdens_2.5km"]



cor(Koalas[,c("phss_5km","phss_2.5km_new",
        "phss_1km_new")],method="spearman")
        
        
        
Glm_5km<- glm(presence~pprim_ssite+psec_ssite+phss_5km+
            phss_2.5km_new+phss_1km_new+pm_5km+
            pm_2.5km_new+pm_1km_new+pdens_5km+
            pdens_2.5km_new+pdens_1km_new+
            rdens_5km+rdens_2.5km_new+rdens_1km_new,
           data=Koalas,family=binomial)
library(Design)
vif(Glm_5km)



Glm_2.5km <- glm(presence ~ pprim_ssite + psec_ssite +
               phss_2.5km + phss_1km_new + pm_2.5km +
               pm_1km_new + pdens_2.5km +
               pdens_1km_new + rdens_2.5km +
               rdens_1km_new,
              data = Koalas, family = binomial)
vif(Glm_2.5km)




Glm_1km <- glm(presence ~ pprim_ssite + psec_ssite + phss_1km +
             pm_1km + pdens_1km + rdens_1km,
             data = Koalas, family = binomial)
vif(Glm_1km)



library(ncf)
Correlog<-spline.correlog(x=Koalas[,"easting"],
              y=Koalas[,"northing"],
              z=Koalas[,"presence"],xmax=10000)
plot.spline.correlog(Correlog)





Correlog_Glm_5km <-
              spline.correlog(x = Koalas[, "easting"],
              y = Koalas[, "northing"],
              z = residuals(Glm_5km, type = "pearson"),
              xmax = 10000)
plot.spline.correlog(Correlog_Glm_5km)




library(glmmML)
Glmm_5km<-glmmML(presence ~ pprim_ssite + psec_ssite +
          phss_5km + phss_2.5km_new + phss_1km_new +
          pm_5km + pm_2.5km_new + pm_1km_new +
          pdens_5km + pdens_2.5km_new + pdens_1km_new +
          rdens_5km + rdens_2.5km_new + rdens_1km_new,
          cluster = site, data = Koalas, family = binomial)


#Copy and paste the function pres.glmmML at the end of this file into R
Correlog.Glmm_5km <- spline.correlog(
       x = Koalas[, "easting"],
         y = Koalas[, "northing"],
         z = pres.glmmML(model = Glmm_5km, data = Koalas),
         xmax = 10000)
plot.spline.correlog(Correlog.Glmm_5km)




Koalas_St <- cbind(Koalas[, 1:5],
     apply(X = Koalas[, 6:ncol(Koalas)],
     MARGIN = 2, FUN = function(x){(x-mean(x))/sd(x)}))

Koalas_St <- cbind(Koalas[, 1:5],
     apply(X = Koalas[, 6:ncol(Koalas)],
     MARGIN = 2, FUN = scale))


glmmML(presence ~ pprim_ssite + psec_ssite + phss_1km +
       pm_1km, cluster = site,
       data = Koalas_St, family = binomial)








#fit the best linear model
Model_Best<-glmmML(formula=presence~pprim_ssite+psec_ssite+phss_1km+pm_1km,
             cluster=site,data=Koalas_St,family=binomial)

#QUANTILE-QUANTILE PLOTS (note that the functions "fitted.glmmML" and "res.glmmML" are requiured here)

#calculate the fitted values for the model
Fitted<-fitted.glmmML(model=Model_Best,data=Koalas_St)

#calculate the ordered residuals for the model
Resids<-sort(res.glmmML(model=Model_Best,data=Koalas_St))
Resids_0<-Resids[Resids<0]
Resids_1<-Resids[Resids>=0]

#specify the number of replicates
Reps<-1000

#simulate Reps data sets
Sims<-matrix(rbinom(n=length(Fitted)*Reps,size=1,p=Fitted),nrow=length(Fitted),ncol=Reps)

#fit the model to each simulated data set
Models<-apply(Sims,MARGIN=2,FUN=function(X){return(list(X,glmmML(formula=X~pprim_ssite+psec_ssite+phss_1km+pm_1km,cluster=site,data=Koalas_St,family=binomial)))})

#calculate the ordered (interpolated) simulated residuals and the point-wise 95% confidence intervals
Resids_Sim<-matrix(unlist(lapply(Models,FUN=function(X){TempData<-Koalas_St;TempData[,"presence"]<-X[[1]];return(sort(res.glmmML(X[[2]],TempData)))})),ncol=Reps,nrow=nrow(Koalas_St))
Resids_Sim_0<-matrix(apply(Resids_Sim,MARGIN=2,FUN=function(X){quantile(X[X<0],ppoints(Resids_0,a=1))}),ncol=Reps,nrow=length(Resids_0))
Resids_Sim_1<-matrix(apply(Resids_Sim,MARGIN=2,FUN=function(X){quantile(X[X>=0],ppoints(Resids_1,a=1))}),ncol=Reps,nrow=length(Resids_1))
Resids_Sim_0_Median<-apply(Resids_Sim_0,MARGIN=1,FUN=median)
Resids_Sim_1_Median<-apply(Resids_Sim_1,MARGIN=1,FUN=median)
Resids_Sim_0_Lower<-apply(Resids_Sim_0,MARGIN=1,FUN=function(X){quantile(X,0.025)})
Resids_Sim_0_Upper<-apply(Resids_Sim_0,MARGIN=1,FUN=function(X){quantile(X,0.975)})
Resids_Sim_1_Lower<-apply(Resids_Sim_1,MARGIN=1,FUN=function(X){quantile(X,0.025)})
Resids_Sim_1_Upper<-apply(Resids_Sim_1,MARGIN=1,FUN=function(X){quantile(X,0.975)})
#plot the qauntile-quantile plot with 95% confidence intervals and 1:1 line
plot(Resids_Sim_0_Median,Resids_0,xlim=c(-1,1),ylim=c(-1,1)),xlab="simulated quantiles",ylab="fitted quantiles")
points(Resids_Sim_1_Median,Resids_1)
lines(Resids_Sim_0_Median,Resids_Sim_0_Lower)
lines(Resids_Sim_0_Median,Resids_Sim_0_Upper)
lines(Resids_Sim_1_Median,Resids_Sim_1_Lower)
lines(Resids_Sim_1_Median,Resids_Sim_1_Upper)
abline(0,1,lty=3)

#PARTIAL RESIDUAL PLOTS (note that the functions "fitted.glmmML" and "res.glmmML" are requiured here)

#calculate the partial residuals for each covariate
Part_Res1<-res.glmmML(Model_Best,Koalas_St)/(fitted.glmmML(Model_Best,Koala_Data_St)*(1-fitted.glmmML(Model_Best,Koalas_St)))+Model_Best$coefficients[2]*Koalas_St[,"pprim_ssite"]
Part_Res2<-res.glmmML(Model_Best,Koalas_St)/(fitted.glmmML(Model_Best,Koala_Data_St)*(1-fitted.glmmML(Model_Best,Koalas_St)))+Model_Best$coefficients[3]*Koalas_St[,"psec_ssite"]
Part_Res3<-res.glmmML(Model_Best,Koalas_St)/(fitted.glmmML(Model_Best,Koala_Data_St)*(1-fitted.glmmML(Model_Best,Koalas_St)))+Model_Best$coefficients[4]*Koalas_St[,"phss_1km"]
Part_Res4<-res.glmmML(Model_Best,Koalas_St)/(fitted.glmmML(Model_Best,Koala_Data_St)*(1-fitted.glmmML(Model_Best,Koalas_St)))+Model_Best$coefficients[5]*Koalas_St[,"pm_1km"]

#plot the partial residuals and the smoothed plot
split.screen(c(2,2))
screen(1)
plot(Koalas_St[,"pprim_ssite"],Part_Res1,ylab="partial residuals")
lines(seq(min(Koalas_St[,"pprim_ssite"]),max(Koalas_St[,"pprim_ssite"]),length.out=100),predict(loess(formula=Part_Res1~pprim_ssite,data=Koala_Data_St),newdata=data.frame(pprim_ssite=seq(min(Koala_Data_St[,"pprim_ssite"]),max(Koalas_St[,"pprim_ssite"]),length.out=100))))
screen(2)
plot(Koalas_St[,"psec_ssite"],Part_Res2,ylab="partial residuals")
lines(seq(min(Koalas_St[,"psec_ssite"]),max(Koalas_St[,"psec_ssite"]),length.out=100),predict(loess(formula=Part_Res2~psec_ssite,data=Koala_Data_St),newdata=data.frame(psec_ssite=seq(min(Koala_Data_St[,"psec_ssite"]),max(Koala_Data_St[,"psec_ssite"]),length.out=100))))
screen(3)
plot(Koalas_St[,"phss_1km"],Part_Res3,ylab="partial residuals")
lines(seq(min(Koalas_St[,"phss_1km"]),max(Koalas_St[,"phss_1km"]),length.out=100),predict(loess(formula=Part_Res3~phss_1km,data=Koalas_St),newdata=data.frame(phss_1km=seq(min(Koalas_St[,"phss_1km"]),max(Koalas_St[,"phss_1km"]),length.out=100))))
screen(4)
plot(Koalas_St[,"pm_1km"],Part_Res4,ylab="partial residuals")
lines(seq(min(Koalas_St[,"pm_1km"]),max(Koalas_St[,"pm_1km"]),length.out=100),predict(loess(formula=Part_Res4~pm_1km,data=Koalas_St),newdata=data.frame(pm_1km=seq(min(Koalas_St[,"pm_1km"]),max(Koalas_St[,"pm_1km"]),length.out=100))))
close.screen(all = TRUE)


###########################################################
#Copy and paste all code below
###########################################################


#THIS FILE CONTAINS THE R CODE TO GENERATES THE EXAMPLE QUANTILE-QUANTILE
#PLOTS AND PARTIAL RESIDUAL PLOTS.

pres.glmmML<-function(model,data)
#This function outputs a vector of
#the Pearson residuals for a glmmML object
#- "model" is a glmmML object, and "data"
#is the data frame of the data to which "model"
#was fitted. This function is based on
#the glmmML package version 0.65-5 and
#may not work properly in earlier or later
#versions. Please report any bugs or
#errors to Jonathan Rhodes (j.rhodes@uq.edu.au)
{
	Model<-model
	Data<-data

	#get the model frame
	ModelFrame<-model.frame(formula(terms(Model)),data=Data)

	#get the design matrix
	DesignMatrix<-model.matrix(terms(Model),data=ModelFrame)

	#get the linear predictor
	LP<-DesignMatrix %*% as.matrix(Model$coefficients)

	#get the cluster
	Cluster<-Data[,toString(Model$call$cluster)]
	Cluster<-as.factor(Cluster)

	#add the random-effect posterior modes
	for (i in 1:length(levels(Cluster)))
	{
		LP[Cluster == levels(Cluster)[i]] = LP[Cluster == levels(Cluster)[i]] + Model$posterior.modes[i]
	}

	#get the response
	Response<-as.matrix(model.response(ModelFrame))

	if (is.null(Model$call$family) || (Model$call$family==as.name(paste("binomial(link = ",dQuote("logit"),")",sep=""))) || (Model$call$family==as.name("binomial")))
	#binomial response with logit link function
	{
		#get the fitted p values
		P<-exp(LP) / (1 + exp(LP))

		#get the pearson residuals
		if (ncol(Response) == 1)
		#bernoulli data
		{
			PRes<-(Response - P) / sqrt(P * (1 - P))
		}
		else if (ncol(Response) == 2)
		#binomial data
		{
			PRes<-(Response[,1] - (Response[,1] + Response[,2]) * P) / sqrt((Response[,1] + Response[,2]) * P * (1 - P))
		}
		else
		#response matrix the wrong size
		{
			stop("wrong number of response variables")
		}

	}
	else if (Model$call$family==as.name(paste("binomial(link = ",dQuote("cloglog"),")",sep="")))
	#binomial response with complementary log-log link
	{

		#get the fitted p values
		P<-1 - exp(-exp(LP))

		#get the pearson residuals
		if (ncol(Response) == 1)
		#bernoulli data
		{
			PRes<-(Response - P) / sqrt(P * (1 - P))
		}
		else if (ncol(Response) == 2)
		#binomial data
		{
			PRes<-(Response[,1] - (Response[,1] + Response[,2]) * P) / sqrt((Response[,1] + Response[,2]) * P * (1 - P))
		}
		else
		{
			stop("wrong number of response variables")
		}
	}
	else if ((Model$call$family==as.name(paste("poisson(link = ",dQuote("log"),")",sep=""))) || (Model$call$family==as.name("poisson")))
	#poisson response with log link function
	{

		#get the fitted count values
		C<-exp(LP)

		#get the pearson residuals
		PRes<-(Response - C) / sqrt(C)
	}
	else
	#response not binomial or poisson
	{
		stop("incompatible family specified")
	}

	return(as.vector(PRes))
}



fitted.glmmML<-function(model,data)
#This function outputs a vector of
#the fitted values for a glmmML object
#- "model" is a glmmML object, and "data"
#is the data frame of the data to which "model"
#was fitted. This function is based on
#the glmmML package version 0.65-5 and
#may not work properly in earlier or later
#versions. Please report any bugs or
#errors to Jonathan Rhodes (j.rhodes@uq.edu.au)
{
	Model<-model
	Data<-data

	#get the model frame
	ModelFrame<-model.frame(formula(terms(Model)),data=Data)

	#get the design matrix
	DesignMatrix<-model.matrix(terms(Model),data=ModelFrame)

	#get the linear predictor
	LP<-DesignMatrix %*% as.matrix(Model$coefficients)

	#get the cluster
	Cluster<-Data[,toString(Model$call$cluster)]
	Cluster<-as.factor(Cluster)

	#add the random-effect posterior modes
	for (i in 1:length(levels(Cluster)))
	{
		LP[Cluster == levels(Cluster)[i]] = LP[Cluster == levels(Cluster)[i]] + Model$posterior.modes[i]
	}

	#get the response
	Response<-as.matrix(model.response(ModelFrame))

	if (is.null(Model$call$family) || (Model$call$family==as.name(paste("binomial(link = ",dQuote("logit"),")",sep=""))) || (Model$call$family==as.name("binomial")))
	#binomial response with logit link function
	{
		#get the fitted values
		if (ncol(Response) == 1)
		#bernoulli data
		{
			Fit<-exp(LP) / (1 + exp(LP))
		}
		else if (ncol(Response) == 2)
		#binomial data
		{
			Fit<-(Response[,1] + Response[,2]) * exp(LP) / (1 + exp(LP))
		}
		else
		#response matrix the wrong size
		{
			stop("wrong number of response variables")
		}

	}
	else if (Model$call$family==as.name(paste("binomial(link = ",dQuote("cloglog"),")",sep="")))
	#binomial response with complementary log-log link
	{
		#get the fitted values
		if (ncol(Response) == 1)
		#bernoulli data
		{
			Fit<-1 - exp(-exp(LP))
		}
		else if (ncol(Response) == 2)
		#binomial data
		{
			Fit<-(Response[,1] + Response[,2]) * (1 - exp(-exp(LP)))
		}
		else
		{
			stop("wrong number of response variables")
		}
	}
	else if ((Model$call$family==as.name(paste("poisson(link = ",dQuote("log"),")",sep=""))) || (Model$call$family==as.name("poisson")))
	#poisson response with log link function
	{

		#get the fitted count values
		Fit<-exp(LP)
	}
	else
	#response not binomial or poisson
	{
		stop("incompatible family specified")
	}

	return(as.vector(Fit))
}




res.glmmML<-function(model,data)
#This function outputs a vector of
#the ordinary residuals for a glmmML object
#- "model" is a glmmML object, and "data"
#is the data frame of the data to which "model"
#was fitted. This function is based on
#the glmmML package version 0.65-5 and
#may not work properly in earlier or later
#versions. Please report any bugs or
#errors to Jonathan Rhodes (j.rhodes@uq.edu.au)
{
	Model<-model
	Data<-data

	#get the model frame
	ModelFrame<-model.frame(formula(terms(Model)),data=Data)

	#get the design matrix
	DesignMatrix<-model.matrix(terms(Model),data=ModelFrame)

	#get the linear predictor
	LP<-DesignMatrix %*% as.matrix(Model$coefficients)

	#get the cluster
	Cluster<-Data[,toString(Model$call$cluster)]
	Cluster<-as.factor(Cluster)

	#add the random-effect posterior modes
	for (i in 1:length(levels(Cluster)))
	{
		LP[Cluster == levels(Cluster)[i]] = LP[Cluster == levels(Cluster)[i]] + Model$posterior.modes[i]
	}

	#get the response
	Response<-as.matrix(model.response(ModelFrame))

	if (is.null(Model$call$family) || (Model$call$family==as.name(paste("binomial(link = ",dQuote("logit"),")",sep=""))) || (Model$call$family==as.name("binomial")))
	#binomial response with logit link function
	{
		#get the fitted p values
		P<-exp(LP) / (1 + exp(LP))

		#get the residuals
		if (ncol(Response) == 1)
		#bernoulli data
		{
			Res<-Response - P
		}
		else if (ncol(Response) == 2)
		#binomial data
		{
			Res<-Response[,1] - (Response[,1] + Response[,2]) * P
		}
		else
		#response matrix the wrong size
		{
			stop("wrong number of response variables")
		}

	}
	else if (Model$call$family==as.name(paste("binomial(link = ",dQuote("cloglog"),")",sep="")))
	#binomial response with complementary log-log link
	{

		#get the fitted p values
		P<-1 - exp(-exp(LP))

		#get the residuals
		if (ncol(Response) == 1)
		#bernoulli data
		{
			Res<-Response - P
		}
		else if (ncol(Response) == 2)
		#binomial data
		{
			Res<-Response[,1] - (Response[,1] + Response[,2]) * P
		}
		else
		{
			stop("wrong number of response variables")
		}
	}
	else if ((Model$call$family==as.name(paste("poisson(link = ",dQuote("log"),")",sep=""))) || (Model$call$family==as.name("poisson")))
	#poisson response with log link function
	{

		#get the fitted count values
		C<-exp(LP)

		#get the residuals
		Res<-Response - C
	}
	else
	#response not binomial or poisson
	{
		stop("incompatible family specified")
	}

	return(as.vector(Res))
}
