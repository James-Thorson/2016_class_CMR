
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2016_classes_private/CMR models/2-2 -- GLMs/Afternoon/" )

###### Example GLm code

#### Method 1 -- Nonlinear optimization
NegLogLike_Fn = function(Par, Data){
  # Parameters
  Mean_hat = Par[1]
  Size_hat = Par[2]
  # Log-likelihood
  LogLike_i = dnbinom( Data$Counts, mu=exp(Mean_hat), size=Size_hat, log=TRUE )
  NegLogLike = -1 * sum(LogLike_i)
  return( NegLogLike )
}

DF = NULL
for(i in 1:1000){
  # simulate example data
  log_TrueMean = log(3)
  TrueSize = 1
  Counts = rnbinom(100, mu=exp(log_TrueMean), size=TrueSize) # Var = mu + mu^2/size

  # Run model
  Data = list( 'Counts'=Counts )
  Start = c(1,1)
  NegLogLike_Fn( Par=Start, Data=Data)
  Opt = optim( par=Start, fn=NegLogLike_Fn, Data=Data, lower=c(-Inf,0.01), upper=Inf, method="L-BFGS-B", hessian=TRUE )
  # Estimated parameters
  Hat =  Opt$par[1]
  SE = sqrt(diag( solve(Opt$hessian) ))[1] # square root of diagonal elements of the inverse-hessian matrix
  CI_mult = qnorm(1 - 0.05/2)
  InsideTF = ifelse( (Hat + CI_mult*SE)>log_TrueMean & (Hat - CI_mult*SE)<log_TrueMean, TRUE, FALSE)
  DF = rbind(DF, c(Hat, SE, InsideTF))
}
mean(DF[,3])
