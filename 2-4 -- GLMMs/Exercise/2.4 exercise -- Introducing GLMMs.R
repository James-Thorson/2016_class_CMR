
# Function to generate data
Sim_Fn = function( n_per_group, n_groups=10, logSD=1, logmean=2, beta_u=1){
  # Simulate measured and unmeasured random components
  u_s = rnorm(n_groups)
  r_s = rnorm(n_groups, mean=logmean, sd=logSD)

  # Simulate samples
  s_i = rep( 1:n_groups, each=n_per_group)
  mean_s = exp( beta_u*u_s + r_s )
  c_i = rpois( n=n_per_group*n_groups, lambda=mean_s[s_i] )

  # Bundle and return
  DF = data.frame( "site"=s_i, "count"=c_i, "covariate"=u_s[s_i])
  return( DF )
}

# Illustrate function
Sim_Fn( n_per_group=10, n_groups=10, beta_u=1)
