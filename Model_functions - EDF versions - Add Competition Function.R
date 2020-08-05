#########################################################################
#
#   Central Europe Tree Growth Project
#
#   Functions for full models with climate and competition
#
### NOTE:  these models use the simpler 1 and 2 parameter functions for size and age effects
### NOTE:   RESCALING OF SIZEX0 AND AGEX0 PARAMETERS
#
#   Date: 27 July 2020 Let competition effect vary by temperature in best models
#
#########################################################################
#
# Model SD as a linear function of mean (with intercept) if residuals are heteroscadistic
linear_dnorm_with_intercept <- function(x,mean,int,sigma,log=T)
{ sd <- int + (mean*sigma)
  dnorm(x,mean,sd,log)
}

# A power function version - seems to fall into a likelihood surface trap
power_dnorm <- function(x,mean,int,sigma,log=T)    # a new pdf to allow stabilization of variance at high
{ sd <- int*(mean^sigma)                           #   predicted values
  dnorm(x,mean,sd,log)
}


####################  Model 4BN - Local differentiation in both mode and breadth of climate response (b and c parameters)
# Competition response depends on temperature

model_4BN_vary_comp <- function(PG,siteplot,region,sizeX0,sizeXb,ageX0,
					  CX0,CXb, alpha,beta,gamma,D,lambda,
					  n.b,n.c,
					  temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
					  prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c )
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   C <- CX0 + CXb * working$tave_ann_hydro_k
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   nitr.effect <- exp(-0.5*((working$N_all_MATCH2_wat_yr_kg_ha_yr - n.c)/n.b)^2)
   temp.effect <- temp.a*exp(-0.5*((working$temp_k - temp.c[region])/temp.b[region])^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k - temp.lag.c[region])/temp.lag.b[region])^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip - prec.c[region])/prec.b[region])^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1 -prec.lag.c[region])/prec.lag.b[region])^2)

   PG[siteplot] * size.effect * age.effect * competition.effect *  nitr.effect *
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}


####################  Model 5BN - Climate response to interannual variation around local mean climate
#                          but with variation in magnitude of climate effect ("a" parameters) as a function
#                          of variation in tree-specific long-term mean conditions (rather than regional means)
#
#                    This will require a new set of variables to set the a parameter as a gaussian function of mean climate
#					 Two year climate model
#
#                    Competition response depends on temperature

model_5BN_vary_comp <- function(PG,siteplot,sizeX0,sizeXb,ageX0,
					  CX0,CXb,alpha,beta,gamma,D,lambda,
					  n.b,n.c,
					  temp.a.prime, temp.b.prime,temp.c.prime,
					  temp.lag.a.prime, temp.lag.b.prime, temp.lag.c.prime,
					  temp.b, temp.c, temp.lag.b, temp.lag.c,
					  prec.a.prime, prec.b.prime,prec.c.prime,
					  prec.lag.a.prime, prec.lag.b.prime, prec.lag.c.prime,
					  prec.b, prec.c, prec.lag.b, prec.lag.c )
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   C <- CX0 + CXb * working$tave_ann_hydro_k
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   nitr.effect <- exp(-0.5*((working$N_all_MATCH2_wat_yr_kg_ha_yr - n.c)/n.b)^2)
   temp.a <- temp.a.prime * exp(-0.5*((working$mean_temp - temp.c.prime)/temp.b.prime)^2)
   temp.effect <- temp.a*exp(-0.5*((working$temp_k_diff - temp.c)/temp.b)^2)
   temp.lag.a <- temp.lag.a.prime * exp(-0.5*((working$mean_temp - temp.lag.c.prime)/temp.lag.b.prime)^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k_diff - temp.lag.c)/temp.lag.b)^2)

   prec.a <- prec.a.prime * exp(-0.5*((working$mean_precip - prec.c.prime)/prec.b.prime)^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip_diff - prec.c)/prec.b)^2)
   prec.lag.a <- prec.lag.a.prime * exp(-0.5*((working$mean_precip - prec.c.prime)/prec.b.prime)^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1_diff -prec.lag.c)/prec.lag.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect *  nitr.effect *
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}


####################  Model 5BN Current Year - Climate response to interannual variation around local mean climate
#                          but with variation in magnitude of climate effect ("a" parameters) as a function
#                          of variation in tree-specific long-term mean conditions (rather than regional means)
#
#                    This will require a new set of variables to set the a parameter as a gaussian function of mean climate
#					 Current year clinate effects only

model_5BN_1yr_vary_comp <- function(PG,siteplot,sizeX0,sizeXb,ageX0,
                          alpha,beta,gamma,D,lambda,CX0,CXb,
						  n.b,n.c,
						  temp.a.prime, temp.b.prime,temp.c.prime,
						  temp.b, temp.c,
						  prec.a.prime, prec.b.prime,prec.c.prime,
						  prec.b, prec.c  )
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   C <- CX0 + CXb * working$tave_ann_hydro_k
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   nitr.effect <- exp(-0.5*((working$N_all_MATCH2_wat_yr_kg_ha_yr - n.c)/n.b)^2)
   temp.a <- temp.a.prime * exp(-0.5*((working$mean_temp - temp.c.prime)/temp.b.prime)^2)
   temp.effect <- temp.a*exp(-0.5*((working$temp_k_diff - temp.c)/temp.b)^2)
   prec.a <- prec.a.prime * exp(-0.5*((working$mean_precip - prec.c.prime)/prec.b.prime)^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip_diff - prec.c)/prec.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect *  nitr.effect * temp.effect * prec.effect
}






















