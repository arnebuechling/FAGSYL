#   Central Europe Tree Growth Project
#
#   Functions for full models with climate and competition
#
### NOTE:  these models use the simpler 1 and 2 parameter functions for size and age effects
### NOTE:   RESCALING OF SIZEX0 AND AGEX0 PARAMETERS
#
#   Date: 3 April 2020
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

####################  Model 1 - just size, age, and PG terms

pg_size_age_model <- function(PG,siteplot,sizeX0,sizeXb,ageX0)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)     # sizeXb bounded between 0 1n 1,  sizeX0 positive
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)                      # ageX0 positive
   PG[siteplot] * size.effect * age.effect
}

####################  Model 2 - PG, size, age, and neighborhood competition (with lambdas)

pg_size_age_nci_model <- function(PG,siteplot,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   PG[siteplot] * size.effect * age.effect * competition.effect
}

####################  Model 3 - PG, size, age, neighborhood competition (with lambdas), and climate

pg_size_age_nci_climate_model <- function(PG,siteplot,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   temp.effect <- temp.a*exp(-0.5*((working$temp_k - temp.c)/temp.b)^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k - temp.lag.c)/temp.lag.b)^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip - prec.c)/prec.b)^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1 - prec.lag.c)/prec.lag.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect *
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}


####################  Model 3N - Adds nitrogen effect to Model 3

model_3N <- function(PG,siteplot,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          n.b,n.c,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   nitr.effect <- exp(-0.5*((working$N_all_MATCH2_wat_yr_kg_ha_yr - n.c)/n.b)^2)
   temp.effect <- temp.a*exp(-0.5*((working$temp_k - temp.c)/temp.b)^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k - temp.lag.c)/temp.lag.b)^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip - prec.c)/prec.b)^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1 - prec.lag.c)/prec.lag.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect * nitr.effect *
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 4 - Local differentiation in mode (c parameter) of climate response

climate_mode_ld_model <- function(PG,siteplot,region,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   temp.effect <- temp.a*exp(-0.5*((working$temp_k - temp.c[region])/temp.b)^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k - temp.lag.c[region])/temp.lag.b)^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip - prec.c[region])/prec.b)^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1 -prec.lag.c[region])/prec.lag.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect * 
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 4N - Local differentiation in mode (c parameter) of climate response

model_4N <- function(PG,siteplot,region,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          n.b,n.c,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   nitr.effect <- exp(-0.5*((working$N_all_MATCH2_wat_yr_kg_ha_yr - n.c)/n.b)^2)
   temp.effect <- temp.a*exp(-0.5*((working$temp_k - temp.c[region])/temp.b)^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k - temp.lag.c[region])/temp.lag.b)^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip - prec.c[region])/prec.b)^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1 -prec.lag.c[region])/prec.lag.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect * nitr.effect *
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 4B - Local differentiation in both mode and breadth of climate response (b and c parameters)

climate_mode_ld_model_4B <- function(PG,siteplot,region,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   temp.effect <- temp.a*exp(-0.5*((working$temp_k - temp.c[region])/temp.b[region])^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k - temp.lag.c[region])/temp.lag.b[region])^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip - prec.c[region])/prec.b[region])^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1 -prec.lag.c[region])/prec.lag.b[region])^2)

   PG[siteplot] * size.effect * age.effect * competition.effect * 
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 4BN - Local differentiation in both mode and breadth of climate response (b and c parameters)

model_4BN <- function(PG,siteplot,region,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          n.b,n.c,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   nitr.effect <- exp(-0.5*((working$N_all_MATCH2_wat_yr_kg_ha_yr - n.c)/n.b)^2)
   temp.effect <- temp.a*exp(-0.5*((working$temp_k - temp.c[region])/temp.b[region])^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k - temp.lag.c[region])/temp.lag.b[region])^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip - prec.c[region])/prec.b[region])^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1 -prec.lag.c[region])/prec.lag.b[region])^2)

   PG[siteplot] * size.effect * age.effect * competition.effect *  nitr.effect *
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 4ABC - Local differentiation in all 3 parameters  (a, b, and c)

climate_mode_ld_model_4ABC <- function(PG,siteplot,region,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   temp.effect <- temp.a[region]*exp(-0.5*((working$temp_k - temp.c[region])/temp.b[region])^2)
   temp.lag1.effect <- temp.lag.a[region]*exp(-0.5*((working$temp_lag1_k - temp.lag.c[region])/temp.lag.b[region])^2)
   prec.effect <- prec.a[region] * exp(-0.5*((working$precip - prec.c[region])/prec.b[region])^2)
   prec.lag1.effect <- prec.lag.a[region] * exp(-0.5*((working$precip_lag1 -prec.lag.c[region])/prec.lag.b[region])^2)

   PG[siteplot] * size.effect * age.effect * competition.effect * 
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 4ABCN - Local differentiation in all 3 parameters  (a, b, and c)

model_4ABCN <- function(PG,siteplot,region,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          n.b,n.c,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   nitr.effect <- exp(-0.5*((working$N_all_MATCH2_wat_yr_kg_ha_yr - n.c)/n.b)^2)
   temp.effect <- temp.a[region]*exp(-0.5*((working$temp_k - temp.c[region])/temp.b[region])^2)
   temp.lag1.effect <- temp.lag.a[region]*exp(-0.5*((working$temp_lag1_k - temp.lag.c[region])/temp.lag.b[region])^2)
   prec.effect <- prec.a[region] * exp(-0.5*((working$precip - prec.c[region])/prec.b[region])^2)
   prec.lag1.effect <- prec.lag.a[region] * exp(-0.5*((working$precip_lag1 -prec.lag.c[region])/prec.lag.b[region])^2)

   PG[siteplot] * size.effect * age.effect * competition.effect *  nitr.effect *
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}


####################  Model 5 - Climate response to interannual variation around local mean climate

#  NOTE:  uses the difference between absolute climate and 20-year mean as the independent variable

climate_ld_diff_model <- function(PG,siteplot,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   temp.effect <- temp.a*exp(-0.5*((working$temp_k_diff - temp.c)/temp.b)^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k_diff - temp.lag.c)/temp.lag.b)^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip_diff - prec.c)/prec.b)^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1_diff -prec.lag.c)/prec.lag.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect * 
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 5N - Climate response to interannual variation around local mean climate

#  NOTE:  uses the difference between absolute climate and 20-year mean as the independent variable

model_5N <- function(PG,siteplot,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          n.b,n.c,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   nitr.effect <- exp(-0.5*((working$N_all_MATCH2_wat_yr_kg_ha_yr - n.c)/n.b)^2)
   temp.effect <- temp.a*exp(-0.5*((working$temp_k_diff - temp.c)/temp.b)^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k_diff - temp.lag.c)/temp.lag.b)^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip_diff - prec.c)/prec.b)^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1_diff -prec.lag.c)/prec.lag.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect * nitr.effect *
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 5a - Climate response to interannual variation around local mean climate
#                          but with regional differentiation in magnitude of climate effect ("a" parameters)

#  NOTE:  uses the difference between absolute climate and 20-year mean as the independent variable

climate_ld_diff_model_5A <- function(PG,siteplot,region,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   temp.effect <- temp.a[region]*exp(-0.5*((working$temp_k_diff - temp.c)/temp.b)^2)
   temp.lag1.effect <- temp.lag.a[region]*exp(-0.5*((working$temp_lag1_k_diff - temp.lag.c)/temp.lag.b)^2)
   prec.effect <- prec.a[region] * exp(-0.5*((working$precip_diff - prec.c)/prec.b)^2)
   prec.lag1.effect <- prec.lag.a[region] * exp(-0.5*((working$precip_lag1_diff -prec.lag.c)/prec.lag.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect * 
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 5AN - Climate response to interannual variation around local mean climate
#                          but with regional differentiation in magnitude of climate effect ("a" parameters)

#  NOTE:  uses the difference between absolute climate and 20-year mean as the independent variable

model_5AN <- function(PG,siteplot,region,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          n.b,n.c,
                                          temp.a, temp.b, temp.c, temp.lag.a, temp.lag.b, temp.lag.c,
                                          prec.a, prec.b, prec.c, prec.lag.a, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   nitr.effect <- exp(-0.5*((working$N_all_MATCH2_wat_yr_kg_ha_yr - n.c)/n.b)^2)
   temp.effect <- temp.a[region]*exp(-0.5*((working$temp_k_diff - temp.c)/temp.b)^2)
   temp.lag1.effect <- temp.lag.a[region]*exp(-0.5*((working$temp_lag1_k_diff - temp.lag.c)/temp.lag.b)^2)
   prec.effect <- prec.a[region] * exp(-0.5*((working$precip_diff - prec.c)/prec.b)^2)
   prec.lag1.effect <- prec.lag.a[region] * exp(-0.5*((working$precip_lag1_diff -prec.lag.c)/prec.lag.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect *  nitr.effect *
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 5B - Climate response to interannual variation around local mean climate
#                          but with variation in magnitude of climate effect ("a" parameters) as a function
#                          of variation in tree-specific long-term mean conditions (rather than regional means)
#
#                    This will require a new set of variables to set the a parameter as a gaussian function of mean climate

climate_ld_diff_model_5B <- function(PG,siteplot,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          temp.a.prime, temp.b.prime,temp.c.prime,
                                          temp.lag.a.prime, temp.lag.b.prime, temp.lag.c.prime,
                                          temp.b, temp.c, temp.lag.b, temp.lag.c,
                                          prec.a.prime, prec.b.prime,prec.c.prime,
                                          prec.lag.a.prime, prec.lag.b.prime, prec.lag.c.prime,
                                          prec.b, prec.c, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
   competition.effect <- exp(-(C/100) * (working$size_ratio^gamma) * (nci^D))
   temp.a <- temp.a.prime * exp(-0.5*((working$mean_temp - temp.c.prime)/temp.b.prime)^2)
   temp.effect <- temp.a*exp(-0.5*((working$temp_k_diff - temp.c)/temp.b)^2)
   temp.lag.a <- temp.lag.a.prime * exp(-0.5*((working$mean_temp - temp.lag.c.prime)/temp.lag.b.prime)^2)
   temp.lag1.effect <- temp.lag.a*exp(-0.5*((working$temp_lag1_k_diff - temp.lag.c)/temp.lag.b)^2)

   prec.a <- prec.a.prime * exp(-0.5*((working$mean_precip - prec.c.prime)/prec.b.prime)^2)
   prec.effect <- prec.a * exp(-0.5*((working$precip_diff - prec.c)/prec.b)^2)
   prec.lag.a <- prec.lag.a.prime * exp(-0.5*((working$mean_precip - prec.c.prime)/prec.b.prime)^2)
   prec.lag1.effect <- prec.lag.a * exp(-0.5*((working$precip_lag1_diff -prec.lag.c)/prec.lag.b)^2)

   PG[siteplot] * size.effect * age.effect * competition.effect * 
       (temp.effect + temp.lag1.effect) * (prec.effect + prec.lag1.effect)
}

####################  Model 5BN - Climate response to interannual variation around local mean climate
#                          but with variation in magnitude of climate effect ("a" parameters) as a function
#                          of variation in tree-specific long-term mean conditions (rather than regional means)
#
#                    This will require a new set of variables to set the a parameter as a gaussian function of mean climate

model_5BN <- function(PG,siteplot,sizeX0,sizeXb,ageX0,alpha,beta,C,gamma,D,lambda,
                                          n.b,n.c,
                                          temp.a.prime, temp.b.prime,temp.c.prime,
                                          temp.lag.a.prime, temp.lag.b.prime, temp.lag.c.prime,
                                          temp.b, temp.c, temp.lag.b, temp.lag.c,
                                          prec.a.prime, prec.b.prime,prec.c.prime,
                                          prec.lag.a.prime, prec.lag.b.prime, prec.lag.c.prime,
                                          prec.b, prec.c, prec.lag.b, prec.lag.c)
{
   size.effect <-  1 - sizeXb * exp(-1.0*(sizeX0/100)*working$dbh_cm)
   age.effect  <-  exp(-1.0*(ageX0/100)*working$age)
   nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
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
