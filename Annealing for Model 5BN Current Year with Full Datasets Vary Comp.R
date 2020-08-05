#############################################################################
#
#	Central Europe Tree Growth Project
#
#       Annealing with Full Datasets
#
#       Model 5B
#
#       CDC,  April 28, 2020
# 
#############################################################################
library(likelihood)

computer<-"azure"

if(computer == "Charlie")
{
  dat_dir <- "\\\\canhamC6\\c\\Users\\canhamc\\Documents\\Current Manuscripts\\Arne - Carpathians\\Earth Day Files"
  out_dir <- "\\\\canhamC6\\c\\Users\\canhamc\\Documents\\Current Manuscripts\\Arne - Carpathians\\Earth Day Files\\Output"
  code_dir <- "\\\\canhamC6\\c\\Users\\canhamc\\Documents\\Current Manuscripts\\Arne - Carpathians\\Earth Day Files\\Code"
} else {
  dat_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/ABIALB/"
  out_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/ABIALB/"
  code_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/ABIALB/"
}


### Get model functions
setwd(code_dir)
source("Model Functions - EDF versions.R")
source("Model_functions - EDF versions - Add Competition Function.R")

# focal species
spp_list <- c("PICABI","FAGSYL","ABIALB","ACEPSE")

index <- 2

### load working target data file
setwd(dat_dir)

load(paste(spp_list[index],"Full Working Dataset with All Precip Variables Nitrogen and Mean Climate.Rdata"))

### specify climate variables
working$temp_k <- working$tave_ann_hydro_k
working$temp_lag1_k <- working$tave_ann_hydro_k_lag1

### best precip variables

if (index %in% c(1,4))     # Picea and Acer use annual precip
{  working$precip <- working$ppt_ann_hydro_mm
   working$precip_lag1 <- working$ppt_ann_lag1_hydro_mm
}

if (index == 2)      # Fagus uses water deficit
{  working$precip <- working$hydro_yr_water_deficit_mm
   working$precip_lag1 <- working$hydro_yr_water_deficit_lag1_mm
}
if (index == 3)      # Abies uses seasonal precip
{  working$precip <- working$seas_prec_hydro_yr_mm
   working$precip_lag1 <- working$seas_prec_hydro_yr_lag1_mm
}


working$temp_k_diff <- working$temp_k - working$mean_temp
working$temp_lag1_k_diff <- working$temp_lag1_k - working$mean_temp

### best precip variables

if (index %in% c(1,4))     # Picea and Acer use annual precip
{  working$precip_diff <- working$precip - working$mean_ann_precip
   working$precip_lag1_diff <- working$precip_lag1 - working$mean_ann_precip
   working$mean_precip <- working$mean_ann_precip
}

if (index == 2)      # Fagus uses water deficit
{  working$precip_diff <- working$precip - working$mean_wd
   working$precip_lag1_diff <- working$precip_lag1 - working$mean_wd
   working$mean_precip <- working$mean_wd
}

if (index == 3)      # Abies uses seasonal precip
{  working$precip_diff <- working$precip - working$mean_seas_precip
   working$precip_lag1_diff <- working$precip_lag1 - working$mean_seas_precip
   working$mean_precip <- working$mean_seas_precip
}


###################################

iterations <- 15000

setwd(out_dir)

if (index == 4) { load(paste(spp_list[index],"Model 5BN Full Azure Results.Rdata")) } else
                { load(paste(spp_list[index],"Model 5BN Full Azure Results.Rdata")) }


var <- list(mean = "predicted", x = "incr_mm", siteplot = "stdcode", log=T)

# set parameter limits
if (index == 4) { par <- model_5BN_full_azure_results_2$best_pars } else
                { par <- model_5BN_full_azure_results$best_pars  }


if (index == 1) { nspp <- 7 } else
                { nspp <- 5 } 

# remove C parameter and estimate with a linear function
par <- par[!names(par) %in% "C" ]

# set parameter limits
par <- list(PG = rep(50,num_site_plot),
		sizeX0 = 0.1, sizeXb = 0.5,
		ageX0 = 1,
 		int = 0, sigma = 1,
         	alpha = 2, beta = 1, gamma = 1, D = 2, lambda = rep(0.5,nspp),

		temp.a.prime = 0.5,   # the climate "a" parameters should be bounded between 0 and 1
                temp.b.prime = 5,
                temp.c.prime = mean(working$temp_k),

		temp.b = 5,    # the climate "b" parameters are the variances of the gaussian - 0 to large
                temp.c = 0,  # the "c" parameters are the mode of the gaussian

		prec.a.prime = 0.5,   # the climate "a" parameters should be bounded between 0 and 1
                prec.b.prime = 25,
                prec.c.prime = mean(working$precip),

                prec.b = 10, 
                prec.c = 0
	   )
	   
par_lo <- list(PG = rep(0,num_site_plot),
		sizeX0 = 0, sizeXb = 0,
		ageX0 = 0,
		int = 0, sigma = 0,
        	alpha = 0, beta = 0, gamma = -2, D = 1, lambda = rep(0,nspp),

		temp.a.prime = 0,   # the climate "a" parameters should be bounded between 0 and 1
                temp.b.prime = 0.5,
                temp.c.prime = min(working$temp_k),

		temp.b = 0.5,     # the climate "b" parameters are the variances of the gaussian - small to large
                temp.c = min(working$temp_k_diff),  # the "c" parameters are the mode of the gaussian

		prec.a.prime = 0,   # the climate "a" parameters should be bounded between 0 and 1
                prec.b.prime = 0.5, 
                prec.c.prime = min(working$precip),

                prec.b = 0.5, 
                prec.c = min(working$precip_diff)
	   )

par_hi <- list(PG = rep(250,num_site_plot),
		sizeX0 = 100, sizeXb = 1,
		ageX0 = 1,
                int = 100, sigma = 2,
        	alpha = 4, beta = 3, gamma = 2, D = 3, lambda = rep(1,nspp),

		temp.a.prime = 1,   # the climate "a" parameters should be bounded between 0 and 1
                temp.b.prime = 5000, 
                temp.c.prime = max(working$temp_k), 

		temp.b = 5000,      # the climate "b" parameters are the variances of the gaussian - 0 to large
                temp.c = max(working$temp_k_diff),  # the "c" parameters are the mode of the gaussian

		prec.a.prime = 1,    # the climate "a" parameters should be bounded between 0 and 1
                prec.b.prime = 10000, 
                prec.c.prime = max(working$precip),

                prec.b = 10000, 
                prec.c = max(working$precip_diff)
                )

par$PG <- model_5BN_full_azure_results$best_pars$PG
par$sizeX0 <- model_5BN_full_azure_results$best_pars$sizeX0
par$sizeXb <- model_5BN_full_azure_results$best_pars$sizeXb
par$ageX0 <- model_5BN_full_azure_results$best_pars$ageX0
par$int <- model_5BN_full_azure_results$best_pars$int
par$sigma <- model_5BN_full_azure_results$best_pars$sigma
par$alpha <- model_5BN_full_azure_results$best_pars$alpha
par$beta <- model_5BN_full_azure_results$best_pars$beta
#par$C <- model_5BN_full_azure_results$best_pars$C
par$gamma <- model_5BN_full_azure_results$best_pars$gamma
par$D <- model_5BN_full_azure_results$best_pars$D
par$lambda <- model_5BN_full_azure_results$best_pars$lambda[1:nspp]
par$temp.a.prime <- model_5BN_full_azure_results$best_pars$temp.a.prime
par$temp.b.prime <- model_5BN_full_azure_results$best_pars$temp.b.prime
par$temp.c.prime <- model_5BN_full_azure_results$best_pars$temp.c.prime
par$temp.b <- model_5BN_full_azure_results$best_pars$temp.b
par$temp.c <- model_5BN_full_azure_results$best_pars$temp.c
par$prec.a.prime <- model_5BN_full_azure_results$best_pars$prec.a.prime
par$prec.b.prime <- model_5BN_full_azure_results$best_pars$prec.b.prime
par$prec.c.prime <- model_5BN_full_azure_results$best_pars$prec.c.prime
par$prec.b <- model_5BN_full_azure_results$best_pars$prec.b
par$prec.c <- model_5BN_full_azure_results$best_pars$prec.c
par$n.c <- model_5BN_full_azure_results$best_pars$n.c
par$n.b <- model_5BN_full_azure_results$best_pars$n.b
par_lo$n.c <- 0
par_lo$n.b <- 0
par_hi$n.c <- 40
par_hi$n.b <- 200

par$CX0 <- 1
par$CXb <- 1
par_lo$CX0 <- -200
par_lo$CXb <- -5
par_hi$CX0 <- 200
par_hi$CXb <- 5


model_5BN_full_1year_results <- anneal(model=model_5BN_1yr_vary_comp, par=par, var=var, source_data=working, par_lo=par_lo, par_hi=par_hi,
			  pdf=linear_dnorm_with_intercept, dep_var="incr_mm", max_iter=iterations)

setwd(out_dir)

write_results(model_5BN_full_1year_results, file = paste(spp_list[index],"Model 5BN Full 1year Azure Results.txt", sep=" "), data=F, print_whole_hist=F)
save(model_5BN_full_1year_results, file=paste(spp_list[index],"Model 5BN Full 1year Azure Results.Rdata",sep=" "))






