#############################################################################
#
#	Central Europe Tree Growth Project
#
#       Annealing with Full Datasets
#
#       Model 5BN
#       Allow competition C parameter to vary with temperature
#
#       27 July 2020
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
  dat_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/"
  out_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/Output_azure/"
  code_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/Code/Annealing_code/"
}

### Get model functions
setwd(code_dir)
source("Model Functions - EDF versions - Add current year models.R")
source("Model_functions - EDF versions - Add Competition Function.R")

# focal species
spp_list <- c("PICABI","FAGSYL","ABIALB","ACEPSE")

index <- 4

### load working target data file
setwd(paste(dat_dir,spp_list[index],sep=""))

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

iterations <- 100000

setwd(out_dir)

if (index == 4) { load(paste(spp_list[index],"Model 5BN Full Azure Results 2.Rdata")) } else
                { load(paste(spp_list[index],"Model 5BN Full Azure Results.Rdata")) }


var <- list(mean = "predicted", x = "incr_mm", siteplot = "stdcode", log=T)

# set parameter limits
if (index == 4) { par <- model_5BN_full_azure_results_2$best_pars } else
                { par <- model_5BN_full_azure_results$best_pars  }

if (index == 1) { nspp <- 7 } else
                { nspp <- 5 } 
par$lambda <- par$lambda[1:nspp]

# remove C parameter and estimate with a linear function
par <- par[!names(par) %in% "C" ]

par_lo <- list(PG = rep(0,num_site_plot),
		sizeX0 = 0, sizeXb = 0,
		ageX0 = 0,
		int = 0, sigma = 0,
        	alpha = 0, beta = 0, gamma = -2, D = 1, lambda = rep(0,nspp),

		temp.a.prime = 0, temp.lag.a.prime = 0,   # the climate "a" parameters should be bounded between 0 and 1
                temp.b.prime = 0.5, temp.lag.b.prime = 0.5,
                temp.c.prime = min(working$temp_k), temp.lag.c.prime =min(working$temp_k),

		temp.b = 5, temp.lag.b = 5,     # the climate "b" parameters are the variances of the gaussian - small to large
                temp.c = min(working$temp_k_diff), temp.lag.c = min(working$temp_k_diff),  # the "c" parameters are the mode of the gaussian

		prec.a.prime = 0, prec.lag.a.prime = 0,   # the climate "a" parameters should be bounded between 0 and 1
                prec.b.prime = 0.5, prec.lag.b.prime = 0.5,
                prec.c.prime = min(working$precip), prec.lag.c.prime =min(working$precip),

                prec.b = 5, prec.lag.b = 5,
                prec.c = min(working$precip_diff), prec.lag.c = min(working$precip_diff)
	   )

par_hi <- list(PG = rep(250,num_site_plot),
		sizeX0 = 100, sizeXb = 1,
		ageX0 = 1,
                int = 100, sigma = 2,
        	alpha = 4, beta = 3,  gamma = 2, D = 3, lambda = rep(1,nspp),

		temp.a.prime = 1, temp.lag.a.prime = 1,   # the climate "a" parameters should be bounded between 0 and 1
                temp.b.prime = 5000, temp.lag.b.prime = 5000,
                temp.c.prime = max(working$temp_k), temp.lag.c.prime =max(working$temp_k),

		temp.b = 5000, temp.lag.b = 5000,     # the climate "b" parameters are the variances of the gaussian - 0 to large
                temp.c = max(working$temp_k_diff), temp.lag.c = max(working$temp_k_diff),  # the "c" parameters are the mode of the gaussian

		prec.a.prime = 1, prec.lag.a.prime = 1,   # the climate "a" parameters should be bounded between 0 and 1
                prec.b.prime = 10000, prec.lag.b.prime = 10000,
                prec.c.prime = max(working$precip), prec.lag.c.prime =max(working$precip),

                prec.b = 10000, prec.lag.b = 10000,
                prec.c = max(working$precip_diff), prec.lag.c = max(working$precip_diff) 	   )

#par$n.c <- 15
#par$n.b <- 20
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


model_5BN_full_azure_var_comp_results <- anneal(model=model_5BN_vary_comp, par=par, var=var, source_data=working, par_lo=par_lo, par_hi=par_hi,
			  pdf=linear_dnorm_with_intercept, dep_var="incr_mm", max_iter=iterations)

setwd(out_dir)

write_results(model_5BN_full_azure_var_comp_results, file = paste(spp_list[index],"Model 5BN Full Azure Vary CompetitionResults.txt", sep=" "), data=F, print_whole_hist=F)
save(model_5BN_full_azure_var_comp_results, file=paste(spp_list[index],"Model 5BN Full Azure Vary Competition Results.Rdata",sep=" "))


