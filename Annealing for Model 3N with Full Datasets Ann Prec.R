#############################################################################
#
#	Central Europe Tree Growth Project
#
#       Annealing with Full Datasets
#
#       Model 3
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
  dat_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/"
  out_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/Output_azure/"
  code_dir <- "/Users/arnebuechling/Library/Mobile Documents/com~apple~CloudDocs/CULS/Projects/Carpathian_tree_growth/Azure_model_files/Code/Annealing_code/"
}


### Get model functions
setwd(code_dir)
source("Model Functions - EDF versions.R")

# focal species
spp_list <- c("PICABI","FAGSYL","ABIALB","ACEPSE")

index <- 2

### load working target data file
setwd(paste(dat_dir,spp_list[index],sep=""))

load(paste(spp_list[index],"Full Working Dataset with All Precip Variables Nitrogen and Mean Climate.Rdata"))


### specify climate variables - use annual precipitation for species comparisons

working$temp_k <- working$tave_ann_hydro_k
working$temp_lag1_k <- working$tave_ann_hydro_k_lag1

### best precip variables

# if (index %in% c(1,4))     # Picea and Acer use annual precip
# {  working$precip <- working$ppt_ann_hydro_mm
#    working$precip_lag1 <- working$ppt_ann_lag1_hydro_mm
# }
# 
# if (index == 2)      # Fagus uses water deficit
# {  working$precip <- working$hydro_yr_water_deficit_mm
#    working$precip_lag1 <- working$hydro_yr_water_deficit_lag1_mm
# }
# if (index == 3)      # Abies uses seasonal precip
# {  working$precip <- working$seas_prec_hydro_yr_mm
#    working$precip_lag1 <- working$seas_prec_hydro_yr_lag1_mm
# }


working$precip <- working$ppt_ann_hydro_mm
working$precip_lag1 <- working$ppt_ann_lag1_hydro_mm


### Set annealing parameters

iterations <- 25000

setwd(paste(dat_dir,spp_list[index],sep=""))

if (index == 4) { load(paste(spp_list[index],"Model 3 Results.Rdata")) } else
                { load(paste(spp_list[index],"Model 3 Full Annual Prec Results.Rdata")) }


var <- list(mean = "predicted", x = "incr_mm", siteplot = "stdcode", log=T)

# set parameter limits

if (index == 4) { par <- model_3_results$best_pars } else
                { par <- model_3_full_annual_prec_results$best_pars  }

par_lo <- list(PG = rep(0,num_site_plot),
		sizeX0 = 0, sizeXb = 0,
		ageX0 = 0,
		int = 0, sigma = 0,
                alpha = 0, beta = 0, C = 0, gamma = -2, D = 1, lambda = rep(0,nspp),
		
                temp.a = 0, temp.lag.a = 0,   # the climate "a" parameters should be bounded between 0 and 1
		temp.b = 1, temp.lag.b = 1,     # the climate "b" parameters are the variances of the gaussian - small to large
                temp.c = min(working$temp_k), temp.lag.c = min(working$temp_k),  # the "c" parameters are the mode of the gaussian

                prec.a = 0, prec.lag.a = 0,
                prec.b = 10, prec.lag.b = 10,
                prec.c = min(working$precip), prec.lag.c = min(working$precip)   )

par_hi <- list(PG = rep(250,num_site_plot),
		sizeX0 = 100, sizeXb = 1,
		ageX0 = 1,
                int = 100, sigma = 2,
        	alpha = 4, beta = 3, C = 1000, gamma = 2, D = 3, lambda = rep(1,nspp),

		temp.a = 1, temp.lag.a = 1,   # the climate "a" parameters should be bounded between 0 and 1
		temp.b = 5000, temp.lag.b = 5000,     # the climate "b" parameters are the variances of the gaussian - 0 to large
                temp.c = max(working$temp_k), temp.lag.c = max(working$temp_k),  # the "c" parameters are the mode of the gaussian

                prec.a = 1, prec.lag.a = 1,
                prec.b = 10000, prec.lag.b = 10000,
                prec.c = max(working$precip), prec.lag.c = max(working$precip)   )

par$prec.b <- 20
par$n.c <- 15
par$n.b <- 20
par_lo$n.c <- 0
par_lo$n.b <- 0
par_hi$n.c <- 40
par_hi$n.b <- 200


model_3N_full_annual_prec_results <- anneal(model=model_3N, par=par, var=var, source_data=working, par_lo=par_lo, par_hi=par_hi,
			  pdf=linear_dnorm_with_intercept, dep_var="incr_mm", max_iter=iterations)

setwd(out_dir)

write_results(model_3N_full_annual_prec_results, file = paste(spp_list[index],"Model 3N Full Annual Prec Results.txt", sep=" "), data=F, print_whole_hist=F)
save(model_3N_full_annual_prec_results, file=paste(spp_list[index],"Model 3N Full Annual Prec Results.Rdata",sep=" "))


###################################