#==================================================================================================
# STAN version of age-structured state-space spawner-recruitment model with AR-1 process variation (Fleischman et al. CJFAS. 2013; Connors et al. MCF 2019)
# Creator: Brendan Connors, Fisheries and Oceans Canada 
# Date: 16.11.2020

# load required packages, settings, working directory =============================================
library(rstan)
library(here)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd(here())

# wrangle data =====================================================================================
Esc <- read.csv('DATA/Data_report_tables/Table_3_Expanded_estimates.csv')[,2];Esc <- Esc[!is.na(Esc)]  # read in escapement as median of expanded estimates, drop nas
Har <- read.csv('DATA/Data_report_tables/Table_4_harvest.csv')[,2] # removed commas from some entries in csv file before reading
S_cv <- read.csv('DATA/Data_report_tables/Table_3_Expanded_estimates.csv')[,10] # removed % from values in csv file before reading
H_cv <- rep(0.15,length(Har)) # placeholder for CV on harvest observations
a_min <- 3
a_max <- 6
nyrs <- length(Esc)
A <- a_max - a_min + 1
nRyrs <- nyrs + A - 1

# deal with SPAWNER age comps
raw_age_comps <- read.csv('DATA/Data_report_tables/Table_6_esc_age_comp.csv')
A_esc_comp <- matrix(NA, nrow = length(Esc), ncol = 4) # matrix to store age comps in

# years with data 
A_esc_comp[42:72,1] <- rowSums(raw_age_comps[30:60,2:3]) # combine age 2 and 3 into age 3 for years with data
A_esc_comp[13:41,1] <- raw_age_comps[1:29,8] # assign "other" ages comps to 3 yr olds for years with data
A_esc_comp[13:72,2] <- raw_age_comps[1:60,4] # age 4 for years with data
A_esc_comp[13:72,3] <- raw_age_comps[1:60,5] # age 5 for years with data
A_esc_comp[13:72,4] <- raw_age_comps[1:60,6] # age 6 for years with data, replaced nas with 0 in years with data for other age comps

A_esc_comp <- as.data.frame(A_esc_comp)

# which years are missing age comps?
missing <- which(is.na(A_esc_comp[,1]))

# assign average age comps to missing years
for(i in missing){
  A_esc_comp[i,] <- colMeans(A_esc_comp[30:60,],na.rm=T)
}

# age comp sample sizes
#A_esc_comp[42:72,5] <- raw_age_comps$Sample.size[30:60] # assume ESS equal to true sample sizes for years with data
A_esc_comp[42:72,5] <- 60 # assume ESS equal to 100 for years with data
A_esc_comp[1:41,5] <- 25 # assume ESS of 50 for years without age comp sample size data

# write to file in case needed elsewhere
A_comp_print <- cbind(seq(1948,2019),A_esc_comp)
                      
colnames(A_esc_comp_print) <-  c("year","3_yr_olds","4_yr_olds","5_yr_olds","6_yr_olds","ESS")
write.csv(A_esc_comp_print,
          file = "DATA/derived_esc_age_comps.csv",
          row.names = F)

# age observations in integers, sum of ages by year is assumed to be equivalent to ESS
S_comps <- round(A_esc_comp[,1:4]*A_esc_comp[,5])


# deal with HARVEST spawner age comps
raw_har_age_comps <- read.csv('DATA/Data_report_tables/Table_5_catch_age_comp.csv')
A_har_comp <- matrix(NA, nrow = length(Esc), ncol = 4) # matrix to store age comps in

# years with data 
A_har_comp[1:41,1] <- raw_har_age_comps[1:41,6]*0.5 # assign half of "other" ages comps to 3 yr olds for years with data
A_har_comp[1:41,4] <- raw_har_age_comps[1:41,6]*0.5 # assign half of "other" ages comps to 3 yr olds for years with data
A_har_comp[42:72,1] <- raw_har_age_comps[42:72,2] # age 3 for years with data
A_har_comp[42:72,4] <- raw_har_age_comps[42:72,5] # age 6 for years with data
A_har_comp[1:72,2] <- raw_har_age_comps[1:72,3] # age 4 for years with data
A_har_comp[1:72,3] <- raw_har_age_comps[1:72,4] # age 5 for years with data

A_har_comp <- as.data.frame(A_har_comp)

# harvest age comp sample sizes
A_har_comp[1:72,5] <- 100 # assume ESS equal to 100 for all years with data

# write to file in case needed elsewhere
A_har_comp_print <- cbind(seq(1948,2019),A_har_comp)

colnames(A_har_comp_print) <-  c("year","3_yr_olds","4_yr_olds","5_yr_olds","6_yr_olds","ESS")
write.csv(A_har_comp_print,
          file = "DATA/derived_harvest_age_comps.csv",
          row.names = F)

# age observations in integers, sum of ages by year is assumed to be equivalent to ESS
H_comps <- round(A_har_comp[,1:4]*A_har_comp[,5])

# Initialize Model Parameters ======================================================================

simdata <- data.frame(id = rep(1:nRyrs, each = 1),
                      x = runif(1,min=0.01, max=0.99),
                      y = runif(1,min=0.01, max=0.99),
                      z = runif(1,min=0.01, max=0.99),
                      zz = runif(1,min=0.01, max=0.99))

init_fn <- function(chain_id=1) {
  list(
    "lnR"=abs(rnorm(nRyrs, mean=0, sd=5)),
    "lnalpha"=abs(rnorm(1, mean=0, sd=1)),
    "beta"=runif(1, min=0.01, max=9.99),
    "sigma_R"=abs(rnorm(1, mean=0, sd=1)),
    "sigma_R0"=abs(rnorm(1, mean=0, sd=1)),
    "phi"=runif(1, min=-1, max=1),
    "lnresid_0"=runif(1, min=-1, max=1),
    "mean_ln_R0"=abs(rnorm(1, mean=0, sd=1)),
    "U"=runif(nyrs, min=0.01, max=0.99),
    "g"=plyr::daply(simdata %>% mutate(id = as.integer(id)), "id",
              function(df) df[1,c("x", "y", "z", "zz")]) %>% as.numeric %>% matrix(ncol=4)
  )
}
# Initial List of Lists for Multiple Chains
init_ll <- lapply(1:4, function(id) init_fn(chain_id = id))

# Run Stan Model ===================================================================================

# data
stan.data <- list("nyrs" = nyrs,
                  "a_min" = a_min,
                  "a_max" = a_max,
                  "A" = A,
                  "nRyrs" = nyrs + A - 1,
                  "H_comps" = as.matrix(H_comps),
                  "S_comps" = as.matrix(S_comps),
                  "S_obs" = Esc/100000,
                  "H_obs" = Har/100000,
                  "S_cv" = S_cv,
                  "H_cv" = H_cv)

# fit 
stan.fit <- stan(file = "SRA/SSSR_AR1_v2.stan",
                 model_name = "Rivers-SSSR-AR1-v2",
                 data = stan.data,
                 chains = 4,
                 iter = 2000,
                 seed = 42,
                 init = init_ll,
                 control = list(adapt_delta = 0.99, max_treedepth = 20))

  
# Summarize model ===================================================================================
#shinystan::launch_shinystan(stan.fit)  

#pairs(stan.fit, pars = c("lnalpha", "beta", "sigma_R", "sigma_R0", "phi", "mean_ln_R0","lnR[1]", "lnR[2]", "log_posterior"))

#summary(stan.fit, pars= c("lnalpha", "beta", "sigma_R", "sigma_R0", "phi", "mean_ln_R0","S_max", "S_eq", "S_msy", "U_msy"))$summary[,4:10]  
  
# save fitted model object 
saveRDS(stan.fit, file="./SRA/output/stan.fit.v2.rds")

# save model data object 
saveRDS(stan.data, file="./SRA/output/stan.data.v2.rds")

