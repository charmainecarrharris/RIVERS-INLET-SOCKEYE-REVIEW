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
H_cv <- rep(0.05,length(Har)) # placeholder for CV on harvest observations
a_min <- 3
a_max <- 6
nyrs <- length(Esc)
A <- a_max - a_min + 1
nRyrs <- nyrs + A - 1

# deal with age comps
raw_age_comps <- read.csv('DATA/Data_report_tables/Table_6_esc_age_comp.csv')
A_comp <- matrix(NA, nrow = length(Esc), ncol = 5) # matrix to store age comps in

# years with data 
A_comp[42:72,1] <- rowSums(raw_age_comps[30:60,2:3]) # combine age 2 and 3 into age 3 for years with data
A_comp[13:41,1] <- raw_age_comps[1:29,8] # assign "other" ages comps to 3 yr olds for years with data
A_comp[13:72,2] <- raw_age_comps[1:60,4] # age 4 for years with data
A_comp[13:72,3] <- raw_age_comps[1:60,5] # age 5 for years with data
A_comp[13:72,4] <- raw_age_comps[1:60,6] # age 6 for years with data, replaced nas with 0 in years with data for other age comps

# years without data
A_comp <- as.data.frame(A_comp)
comp_means <- colMeans(A_comp[30:60,],na.rm=T)
A_comp <- A_comp %>% replace_na(list(V1=comp_means[1],
                           V2=comp_means[2],
                           V3=comp_means[3],
                           V4=comp_means[4]))

# age comp sample sizes
#A_comp[42:72,5] <- raw_age_comps$Sample.size[30:60] # assume ESS equal to true sample sizes for years with data
A_comp[42:72,5] <- 100 # assume ESS equal to 100 for years with data
A_comp[1:41,5] <- 50 # assume ESS of 50 for years without age comp sample size data

# write to file in case needed elsewhere
A_comp_print <- cbind(seq(1948,2019),A_comp)
                      
colnames(A_comp_print) <-  c("year","3_yr_olds","4_yr_olds","5_yr_olds","6_yr_olds","ESS")
write.csv(A_comp_print,
          file = "DATA/derived_age_comps.csv",
          row.names = F)

# age observations in integers, sum of ages by year is assumed to be equivalent to ESS
A_obs <- round(A_comp[,1:4]*A_comp[,5])

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
                  "A_obs" = as.matrix(A_obs),
                  "S_obs" = Esc/100000,
                  "H_obs" = Har/100000,
                  "S_cv" = S_cv,
                  "H_cv" = H_cv)

  
# fit 
stan.fit <- stan(file = "SRA/SSSR_AR1.stan",
                 model_name = "Rivers-SSSR-AR1",
                 data = stan.data,
                 chains = 4,
                 iter = 1000,
                 seed = 42,
                 init = init_ll)
  
# Summarize model ===================================================================================
shinystan::launch_shinystan(stan.fit)  

pairs(stan.fit, pars = c("lnalpha", "beta", "sigma_R", "sigma_R0", "phi", "mean_ln_R0","lnR[1]", "lnR[2]"))

summary(stan.fit, pars= c("lnalpha", "beta", "sigma_R", "sigma_R0", "phi", "mean_ln_R0","S_max", "S_eq", "S_msy", "U_msy"))$summary[,4:10]  
  

# Spawner-recruit figure =============================================================================

# brood table
fy <- 1948 # first year of spawners
no_by <- 69 # no. of brood years
byrs <- seq(fy,fy+no_by-1)

spwn <- summary(stan.fit, pars= c("S"),probs = c(0.025, 0.50, 0.975))$summary[,4:6]  
rec <- summary(stan.fit, pars= c("R"),probs = c(0.025, 0.50, 0.975))$summary[,4:6] 

bt <- matrix(NA,no_by,7)
bt[,1] <- byrs
bt[,2:4] <- spwn[1:no_by,]
bt[,5:7] <- rec[(a_max + 1):nRyrs,]
colnames(bt) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
bt <- as.data.frame(bt)

# SR relationship
spw <- seq(0,max(spwn[,2]*1.5),length.out=100)
SR_pred <- matrix(NA,100,10000)
posteriors <- as.matrix(stan.fit)

for(i in 1:10000){
  r <- sample(seq(1,1000),1,replace=T)
  a <- posteriors[r,"lnalpha"]
  b <- posteriors[r,"beta"]
  SR_pred[,i] <- (exp(a)*spw*exp(-b*spw))
}

SR_pred <- cbind(spw,t(apply(SR_pred,c(1),quantile,probs=c(0.05,0.5,0.95),na.rm=T)))
colnames(SR_pred) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr")
SR_pred <- as.data.frame(SR_pred)

ggplot() +
  geom_ribbon(data = SR_pred, aes(x = Spawn, ymin = Rec_lwr, ymax = Rec_upr), 
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_line(data = SR_pred, aes(x = Spawn, y = Rec_med), color="black", size = 1) +
  geom_errorbar(data = bt, aes(x= S_med, y = R_med, ymin = R_lwr, ymax = R_upr), width=0,colour="light grey", width=0.2, size=0.3) +
  geom_errorbarh(data = bt, aes(x= S_med, y = R_med, xmin = S_lwr, xmax = S_upr), 
                 height=0,colour = "light grey", width=0.2, size = 0.3) + 
  geom_point(data = bt, aes(x = S_med, y = R_med, color=BroodYear, width=0.9), size=2)+
  coord_cartesian(xlim=c(0, 33), ylim=c(0,50))+
  scale_colour_viridis_c()+
  xlab("Spawners (100,000s)") +
  ylab("Recruits (100,000s)") +
  theme_bw() +
  theme(strip.text.x = element_text(size=8),
        axis.title = element_text(size=10),
        axis.text = element_text(size=7))+
  geom_abline(intercept = 0, slope = 1)
ggsave("./FIGURES/SS-SRA_spawner-recruit.jpg",height = 4.5, width = 6)

# AR-1 prod figure =============================================================================

resid <- summary(stan.fit, pars= c("lnalpha_time"),probs = c(0.025, 0.50, 0.975))$summary[(a_max + 1):nRyrs,4:6]
resid_exp <- exp(summary(stan.fit, pars= c("lnalpha_time"),probs = c(0.025, 0.50, 0.975))$summary[(a_max + 1):nRyrs,4:6])

alpha_resid <- matrix(NA,no_by,4)
alpha_resid[,1] <- byrs
alpha_resid[,2:4] <- resid[1:no_by,]
colnames(alpha_resid) <- c("BroodYear","lnResid_lwr","lnResid_med","lnResid_upr")
alpha_resid <- as.data.frame(alpha_resid)


ggplot(alpha_resid, aes(x=BroodYear, y = lnResid_med ), show.legend = F) +
  geom_point(size=1,show.legend = F)+
  geom_line(show.legend = F) + 
  geom_ribbon(aes(ymin = lnResid_lwr, ymax = lnResid_upr), show.legend = F, alpha=0.4) +
  xlab("Brood year") +
  ylab("Productivity index") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, lty = "dotted") +
  theme_bw()
ggsave("./FIGURES/SS-SRA_prod_index.jpg",height = 3, width = 8.5)


  
