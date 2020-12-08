#==================================================================================================
# Dynamic linear models (e.g., fit via Kalman filter) for linearized stock recruit curves for River Inlet sockeye. 
# Depending on formulation this allows the parameters of the linear model to vary over time 
# Creator: Brendan Connors, Fisheries and Oceans Canada, adapted from code provide by Cody Szuwalski and Brigitte Dorner 
# Date: 16.11.2020

# load required packages, settings, working directory =============================================
library(dlm)
library(viridis)
library(here)
library(tidyverse)

setwd(here())

source("Kalman_code.R")

# create brood table ===============================================================================
Esc <- read.csv('DATA/Data_report_tables/Table_3_Expanded_estimates.csv')[,2];Esc <- Esc[!is.na(Esc)]
Har <- read.csv('DATA/Data_report_tables/Table_4_harvest.csv')[,2] # removed commas from some entries in csv file before reading
Age<-read.csv("DATA/derived_age_comps.csv")

brood <- matrix(NA,72,8)
colnames(brood) <- c("year",
                     "spawners",
                     "tot_returns",
                     "age3_returns",
                     "age4_returns",
                     "age5_returns",
                     "age6_returns",
                     "recruits")

brood[,1]<- Age$year # brood year
brood[,2]<- Esc # spawners
brood[,3]<- Esc+Har # returns
brood[,4]<- brood[,3]* Age[,2]  # age 3 returns
brood[,5]<- brood[,3]* Age[,3] # age 4 returns
brood[,6]<- brood[,3]* Age[,4] # age 5 returns
brood[,7]<- brood[,3]* Age[,5] # age 6 returns
for(i in 1:66){
  brood[i,8] <- brood[i+3,4] + brood[i+4,5] + brood[i+5,6] + brood[i+6,7] 
}

btf <- as.data.frame(brood[complete.cases(brood), ]) # drop incomplete brood years

# pre process for DLM ===============================================================================

recr		<-btf$recruits/100000
Tssb		<-btf$spawners/100000
BY <- btf$year

lnRS	<-log(btf$recruits/btf$spawners)
alpha		<-NULL
beta		<-NULL
plotIn	<-lnRS[-length(lnRS)]

# spawner recruit data
mod 		<- dlmModReg(Tssb) # this specifies a linear model

# this sets up output for models that consider all potential formulations of the variance structure
lls 		<- numeric(4)
dlm_out <- list()
AICc    <-rep(0,4)

# this determines the number of parameters used in the AICc calculation
dlmPars	<-c(3,4,4,5)

#
for(i in 1:4)
{
  # this changes the variance structure for the linear regression model built above
  build_mod <- function(parm) 
  {
    mod$V <- exp(parm[1]) 	
    if(i==2){mod$W[1,1]=exp(parm[2]); mod$W[2,2]=0}
    if(i==3){mod$W[1,1]=0; mod$W[2,2]=exp(parm[2])}
    if(i==4){mod$W[1,1]=exp(parm[2]); mod$W[2,2]=exp(parm[3])}
    return(mod)
  }
  
  # this does maximum likelihood optimization of the variance
  dlm_out[[i]] <- dlmMLE(y=lnRS, build=build_mod, parm=c(-.1,-.1,-.1), method="Nelder-Mead")
  lls[i] 	 <- dlm_out[[i]]$value
  
  # specifies the model based on 
  dlmMod	<-build_mod(dlm_out[[i]]$par)
  
  # applies Kalman filter 
  outsFilter	<-dlmFilter(y=lnRS,mod=dlmMod)
  filtered	<-outsFilter$m[-1,1]+Tssb*outsFilter$m[-1,2] # estimated log(r/s)
  
  # backward recurvise smoothing 
  outsSmooth	<-dlmSmooth(outsFilter)
  smoothed	<-outsSmooth$s[-1,1]+Tssb*outsSmooth$s[-1,2] # estimated log(r/s)
  
  # collecting the parameters
  alpha 	<- cbind(alpha,outsSmooth$s[-1,1,drop=FALSE])
  beta  	<- cbind(beta,outsSmooth$s[-1,2,drop=FALSE])
  
  # calculation of AICc for model comparison (I'm a little suspicious of this, but that's what we did)
  AICc[i]	<- 2*lls[i] + 2*dlmPars[i] +(2*dlmPars[i]*(dlmPars[i]+1)/(length(recr)-dlmPars[i]-1))
}

# print AICc for each model
AICc


# plot estimated stock recruit curves
# each line in top row represents the yearly 'realized' SR curve based on 
# estimated parameters
# alpha and beta rows are the estimates of the parameters from the dlm

jpeg("./DLM/figures/Rivers_DLM_comparision.jpg",width=7, height=5, units="in",res=800)

par(mfcol=c(3,4),mar=c(.1,.1,.1,.1),oma=c(4,5,5,4),las=1)
plotSB<-seq(0,max(Tssb),1)
time_col <- viridis(nrow(alpha))

legend_txt<-c("Static","Alpha varies","Beta varies","Both vary")
for(x in 1:4)
{
  plot(recr~Tssb,xlim=c(0,max(Tssb)),ylim=c(-0.5,max(recr)),xaxt='n',yaxt='n',pch=16)
  for(y in 1:nrow(alpha))
    lines(plotSB*exp(alpha[y,x]+plotSB*beta[y,x])~plotSB,col=time_col[y])
  points(recr~Tssb,xlim=c(0,max(Tssb)),ylim=c(-0.5,max(recr)),xaxt='n',pch=16)
  legend('topleft',bty='n',legend_txt[x])
  box(col="grey")
  if(x==1)
  {
    axis(side=2,las=1,col="grey")
    axis(side=3,las=1,col="grey")
    mtext(side=2,"Recruits (^5)",las=0,line=3.5)
    mtext(side=3,"Spawners (^5)",line=2.25)
  }
  
  plot(alpha[,x]~BY,yaxt='n',xaxt='n',type='l',ylim=c(min(alpha),max(alpha)))
  box(col="grey")
  if(x==1)
  {
    axis(side=2,las=1, col="grey")
    mtext(side=2,"Alpha",las=0,line=3.5)
  }

  plot(beta[,x]~ BY,yaxt='n',xaxt='n',type='l',ylim=c(min(beta),max(beta)))
  box(col="grey")
  if(x==1)
  {
    axis(side=2,las=1, col="grey", at=c(-0.2,-0.15,-0.1,-0.05))
    axis(side=1,las=1,col="grey")
    mtext(side=2,"Beta",las=0,line=3.5)
  }
  }

dev.off()

# Kalman filter al-la Peterman et al. CJFAS 2003 ===============================================================================

# fit model
data <- cbind(btf,lnRS); colnames(data)[2] <- "S"; colnames(data)[1] <- "BY" # add log(R/S) to brood table and change spawners column name to S for Kalman filter function
kalman_fit <- run.kalman.Ricker(data)

# plot estimate of time varying productivity
ggplot(kalman_fit$df, aes(x=BY, y = a.smooth ), show.legend = F) +
  geom_point(size=1,show.legend = F)+
  geom_line(show.legend = F) + 
  geom_ribbon(aes(ymin = a.smooth-(a.smooth.var*2), ymax = a.smooth+(a.smooth.var*2)), show.legend = F, alpha=0.4) +
  xlab("Brood year") +
  ylab("Productivity index") +
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 0 ,col="dark grey", lty=2) +
  theme_bw()+
  coord_cartesian(ylim=c(-2.25,3)) 

ggsave("./DLM/figures/kalman_prod_index.jpg",height = 4, width = 6.5)
