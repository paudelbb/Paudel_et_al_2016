# ===========================================================================================================
# Inferring energy landscape from the rate constants;
# ===========================================================================================================
# Different rate constants that we have: kab, kba, kbc & kcb
# How likely it is to transition to A from A; kab
# How likely it is to transition to B from A; kba
# How likely it is to transition to C from B; kcb
# How likely it is to transition to B from C; kbc
library(gplots)
require(ggplot2)
require(FME)
require(deSolve)
library(devtools)
library(ggbiplot)
library(preprocessCore)
setwd("~/Paudel_et_al_2016/SourceData/Model_Fitting") # Set the working directory;
# Set the output for plots
output <- "/Users/paudelbb/Dropbox/Paudel_et_al_2017_idle_paper" 
fig = "Barrier Height" # Which figure is this in paper?
output = paste0(output, "/", fig)

# ===========================================================================================================
# Calculating the energy of activation from the rate constants
# ===========================================================================================================
# Arrhenius Equation: k = A*exp(-Ea/RT)
# Solving this, we get: Ea = -K*ln(k); 
# ===========================================================================================================
# Get the MCMC estimated parameters
# ===========================================================================================================
set.seed(12345)
files <- grep(".RData", dir(), value=T)
#chain1 <- load("MCMC_Subclones_Mix_150000_Run.RData")
#chain1 <- load("MCMC1_A375_150000_Run.RData" )
#chain1 <- load("MCMC1_SKMEL19_150000_Run.RData")
#chain1 <- load("MCMC1_SKMEL28_150000_Run.RData")
#chain1 <- load("MCMC1_WM164_1e_05_Run.RData")
#chain1 <- load("MCMC1_WM793_150000_Run.RData")
#chain1 <- load("MCMC1_WM88_150000_Run.RData")
chain1 <- load("MCMC1_A2058_150000_Run.RData")

dim(MCMC1$pars)
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
rates <- data.frame(MCMC1$pars[,c(1,2,3,4)])
cell = "A2058"

# ===========================================================================================================
# calculate the energy barrier between states based on Arrhenius equation 
# ===========================================================================================================
randomrows = sample(nrow(rates), 1000)
store = data.frame()
for(i in 1:length(randomrows)){
  pars <- rates[i,]
  dc = data.frame(AB = -log(pars$k_ab), BA = -log(pars$k_ba), BB = -log(1), BC = -log(pars$k_bc), CB = -log(pars$k_cb))
  store = rbind(store, dc)
}
store$cell = cell
# ===========================================================================================================
# ===========================================================================================================
write.csv(store, file=paste0(output, "/energy_barrier_", cell, ".csv"))


#boxplot(store, notch = T, ylim=c(0, 20), main=paste0(cell))


