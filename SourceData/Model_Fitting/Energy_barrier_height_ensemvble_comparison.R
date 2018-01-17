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
setwd("/Users/paudelbb/Dropbox/Paudel_et_al_2017_idle_paper/Barrier Height") # Set the working directory;

# ===========================================================================================================
# Calculating the energy of activation from the rate constants
# ===========================================================================================================
# Arrhenius Equation: k = A*exp(-Ea/RT)
# Solving this, we get: Ea = -K*ln(k); 
# ===========================================================================================================
# Get the MCMC estimated parameters
# ===========================================================================================================
files = dir()
dn = data.frame()
for(i in 1:length(files)){
  temp = read.csv(files[i])
  dn = rbind(dn, temp)
}
melanoma = dn
# ===========================================================================================================
boxplot(AB~cell, data=melanoma, ylim=c(0,15), main=paste0("AB"))
boxplot(BA~cell, data=melanoma, ylim=c(0,15), main=paste0("BA"))
boxplot(BC~cell, data=melanoma, ylim=c(0,15), main=paste0("BC"))
boxplot(CB~cell, data=melanoma, ylim=c(0,15), main=paste0("CB"))
# ===========================================================================================================


