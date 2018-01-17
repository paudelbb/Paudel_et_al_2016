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
require(Hmisc)
# ===========================================================================================================
setwd("~/Paudel_et_al_2016/SourceData/Model_Fitting") # Set the working directory;
# Set the output for plots
output <- "/Users/paudelbb/Dropbox/Paudel_et_al_2017_idle_paper" 
fig = "Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)

# ===========================================================================================================
# Calculating the energy of activation from the rate constants
# ===========================================================================================================
# Arrhenius Equation: k = A*exp(-Ea/RT)
# Solving this, we get: Ea = -K*ln(k); 
# Sinusoidal function to model the attractor states with 2*Amplitude == Energy barrier from rate constants
# ===========================================================================================================
# Get the MCMC estimated parameters
# ===========================================================================================================
# ===========================================================================================================
# ===========================================================================================================
# PLot the energy landscape based on the model fitted transition rate constants
# ===========================================================================================================
# ===========================================================================================================
#                       Ordinary differential equation models and initial parameters
# ===========================================================================================================

solvePop <- function(pars, times, Tot = 10000){
  derivs <- function(t, state, pars){
    with(as.list(c(state, pars)),{
      dA <- (-0.055 - k_ab)*A + k_ba * B
      dB <- ( 0.000 - k_ba - k_bc)*B + k_ab * A + k_cb * C
      dC <- ( 0.015 - k_cb)*C + k_bc * B
      return(list(c(dA, dB, dC), tot = A0+B0, l2=log2(A+B+C), nl2 = (log2(A+B+C)-log2(Tot))))
    })
  }
  state <- with(as.list(pars), c(A = Tot * A0, B = Tot * B0, C= Tot * (1 - A0 - B0)))
  output = ode(y=state, times=times, func=derivs, parm=pars, method="lsoda", hmax=1, jacfunc=mod)
  print(diagnostics(output))
  return(output)
}
# Function to calculate Jacobian of the ODEs.
mod <- function(t, state, pars){
  with(as.list(c(state, pars)),{
    dA <- (-0.055 - k_ab)*A + k_ba * B
    dB <- (0.000 - k_ba - k_bc)*B + k_ab * A + k_cb * C
    dC <- (0.015 - k_cb)*C + k_bc * B
    return(list(c(dA, dB, dC)))
  })
}
# ===========================================================================================================
# 1000 randomly selected MCMC ensembles
# ===========================================================================================================
set.seed(12345)
files <- grep(".RData", dir(), value=T)
chain1 <- load("MCMC1_A2058_150000_Run.RData")
dim(MCMC1$pars)
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
rates <- data.frame(MCMC1$pars[,c(1,2,3,4, 5,6)])
cell = "A2058"
nd <- data.frame(rates[sample(nrow(rates), 10, replace = F), ])
bestFit <- data.frame(solvePop(MCMC1$bestpar, times))

times = seq(0, 10000)
# ===========================================================================================================
# Simulation of MCMC sets
# ===========================================================================================================

dn <- data.frame()
for(i in 1:10){
  pars <- nd[i,]
  out <- as.data.frame(solvePop(pars, times))
  #out <- out[,c("time", "nl2")]
  out$id <- i
  dn <- rbind(dn, out)
}
# 
dn$cell <- cell
cl1 <- rainbow(length(unique(dn$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~time, data=dn, type="n",  xlim=c(0,10000), ylab="", xlab="")
for(w in 1:length(unique(dn$id))){
  temp <- subset(dn, dn$id==w)
  lines(temp$time, temp$nl2, col=cl1[w])
}
lines(bestFit$nl2~bestFit$time, col="black", lwd=2)

# Write MCMC ensembles into a csv file

# ===========================================================================================================
