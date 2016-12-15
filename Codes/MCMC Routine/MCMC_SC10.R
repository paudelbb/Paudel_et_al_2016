# ===========================================================================================================
# Figure 4: SKMEL5 Subclones Mix Model Fitting and MCMC Algorithms
# ===========================================================================================================
library(gplots)
require(ggplot2)
require(FME)
require(deSolve)
library(devtools)
library(ggbiplot)

setwd("~/Paudel_et_al_2016/SourceData/Model_Fitting") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData/Model_Fitting"))), "Figures") 
fig = "Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)
# ===========================================================================================================
#                                         Read experimental data
# ===========================================================================================================
data <- read.csv("SKMEL5_Subclones_Population_Level.csv", header=T, sep=",")
cell <- "SKMEL5 Subclone10"; cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
d_exp <- data.matrix(s1[,c("Time", "nl2")]);
colnames(d_exp) <- c("time", "nl2")
d_exp1 <- data.frame(time=seq(0,500), nl2=c(-1000))

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

# ===========================================================================================================
#               Cost function to find the residuals between the model and the data
# ===========================================================================================================
pars <- c(k_ab    = 1/100,
          k_ba    = 1/200,
          k_cb    = 1/100,
          k_bc    = 1/200,
          A0 = 1/4, 
          B0 = 1/4)
times <- unique(s1$Time)
times <- times[order(times)]
# Objective function to calculate the cost--using modCost of R-package-FME
Objective <- function(x){
  pars <- pars
  pars[names(x)] <- x
  out <- solvePop(pars, times)
  out <- out[,c("time", "nl2")]
  c <- pars["A0"] + pars["B0"]
  print(c)
  if(c <= 1){
    cost1 <- modCost(model = out, obs = d_exp)
    print(cost1$model)
    return(cost1)
  }
  else {
    cost2 = modCost(model = out, obs = d_exp1)
    print(cost2$model)
    return(cost2)
  }
}
# Function to calculate Jacobian of the ODEs.
mod <- function(t, state, pars){
  with(as.list(c(state, pars)),{
    dA <- (-0.055 - k_ab)*A + k_ba * B
    dB <- ( 0.000 - k_ba - k_bc)*B + k_ab * A + k_cb * C
    dC <- ( 0.015 - k_cb)*C + k_bc * B
    return(list(c(dA, dB, dC)))
  })
}
# ===========================================================================================================
#                               Check to see if ODE solver works
# ===========================================================================================================

out <- solvePop(pars, times)
plot(out)

# ===========================================================================================================
#                         Model Fitting algorithm--uses Lavenberg-Marq method
# ===========================================================================================================

Fit <- modFit(p= c(pars[c("k_ab", "k_ba", "k_cb", "k_bc", "A0", "B0")]), f = Objective, 
              lower=0, upper=c(0.06, 0.06, 0.06, 0.06, 1, 1),
              hessian = T)
pars[c("k_ab", "k_ba", "k_cb", "k_bc","A0", "B0")] <- c(Fit$par)
summary(Fit)

# ===========================================================================================================
#  MCMC algorithms--uses modMCMC from FME package
# Using the distributions for Subclones Mix as a prior distribution
# ===========================================================================================================
output <- paste0(substr(getwd(), 1, (nchar(getwd()))))
setwd(output)
set.seed(12345)
files1 <- grep(".RData", dir(), value=T); 
chain1 <- load("MCMC_Subclones_Mix_150000_Run.RData")
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
dm <- data.frame(MCMC1$pars[,c(1,2,3,4)])

# ===========================================================================================================
#  MCMC algorithms--uses modMCMC from FME package
# ===========================================================================================================
pars <- Fit$par
names(pars) <- names(Fit$par)
set.seed(12345) # Define Prior distributions
gprior <- function(pars){
  for(i in 1:length(pars)){
    return(-2*log(dnorm(pars, pars[i], 1)))
  }
}# Number of iterations to run for each MCMC
iter = 150000
var0 = Fit$var_ms_unweighted # Var0 to scale each variable independently, obtained from modFit results.
cov0 = summary(Fit)$cov.scaled * 2.4^2/6 #Parameter covariance is defined as 1/(0.5*H)*(sum of squared residuals)
cov = (1/(0.5*Fit$hessian))*Fit$ssr
# ============================================================================================
# Run MCMC algorithms
# ============================================================================================
print(system.time(MCMC1 <- modMCMC(f = Objective, p = c(pars), niter = iter, 
                                   updatecov = 100,  var0 = var0, ntrydr = 2, prior = gprior, 
                                   lower=c(min(dm$k_ab), min(dm$k_ba), min(dm$k_cb), min(dm$k_bc), 0.0, 0.0), 
                                   upper=c(max(dm$k_ab), max(dm$k_ba), max(dm$k_cb), max(dm$k_bc), 1, 1))))

plot(MCMC1$pars[,5]~MCMC1$pars[,6], xlab="B0", ylab="A0", xlim=c(0,1), ylim=c(0,1))
lines(x = c(1,-0.0), y = c(-0.0,1), col="red", lwd=3)

plot(MCMC1, Full = T)
pairs(MCMC1, nsample = 100)# Plots the distributions for each parameter
save(MCMC1, file=paste0(output, "/MCMC1_", cell, "_", iter,"_Run.RData"))


