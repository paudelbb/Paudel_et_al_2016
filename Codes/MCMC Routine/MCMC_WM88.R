# ===========================================================================================================
# Figure 4C: MCMC sampling for A375 model fitting
# ===========================================================================================================
library(gplots)
require(ggplot2)
require(FME)
require(deSolve)
setwd("~/Paudel_et_al_2016/SourceData/Model_Fitting") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData/Model_Fitting"))), "Figures") 
fig = "Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)
# ===========================================================================================================
#                                     Read experimental data
# ===========================================================================================================

data <- read.csv("Population_Response_CellLines_post_equilibration.csv", header=T, sep=",")
cell <- "WM88"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc & data$Time<=350)
d_exp <- data.matrix(s1[,c("Time", "nl2")]);
colnames(d_exp) <- c("time", "nl2")
d_exp1 <- data.frame(time=seq(0,500), nl2=c(-1000))
d_exp1 <- data.matrix(d_exp1)

# ===========================================================================================================
#                   Ordinary differential equation models and initial parameters
# ===========================================================================================================

solvePop <- function(pars, times, Tot = 10000){
  derivs <- function(t, state, pars){
    with(as.list(c(state, pars)),{
      dA <- (-0.055 - k_ab)*A + k_ba * B
      dB <- (0.000 - k_ba - k_bc)*B + k_ab * A + k_cb * C
      dC <- (0.015 - k_cb)*C + k_bc * B
      return(list(c(dA, dB, dC), tot = A0+B0, l2=log2(A+B+C), nl2 = (log2(A+B+C)-log2(Tot))))
    })
  }
  state <- with(as.list(pars), c(A = Tot * A0, B = Tot * B0, C= Tot * (1 - A0 - B0)))
  output = ode(y=state, times=times, func=derivs, parm=pars, method="lsoda", hmax=1, jacfunc=mod)
  print(diagnostics(output))
  return(output)
}

# ===========================================================================================================
#                       Cost function to find the residuals between the model and the data
# ===========================================================================================================
pars <- c(k_ab    = 1/100,
          k_ba    = 1/200,
          k_cb    = 1/100,
          k_bc    = 1/200,
          A0 = 1/3, 
          B0 = 1/3)
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
    dB <- (0.000 - k_ba - k_bc)*B + k_ab * A + k_cb * C
    dC <- (0.015 - k_cb)*C + k_bc * B
    return(list(c(dA, dB, dC)))
  })
}

# ===========================================================================================================
#                                   Check to see if ODE solver works
# ===========================================================================================================

out <- solvePop(pars, times)
plot(out)

# ===========================================================================================================
#             Model Fitting algorithm--uses Lavenberg-Marq method from R-package FME (function modFit)
# ===========================================================================================================

Fit <- modFit(p= c(pars), f = Objective, lower=0, upper=c(0.06, 0.06, 0.06, 0.06, 1, 1),
              hessian = T)
pars[c("k_ab", "k_ba", "k_cb", "k_bc","A0", "B0")] <- Fit$par
summary(Fit)

# ===========================================================================================================
#                         MCMC algorithms--uses modMCMC from FME package
# ===========================================================================================================
output <- paste0(substr(getwd(), 1, (nchar(getwd()))))
setwd(output)
set.seed(12345)
# Define Prior distributions
gprior <- function(pars){
  for(i in 1:length(pars)){
    return(-2*log(dnorm(pars, pars[i], 1)))
  }
}
# Number of iterations to run for each MCMC
iter = 150000
# Var0 to scale each variable independently, obtained from modFit results.
var0 = Fit$var_ms_unweighted
cov0 = summary(Fit)$cov.scaled * 2.4^2/6
#Parameter covariance is defined as 1/(0.5*H)*(sum of squared residuals)
cov = (1/(0.5*Fit$hessian))*Fit$ssr

# ===========================================================================================================
# Define prior distributions for each parameter
# ===========================================================================================================
n = 75000
p1 = rnorm(n, mean = pars[1], sd = var0)
p2 = rnorm(n, mean = pars[2], sd = var0)
p3 = rnorm(n, mean = pars[3], sd = var0)
p4 = rnorm(n, mean = pars[4], sd = var0)
prior = data.frame(k_ab = p1, k_ba = p2, k_cb = p3, k_bc = p4)
write.csv(prior, file=paste0(output, "/prior", cell, ".csv"))


# ===========================================================================================================
# Run MCMC algorithms
# ===========================================================================================================

print(system.time(MCMC1 <- modMCMC(f = Objective, p = Fit$par, niter = iter, lower=0, 
                                   updatecov = 100,  var0 = var0, ntrydr = 2, prior = gprior, 
                                   upper=c(0.06, 0.06, 0.06, 0.06, 1, 1))))

plot(MCMC1$pars[,5]~MCMC1$pars[,6], xlab="B0", ylab="A0", xlim=c(0,1), ylim=c(0,1))
lines(x = c(1,-0.0), y = c(-0.0,1), col="red", lwd=3)
plot(MCMC1, Full = T)
pairs(MCMC1, nsample = 100)# Plots the distributions for each parameter

save(MCMC1, file=paste0(output, "/MCMC1_", cell, "_", iter,"_Run.RData"))
# ===========================================================================================================