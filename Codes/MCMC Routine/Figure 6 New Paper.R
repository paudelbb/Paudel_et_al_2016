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
fig = "Supp.Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)
source("~//Paudel_et_al_2016//Codes//summarySE.R") # Summarize function to get the summary statistics

# ===========================================================================================================
# ===========================================================================================================
# ===========================================================================================================
#                       Ordinary differential equation models and initial parameters
# ===========================================================================================================
times = seq(0, 300)
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
chain1 <- load("MCMC1_SKMEL5 Subclone07_150000_Run.RData")
dim(MCMC1$pars)
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
rates <- data.frame(MCMC1$pars)
nd <- data.frame(rates[sample(nrow(rates), 1000, replace = F), ])
bestFit <- data.frame(solvePop(MCMC1$bestpar, times))
pars = MCMC1$bestpar
pars[5] = 1/3
pars[6] = 1/3
pars[2] = 0.056586924/9.3
# ===========================================================================================================
# Simulation of MCMC sets
# ===========================================================================================================

out <- as.data.frame(solvePop(pars, times))
dn = out
dn$cell = "model"
# ===========================================================================================================
dn$tot = dn$A + dn$B + dn$C

# ===========================================================================================================
# ===========================================================================================================
pdf(paste0(output, "/Graphical_Schematic_1.pdf"), width=6.2, height=3.5)
plot(A/1000~time, data=dn, ylim=c(-1.5, 12), col="darkorange", ylab="", xlab="", type = "l", yaxt="n", xaxt="n", lwd=3)
points(B/1000~time, data=dn, col="darkgreen", ylab="", xlab="", type = "l", xaxt="n", lwd=3)
points(C/1000~time, data=dn, col="blue", ylab="", xlab="", type = "l", xaxt="n", lwd=3)
dev.off()

pdf(paste0(output, "/Graphical_Schematic_2.pdf"), width=6.2, height=3.5)
plot(tot/1000~time, data=dn, ylim=c(8, 12), ylab="", xlab="", type="l", xaxt="n", yaxt="n", lwd=3)
dev.off()

