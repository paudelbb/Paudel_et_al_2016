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
data <- read.csv("Population_Response_CellLines_post_equilibration.csv", header=T, sep=",")
cell <- "SKMEL28"
cnc <- 16
s1 <- subset(data, data$CellLine==cell & data$conc==cnc & data$Time<=300)
d_exp <- data.matrix(s1[,c("Time", "nl2")]);
colnames(d_exp) <- c("time", "nl2")
d_exp1 <- data.frame(time=seq(0,500), nl2=c(-1000))
times <- unique(s1$Time)
times <- times[order(times)]
times <- seq(1:max(s1$Time)); times <- c(0, times)
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
chain1 <- load("MCMC1_SKMEL28_150000_Run.RData")
dim(MCMC1$pars)
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
rates <- data.frame(MCMC1$pars)
nd <- data.frame(rates[sample(nrow(rates), 1000, replace = F), ])
bestFit <- data.frame(solvePop(MCMC1$bestpar, times))

# ===========================================================================================================
# Simulation of MCMC sets
# ===========================================================================================================

dn <- data.frame()
for(i in 1:1000){
  pars <- nd[i,]
  out <- as.data.frame(solvePop(pars, times))
  out <- out[,c("time", "nl2")]
  out$id <- i
  dn <- rbind(dn, out)
}
# 
cl1 <- rainbow(length(unique(dn$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~time, data=dn, type="n", ylim=c(-2,1), xlim=c(0,300), ylab="", xlab="")
for(w in 1:length(unique(dn$id))){
  temp <- subset(dn, dn$id==w)
  lines(temp$time, temp$nl2, col=cl1[w])
}
lines(bestFit$nl2~bestFit$time, col="black", lwd=2)

# Write MCMC ensembles into a csv file
write.csv(dn, file=paste0(getwd(), "/MCMC_ensembles_", cell,".csv"))
# ===========================================================================================================
