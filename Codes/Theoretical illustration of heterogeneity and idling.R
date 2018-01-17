# ===========================================================================================================
# Figure 4C: SKMEL5 Parental and Subclones No Transitions
# ===========================================================================================================

library(gplots)
require(ggplot2)
require(FME)
require(deSolve)
library(devtools)
library(ggbiplot)
require(dplyr)

setwd("~/Paudel_et_al_2016/SourceData/Model_Fitting") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData/Model_Fitting"))), "Figures") 
fig = "Supp.Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)

# ===========================================================================================================
# ===========================================================================================================
#                       Ordinary differential equation models and initial parameters
# ===========================================================================================================
times = seq(0, 350, by=10)

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
# ===========================================================================================================

load("MCMC_Subclones_Mix_150000_Run.RData") # SKMEL5
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
# ===========================================================================================================
#               Cost function to find the residuals between the model and the data
# ===========================================================================================================
pars = data.frame(k_ab    = rep(0.025,1000),
                  k_ba    = rep(0.00000025,1000),
                  k_cb    = rep(0.025,1000),
                  k_bc    = rep(0.00000025,1000))

A0 = runif(2000, 0, 1)
dn = data.frame()
for(i in 1:length(A0)){
  sc = runif(1, 0, A0[i])
  dc = data.frame(A0 = A0[i], B0 = sc)
  dn = rbind(dn, dc)
}
dn = subset(dn, dn$A0+dn$B0<=1)
dn = sample_n(dn, 1000)
pars$A0 = dn$A0; pars$B0 = dn$B0
times <- times[order(times)]
nd = pars
# Objective function to calculate the cost--using modCost of R-package-FME
# ===========================================================================================================
#                         Compare model output with experimental data
# ===========================================================================================================
dn <- data.frame()
for(i in 1:1000){
  pars <- nd[i,]
  out <- as.data.frame(solvePop(pars, times))
  #out <- out[,c("time", "nl2")]
  out$id <- i
  dn <- rbind(dn, out)
}
# 
dn$cell <- "simulated"
data = dn

# ===========================================================================================================
# Early heterogeneity simulation
# ===========================================================================================================
df = subset(data, data$time<=80)
cl1 <- heat.colors(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~time, data=df, type="n", ylim=c(-2,2), xlim=c(0,80), ylab="", xlab="")
for(w in 1:length(unique(df$id))){
  temp <- subset(df, df$id==w)
  lines(temp$time, temp$nl2, col=cl1[w])
}
#Calculate rates of proliferation
dn = data.frame()
for(i in unique(df$id)){
  temp = subset(df, df$id==i & df$time<=80)
  dc = data.frame(id = unique(temp$id), rates = coef(lm(nl2~time, data=temp))[2])
  dn = rbind(dn, dc)
}
rates1 = dn
rates1$phase = "early"
# # ===========================================================================================================
# # late idling simulation
# # ===========================================================================================================
df = subset(data, data$time>=0)
cl1 <- rainbow(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~time, data=df, type="n", ylim=c(-2,2), xlim=c(0,350), ylab="", xlab="")
for(w in 1:length(unique(df$id))){
  temp <- subset(df, df$id==w)
  lines(temp$time, temp$nl2, col=cl1[w])
}
# Calculate rates of proliferation
dn = data.frame()
for(i in unique(df$id)){
  temp = subset(df, df$id==i & df$time>=10)
  dc = data.frame(id = unique(temp$id), rates = coef(lm(nl2~time, data=temp))[2])
  dn = rbind(dn, dc)
}
rates2 = dn
rates2$phase = "idling"
# # ===========================================================================================================
drates = rbind(rates1, rates2)
boxplot(rates~phase, data=drates, notch=T, ylim=c(-0.05, 0.05))
# # # ===========================================================================================================
# # 
# # # ===========================================================================================================
# # head(nd)
# # boxplot(nd$A0, nd$B0)
# at = subset(df, df$id==900)
# plot(nl2~time, data=at, ylim=c(min(df$nl2),max(df$nl2)), 
#      xlim=c(0,350), ylab="", xlab="", type="l")
# 
# # # ===========================================================================================================
# # # ===========================================================================================================
