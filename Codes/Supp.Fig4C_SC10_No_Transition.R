# ===========================================================================================================
# Figure 4C: SKMEL5 Parental and Subclones No Transitions
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

# ===========================================================================================================
# Read experimental data
# ===========================================================================================================

data <- read.csv("SKMEL5_Subclones_Population_Level.csv", header=T, sep=",")
cell <- "SKMEL5 Subclone10"; cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
data <- s1

# ===========================================================================================================
#                       Ordinary differential equation models and initial parameters
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
#               Cost function to find the residuals between the model and the data
# ===========================================================================================================
pars <- c(k_ab    = 0,
          k_ba    = 0,
          k_cb    = 0,
          k_bc    = 0,
          A0 = 0, 
          B0 = 0.)
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
#                         Compare model output with experimental data
# ===========================================================================================================
out <- as.data.frame(solvePop(pars, times)); 
d1 = out[,c("time", "nl2")]
d1$CellLine = cell
d1 <- d1[,c("CellLine", "time", "nl2")]; 
colnames(d1) <- c("CellLine", "Time", "nl2");
s1 <- s1[,c("CellLine", "Time", "nl2")]
data <- s1
ggplot(data = data, aes(x=data$Time, y=data$nl2, col=CellLine)) +
  theme_bw() + geom_smooth(span=.4, aes(fill=CellLine), fill = "grey", alpha=0.4, data=data, 
  method = "loess", size=.5, level=0.95)+ 
  scale_colour_manual(values=c("black")) + ylim(-2.5, 2.5) + xlim(0,300)+
  theme(legend.position="none") + 
  geom_line(aes(x=d1$Time, y=d1$nl2), data=d1, col="red",size=.6)+ ggtitle("")+ labs(x="", y="") +
  theme(axis.text=element_text(size=12))+theme(text = element_text(size=12)) +
  ggsave(paste0("No_transitions_", cell, ".pdf"), path=output, width=3, height=3)

# ===========================================================================================================
