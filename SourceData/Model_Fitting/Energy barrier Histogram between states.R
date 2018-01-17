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
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData/Model_Fitting"))), "Figures") 
fig = "Supp.Fig7a" # Which figure is this in paper?
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
chain1 <- load("MCMC_Subclones_Mix_150000_Run.RData")
dim(MCMC1$pars)
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
rates <- data.frame(MCMC1$pars[,c(1,2,3,4)])

# ===========================================================================================================
# Calculate the activation energy based on the transition rate constants between the states
# ===========================================================================================================
# Transition from state A to state B
n = 75000
Ea = -log(rates$k_ab)
Eba = -log(rates$k_ba)
Ebc = -log(rates$k_bc)
Ec = -log(rates$k_cb)
# ===========================================================================================================
# Find the energy difference between states
# ===========================================================================================================
td1 = data.frame(A = Ea, B = 0, C = Ec)
pdf(paste0(output, "/Energy barrier to state B from other states.pdf"), width=4, height=4)
boxplot(td1, outline = F, notch = T, xlab ="to state B from:", ylim = c(0, 14))
dev.off()

td2 = data.frame(A = Eba, B = 0, C = Ebc)
pdf(paste0(output, "/Energy barrier from state B to other states.pdf"), width=4, height=4)
boxplot(td2, outline = F, notch = T, xlab ="from state B to:", ylim = c(0, 14))
dev.off()

td3 = data.frame(AB = Ea, BA = Eb1, BB = 0, BC = Eb2, CB = Ec)
boxplot(td3, outline = T, notch = T, xlab ="from state i to state j (ij)", ylim = c(0, 20))

dt1 = data.frame(Energy =  Ea, type = "AB")
dt2 = data.frame(Energy = Eba, type = "BA")
dt3 = data.frame(Energy =  0, type = "BB")
dt4 = data.frame(Energy = Ebc, type = "BC")
dt5 = data.frame(Energy =  Ec, type = "CB")
dtr = rbind(dt1, dt2, dt3, dt4, dt5)

library(ggplot2)
ggplot(dtr, aes(x = factor(type), y = Energy, fill = factor(type))) + ylim(0,20)+
  geom_violin() + ylab("Relative height") + xlab("from state i to j (ij)") + theme_bw()+
  theme(legend.position="none") + theme(axis.text=element_text(size=10)) +
  theme(text = element_text(size=10)) + ggtitle(paste0()) + 
  ggsave(paste0("Energy barrier violin plot.pdf"), path=output, width=3.5, height=3)

# ===========================================================================================================
# ===========================================================================================================

# size between A and B; 
a1 = Eba/Ea # energy forward divided by energy backward
a2 = Ebc/Ec # energy forward divided by energy backward

pdf(paste0(output, "/Ratio of Energy Barrier between states.pdf"), width = 3.5, height = 3.5)
hist(a1, freq = F, breaks = 30, xlim = c(0, 5), ylim = c(0,2.0), col = rgb(0,1,1,0.5), xlab = "Ratio", main="")
hist(a2, freq = F, breaks = 30, add = T, col = rgb(0,0,0,0.5))
legend("topleft", c("E(B<->A)", "E(B<->C)"), fill=c(rgb(0,1,1,0.5), rgb(0,0,0,0.5)), bty = "n")
dev.off()
# ===========================================================================================================
# ===========================================================================================================


