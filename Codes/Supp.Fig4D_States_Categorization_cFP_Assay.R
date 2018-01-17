# ===========================================================================================================
# Figure Supp.Fig 4D: State Categorization cFP Assay
# ===========================================================================================================
library(gplots)
require(ggplot2)
require(FME)
require(deSolve)
library(devtools)
library(ggbiplot)

setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Supp.Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)

# ===========================================================================================================
# Read experimental data
# ===========================================================================================================

d1 <- read.csv("SKMEL5 Subclone1 DIP.csv", header=T, sep=",")
d2 <- read.csv("SKMEL5 Subclone7 DIP.csv", header=T, sep=",")
d3 <- read.csv("SKMEL5 Subclone10 DIP.csv", header=T, sep=",")
colnames(d1) = c("X", "cID", "rates", "R2"); colnames(d2) = c("X", "cID", "rates", "R2");
colnames(d3) = c("X", "cID", "rates", "R2");
d4 <- read.csv("SKMEL19_cFP_rates.csv", header=T, sep=",")
d5 <- read.csv("WM88_cFP_rates.csv", header=T, sep=",")
d6 <- read.csv("SKMEL5_cFP_rates.csv", header=T, sep=",")

# ===========================================================================================================
# State categorization discretization: +- 1 doublings/2weeks
# ===========================================================================================================
# m +-10%*m
m <- 1/360
lower = m - (0.10*m)
upper = m + (0.10*m)
# ===========================================================================================================
# Distribution of cFP DIP rates for different cell lines
# ===========================================================================================================
par(ps = 12, cex = 1, cex.axis = 1)
# SKMEL19
d1A <- subset(d4, d4$rates< -m); d1LA = subset(d4, d4$rates< -lower); d1UA = subset(d4, d4$rates< -upper)
d1B <- subset(d4, d4$rates>= -m & d4$rates<=m); 
d1LB = subset(d4, d4$rates>= -lower & d4$rates<= lower); 
d1UB = subset(d4, d4$rates>= -upper & d4$rates<= upper);
d1C <- subset(d4, d4$rates> m); d1LC = subset(d4, d4$rates> lower); d1UC = subset(d4, d4$rates> upper)


a1 = round(length(d1A$rates)/length(d4$rates), 4)
a2 = round(length(d1LA$rates)/length(d4$rates), 4)
a3 = round(length(d1UA$rates)/length(d4$rates), 4)

b1 = round(length(d1B$rates)/length(d4$rates), 4)
b2 = round(length(d1LB$rates)/length(d4$rates), 4)
b3 = round(length(d1UB$rates)/length(d4$rates), 4)

c1 = round(length(d1C$rates)/length(d4$rates), 4)
c2 = round(length(d1LC$rates)/length(d4$rates), 4)
c3 = round(length(d1UC$rates)/length(d4$rates), 4)

# Plot histogram and divide groups
pdf(paste0(output, "/SKMEL19_cFP_States.pdf"), width=4, height=4)
hist(d4$rates, freq=F, breaks=seq(-0.10, 0.06, by=0.0025), main="SKMEL19", xlab="", ylab="", 
     xlim=c(-0.1, 0.06), ylim=c(0,100))
#abline(v=0, col="black", lty=2, lwd=2)
abline(v=-m, col="red", lwd=2, lty=2)
abline(v=m, col="red", lwd=2, lty=2)
text(-0.05, 80, paste0("R0 = ", a1))
text(-0.05, 70, paste0("S0 = ", b1))
text(-0.05, 60, paste0("E0 = ", c1))
dev.off()

#vioplot(d5$rates, horizontal = T, ylim=c(-0.1, 0.06))

cell1 <- data.frame(cell = rep(unique(d4$Cell), 3), A0=c(a1, a2, a3), B0=c(b1, b2, b3), C0=c(c1, c2, c3))
# ===========================================================================================================
# WM88
d1A <- subset(d5, d5$rates< -m); d1LA = subset(d5, d5$rates< -lower); d1UA = subset(d5, d5$rates< -upper)
d1B <- subset(d5, d5$rates>= -m & d5$rates<=m); 
d1LB = subset(d5, d5$rates>= -lower & d5$rates<= lower); 
d1UB = subset(d5, d5$rates>= -upper & d5$rates<= upper);
d1C <- subset(d5, d5$rates> m); d1LC = subset(d5, d5$rates> lower); d1UC = subset(d5, d5$rates> upper)


a1 = round(length(d1A$rates)/length(d5$rates), 4)
a2 = round(length(d1LA$rates)/length(d5$rates), 4)
a3 = round(length(d1UA$rates)/length(d5$rates), 4)

b1 = round(length(d1B$rates)/length(d5$rates), 4)
b2 = round(length(d1LB$rates)/length(d5$rates), 4)
b3 = round(length(d1UB$rates)/length(d5$rates), 4)

c1 = round(length(d1C$rates)/length(d5$rates), 4)
c2 = round(length(d1LC$rates)/length(d5$rates), 4)
c3 = round(length(d1UC$rates)/length(d5$rates), 4)

# Plot histogram and divide groups
pdf(paste0(output, "/WM88_cFP_States.pdf"), width=4, height=4)
hist(d5$rates, freq=F, breaks=seq(-0.10, 0.06, by=0.0025), main="WM88", xlab="", ylab="", 
     xlim=c(-0.1, 0.06), ylim=c(0,100))
#abline(v=0, col="black", lty=2, lwd=2)
abline(v=-m, col="red", lwd=2, lty=2)
abline(v=m, col="red", lwd=2, lty=2)
text(-0.05, 80, paste0("R0 = ", a1))
text(-0.05, 70, paste0("S0 = ", b1))
text(-0.05, 60, paste0("E0 = ", c1))
dev.off()
cell2 <- data.frame(cell = rep(unique(d5$Cell), 3), A0=c(a1, a2, a3), B0=c(b1, b2, b3), C0=c(c1, c2, c3))
# ===========================================================================================================
# SC01
d1A <- subset(d1, d1$rates< -m); d1LA = subset(d1, d1$rates< -lower); d1UA = subset(d1, d1$rates< -upper)
d1B <- subset(d1, d1$rates>= -m & d1$rates<=m); 
d1LB = subset(d1, d1$rates>= -lower & d1$rates<= lower); 
d1UB = subset(d1, d1$rates>= -upper & d1$rates<= upper);
d1C <- subset(d1, d1$rates> m); d1LC = subset(d1, d1$rates> lower); d1UC = subset(d1, d1$rates> upper)


a1 = round(length(d1A$rates)/length(d1$rates), 4)
a2 = round(length(d1LA$rates)/length(d1$rates), 4)
a3 = round(length(d1UA$rates)/length(d1$rates), 4)

b1 = round(length(d1B$rates)/length(d1$rates), 4)
b2 = round(length(d1LB$rates)/length(d1$rates), 4)
b3 = round(length(d1UB$rates)/length(d1$rates), 4)

c1 = round(length(d1C$rates)/length(d1$rates), 4)
c2 = round(length(d1LC$rates)/length(d1$rates), 4)
c3 = round(length(d1UC$rates)/length(d1$rates), 4)
# Plot histogram and divide groups
pdf(paste0(output, "/SC01_cFP_States.pdf"), width=4, height=4)
hist(d1$rates, freq=F, breaks=seq(-0.10, 0.06, by=0.0025), main="SC01", xlab="", ylab="", 
     xlim=c(-0.1, 0.06), ylim=c(0,100))
#abline(v=0, col="black", lty=2, lwd=2)
abline(v=-m, col="red", lwd=2, lty=2)
abline(v=m, col="red", lwd=2, lty=2)
text(-0.05, 80, paste0("R0 = ", a1))
text(-0.05, 70, paste0("S0 = ", b1))
text(-0.05, 60, paste0("E0 = ", c1))
dev.off()

cell3 <- data.frame(cell = rep("SC01", 3), A0=c(a1, a2, a3), B0=c(b1, b2, b3), C0=c(c1, c2, c3))
# ===========================================================================================================
# SC07
d1A <- subset(d2, d2$rates< -m); d1LA = subset(d2, d2$rates< -lower); d1UA = subset(d2, d2$rates< -upper)
d1B <- subset(d2, d2$rates>= -m & d2$rates<=m); 
d1LB = subset(d2, d2$rates>= -lower & d2$rates<= lower); 
d1UB = subset(d2, d2$rates>= -upper & d2$rates<= upper);
d1C <- subset(d2, d2$rates> m); d1LC = subset(d2, d2$rates> lower); d1UC = subset(d2, d2$rates> upper)


a1 = round(length(d1A$rates)/length(d2$rates), 4)
a2 = round(length(d1LA$rates)/length(d2$rates), 4)
a3 = round(length(d1UA$rates)/length(d2$rates), 4)

b1 = round(length(d1B$rates)/length(d2$rates), 4)
b2 = round(length(d1LB$rates)/length(d2$rates), 4)
b3 = round(length(d1UB$rates)/length(d2$rates), 4)

c1 = round(length(d1C$rates)/length(d2$rates), 4)
c2 = round(length(d1LC$rates)/length(d2$rates), 4)
c3 = round(length(d1UC$rates)/length(d2$rates), 4)
# Plot histogram and divide groups
pdf(paste0(output, "/SC07_cFP_States.pdf"), width=4, height=4)
hist(d2$rates, freq=F, breaks=seq(-0.10, 0.06, by=0.0025), main="SC07", xlab="", ylab="", 
     xlim=c(-0.1, 0.06), ylim=c(0,100))
#abline(v=0, col="black", lty=2, lwd=2)
abline(v=-m, col="red", lwd=2, lty=2)
abline(v=m, col="red", lwd=2, lty=2)
text(-0.05, 80, paste0("R0 = ", a1))
text(-0.05, 70, paste0("S0 = ", b1))
text(-0.05, 60, paste0("E0 = ", c1))
dev.off()

#cell4 <- data.frame(cell = "SC07", R0=a, S0=b, E0=c)
cell4 <- data.frame(cell = rep("SC07", 3), A0=c(a1, a2, a3), B0=c(b1, b2, b3), C0=c(c1, c2, c3))
# ===========================================================================================================
# SC10
d1A <- subset(d3, d3$rates< -m); d1LA = subset(d3, d3$rates< -lower); d1UA = subset(d3, d3$rates< -upper)
d1B <- subset(d3, d3$rates>= -m & d3$rates<=m); 
d1LB = subset(d3, d3$rates>= -lower & d3$rates<= lower); 
d1UB = subset(d3, d3$rates>= -upper & d3$rates<= upper);
d1C <- subset(d3, d3$rates> m); d1LC = subset(d3, d3$rates> lower); d1UC = subset(d3, d3$rates> upper)


a1 = round(length(d1A$rates)/length(d3$rates), 4)
a2 = round(length(d1LA$rates)/length(d3$rates), 4)
a3 = round(length(d1UA$rates)/length(d3$rates), 4)

b1 = round(length(d1B$rates)/length(d3$rates), 4)
b2 = round(length(d1LB$rates)/length(d3$rates), 4)
b3 = round(length(d1UB$rates)/length(d3$rates), 4)

c1 = round(length(d1C$rates)/length(d3$rates), 4)
c2 = round(length(d1LC$rates)/length(d3$rates), 4)
c3 = round(length(d1UC$rates)/length(d3$rates), 4)
# Plot histogram and divide groups
pdf(paste0(output, "/SC10_cFP_States.pdf"), width=4, height=4)
hist(d3$rates, freq=F, breaks=seq(-0.10, 0.06, by=0.0025), main="SC10", xlab="", ylab="", 
     xlim=c(-0.1, 0.06), ylim=c(0,100))
#abline(v=0, col="black", lty=2, lwd=2)
abline(v=-m, col="red", lwd=2, lty=2)
abline(v=m, col="red", lwd=2, lty=2)
text(-0.05, 80, paste0("R0 = ", a1))
text(-0.05, 70, paste0("S0 = ", b1))
text(-0.05, 60, paste0("E0 = ", c1))
dev.off()

#cell5 <- data.frame(cell = "SC10", R0=a, S0=b, E0=c)
cell5 <- data.frame(cell = rep("SC10", 3), A0=c(a1, a2, a3), B0=c(b1, b2, b3), C0=c(c1, c2, c3))
# ===========================================================================================================
# SKMEL5
d6 <- data.frame(rates = c(d1$rates, d2$rates, d3$rates, d6$rates))

d1A <- subset(d6, d6$rates< -m); d1LA = subset(d6, d6$rates< -lower); d1UA = subset(d6, d6$rates< -upper)
d1B <- subset(d6, d6$rates>= -m & d6$rates<=m); 
d1LB = subset(d6, d6$rates>= -lower & d6$rates<= lower); 
d1UB = subset(d6, d6$rates>= -upper & d6$rates<= upper);
d1C <- subset(d6, d6$rates> m); d1LC = subset(d6, d6$rates> lower); d1UC = subset(d6, d6$rates> upper)

a1 = round(length(d1A$rates)/length(d6$rates), 4)
a2 = round(length(d1LA$rates)/length(d6$rates), 4)
a3 = round(length(d1UA$rates)/length(d6$rates), 4)

b1 = round(length(d1B$rates)/length(d6$rates), 4)
b2 = round(length(d1LB$rates)/length(d6$rates), 4)
b3 = round(length(d1UB$rates)/length(d6$rates), 4)

c1 = round(length(d1C$rates)/length(d6$rates), 4)
c2 = round(length(d1LC$rates)/length(d6$rates), 4)
c3 = round(length(d1UC$rates)/length(d6$rates), 4)
# Plot histogram and divide groups
pdf(paste0(output,"/SKMEL5_cFP_States.pdf"), width=4, height=4)
hist(d6$rates, freq=F, breaks=seq(-0.10, 0.06, by=0.0025), main="SKMEL5", xlab="", ylab="", 
     xlim=c(-0.1, 0.06), ylim=c(0,100))
#abline(v=0, col="black", lty=2, lwd=2)
abline(v=-m, col="red", lwd=2, lty=2)
abline(v=m, col="red", lwd=2, lty=2)
text(-0.05, 80, paste0("R0 = ", a1))
text(-0.05, 70, paste0("S0 = ", b1))
text(-0.05, 60, paste0("E0 = ", c1))
dev.off()

#cell6 <- data.frame(cell = "SKMEL5", R0=a, S0=b, E0=c)
cell6 <- data.frame(cell = rep("SKMEL5", 3), A0=c(a1, a2, a3), B0=c(b1, b2, b3), C0=c(c1, c2, c3))

data <- rbind(cell1, cell2, cell3, cell4, cell5, cell6)
write.csv(data, file="State_proportions_cFP_assay.csv")
# ===========================================================================================================

pdf(paste0(output,"/State_categories.pdf"), width=4, height=4)
hist(d6$rates, freq=F, breaks=seq(-0.10, 0.06, by=0.0025), main="", xlab="", ylab="", 
     xlim=c(-0.1, 0.06), ylim=c(0,100), yaxt="n")
abline(v=-m, col="blue", lwd=2, lty=2)
abline(v=m, col="blue", lwd=2, lty=2)
text(-(m+0.009), 60, expression(-m), cex=1.2)
text((m+0.009), 60, expression(m), cex=1.2)
text(0, 80, "S", cex=1.2)
text(-0.03, 80, "R", cex=1.2)
text(0.02, 80, "E", cex=1.2)
dev.off()
# ===========================================================================================================
# Functional relevance of discrete categorization of three states
# ===========================================================================================================
value = 1/360; rate1 <- seq(-value, value, by=0.0001); rate2 <- seq(-0.1, -value, by=0.0001)
rate3 <- seq(value, 0.05, by=0.0001); times <- seq(1,300); times <- c(0, times)
y_in = 10000
# ODE solver
solvePop <- function(pars, times, Tot = 10000){
  derivs <- function(t, state, pars){
    with(as.list(c(state, pars)),{
      dA <- pars * A
      return(list(c(dA), nl2 = (log2(A)-log2(Tot))))
    })
  }
  state <- with(as.list(pars), c(A = Tot))
  output = ode(y=state, times=times, func=derivs, parm=pars, method="lsoda", hmax=1, jacfunc=mod)
  print(diagnostics(output))
  return(output)
}

# ===========================================================================================================
dn <- data.frame()
for(i in 1:length(rate1)){
  pars <- rate1[i]
  out <- as.data.frame(solvePop(pars, times))
  #out <- out[,c("time", "nl2")]
  out$id <- i
  dn <- rbind(dn, out)
}
td1 <- dn
#
dn <- data.frame()
for(i in 1:length(rate2)){
  pars <- rate2[i]
  out <- as.data.frame(solvePop(pars, times))
  #out <- out[,c("time", "nl2")]
  out$id <- i
  dn <- rbind(dn, out)
}
td2 <- dn
#
dn <- data.frame()
for(i in 1:length(rate3)){
  pars <- rate3[i]
  out <- as.data.frame(solvePop(pars, times))
  #out <- out[,c("time", "nl2")]
  out$id <- i
  dn <- rbind(dn, out)
}
td3 <- dn

cl1 <- grey.colors(length(unique(td1$id)))
cl2 <- cm.colors(length(unique(td2$id)))
cl3 <- heat.colors(length(unique(td3$id)))
pdf(paste0(output, "/Functional_categories_relevance.pdf"), width=3., height=3.2)
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~time, data=td1, type="n", ylim=c(-5,5), xlim=c(0,300), ylab="", xlab="")
for(w in 1:length(unique(td1$id))){
  temp <- subset(td1, td1$id==w)
  lines(temp$time, temp$nl2, col=cl1[w])
}
for(w in 1:length(unique(td2$id))){
  temp <- subset(td2, td2$id==w)
  lines(temp$time, temp$nl2, col=cl2[w])
}
for(w in 1:length(unique(td3$id))){
  temp <- subset(td3, td3$id==w)
  lines(temp$time, temp$nl2, col=cl3[w])
}
abline(h=c(-1,1), col="blue", lty=2, lwd=2)
dev.off()
# ===========================================================================================================




