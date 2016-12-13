#=====================================================================================================
# Figure 2: Short-term population response clonal heterogeneity
#=====================================================================================================
# SKMEL5 short-term response in 8uM PLX4720 
# Clonal fractional proliferation assay--single-cell-derived responses
#=====================================================================================================
#=====================================================================================================
library(gplots)
require(ggplot2)
setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Supp.Fig2" # Which figure is this in paper?
output = paste0(output, "/", fig)

#=====================================================================================================
# Read the experimental data for SKMEL5 clonal heterogeneity
#=====================================================================================================

data <- read.csv("SKMEL5_cFP_Assay.csv", header=T, sep=",")
# Select 8uM PLX4720 and short-term response i.e. less that 120 hrs
data <- subset(data, data$conc==8)

#=====================================================================================================
# CALCULATION OF DIP RATES TO GET A DISTRIBUTION
# Define the linear fit model to calculate the DIP rates for each colonies
#=====================================================================================================

lmd <- function(l2){
  timePt <- unique(s$Time)
  coef(lm(l2~timePt))[2]
}

data_dmso <- subset(data, data$Time<=0) # Before the drug is added,
data_plx <- subset(data, data$Time>=0)  # After the drug is added, 

s <- subset(data_dmso)
dn <- data.frame()
for(i in unique(s$Date)){
  temp <- subset(s, s$Date==i)
  for(j in unique(temp$cID)){
    temp1 <- subset(temp, temp$cID==j)
    dc <- data.frame(cID = unique(temp1$cID) ,rates=coef(lm(l2~Time, data=temp1))[2])
    dn <- rbind(dn, dc)
  }
}
dmso_rates <- dn

s <- subset(data_plx, data_plx$Time>24 & data_plx$Time<120)
dn <- data.frame()
for(i in unique(s$Date)){
  temp <- subset(s, s$Date==i)
  for(j in unique(temp$cID)){
    temp1 <- subset(temp, temp$cID==j)
    dc <- data.frame(cID = unique(temp1$cID) ,rates=coef(lm(l2~Time, data=temp1))[2])
    dn <- rbind(dn, dc)
  }
}
plx_rates <- dn

d_rates <- data.frame(cID=dmso_rates$cID, dmso_rates = dmso_rates$rates, plx_rates = plx_rates$rates)
d_rates <- d_rates[!is.na(d_rates[,3]),]

#=====================================================================================================
# Plot scatter plot between pre and post growth rates
#=====================================================================================================
# scatterplot
pdf(paste0(output,"/","SKMEL5 pre vs post DIP rates (8uM PLX).pdf"), width=5, height=5)
par(ps = 12, cex = 1, cex.axis = 1)
d_pos <- subset(d_rates, d_rates$plx_rates>=0)
d_neg <- subset(d_rates, d_rates$plx_rates<0)

plot(dmso_rates~plx_rates, data=d_pos, xlab="", ylab="", 
     xlim=c(-0.055, 0.055), ylim=c(-0.055, 0.055), pch=1, col="black")
points(dmso_rates~plx_rates, data=d_neg,  pch=1, col="black")
abline(h=0, v=0, lty=2)
dev.off()

m <- lm(dmso_rates~plx_rates, data=d_rates)
summary(m)

#=====================================================================================================
# Plot individual clonal responses in the short-term in different concentrations of PLX4720
#=====================================================================================================
data <- read.csv("SKMEL5_cFP_Assay.csv", header=T, sep=",")
data <- subset(data, data$Time >=0 & data$Time <=96)
data1 <- subset(data, data$conc==2)
data2 <- subset(data, data$conc==16)

pdf(paste0(output,"/","SKMEL5 + 2uM PLX4720 cFP Assay.pdf"), width=3, height=4)
cl1 <- rainbow(length(unique(data1$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~Time, data=data1, type="n", ylim=c(-3,3), xlim=c(0,105), ylab="", xlab="")
for(w in 1:length(unique(data1$cID))){
  temp <- subset(data1, data1$cID==w)
  lines(temp$Time, temp$nl2, col=cl1[w])
}
dev.off()

pdf(paste0(output,"/","SKMEL5 + 16uM PLX4720 cFP Assay.pdf"), width=3, height=4)
cl1 <- rainbow(length(unique(data2$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~Time, data=data2, type="n", ylim=c(-3,3), xlim=c(0,105), ylab="", xlab="")
for(w in 1:length(unique(data2$cID))){
  temp <- subset(data2, data2$cID==w)
  lines(temp$Time, temp$nl2, col=cl1[w])
}
dev.off()

#=====================================================================================================
# Plot the heterogeneity for other cell lines: SKMEL19 and WM88
#=====================================================================================================

data <- read.csv("SKMEL19+WM88 cFP Assay.csv", header=T, sep=",")
# Select 8uM PLX4720 and short-term response i.e. less that 120 hrs
data <- subset(data, data$Time >=0 & data$Time <=96)
data1 <- subset(data, data$CellLine=="WM88")
data2 <- subset(data, data$CellLine=="SKMEL19")

pdf(paste0(output,"/","WM88 + 8uM PLX4720 cFP Assay.pdf"), width=3, height=4)
cl1 <- rainbow(length(unique(data1$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~Time, data=data1, type="n", ylim=c(-3,3), xlim=c(0,105), ylab="", xlab="")
for(w in 1:length(unique(data1$cID))){
  temp <- subset(data1, data1$cID==w)
  lines(temp$Time, temp$nl2, col=cl1[w])
}
dev.off()

pdf(paste0(output,"/","SKMEL19 + 8uM PLX4720 cFP Assay.pdf"), width=3, height=4)
cl1 <- rainbow(length(unique(data2$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~Time, data=data2, type="n", ylim=c(-3,3), xlim=c(0,105), ylab="", xlab="")
for(w in 1:length(unique(data2$cID))){
  temp <- subset(data2, data2$cID==w)
  lines(temp$Time, temp$nl2, col=cl1[w])
}
dev.off()

#=====================================================================================================
#=====================================================================================================
