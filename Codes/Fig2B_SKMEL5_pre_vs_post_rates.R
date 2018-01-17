#=====================================================================================================
# Supp. Figure 3A: Single-cell-derived SKMEL5 sublines initial response
# concentration used: 8uM PLX4720
# Scatterplot of untreated proliferation rates of sublines with their drug-treated growth rates
#=====================================================================================================
#=====================================================================================================
library(gplots)
require(ggplot2)

setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Fig2" # Which figure is this in paper?
output = paste0(output, "/", fig)
#=====================================================================================================
#=====================================================================================================

data <- read.csv("SKMEL5_cFP_Assay.csv", header=T, sep=",")
# Select 8uM PLX4720 and short-term response i.e. less that 120 hrs
data <- subset(data, data$conc==8)
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

d_rates1 <- data.frame(cID=dmso_rates$cID, dmso_rates = dmso_rates$rates, plx_rates = plx_rates$rates)
d_rates1 <- d_rates1[!is.na(d_rates1[,3]),]

#=====================================================================================================
# Read experimental data
#=====================================================================================================

d <- read.csv("SKMEL5_Subclones_cFP_assay_PLX4720.csv", header=T, sep=",")
# Select 8uM PLX4720 and short-term response i.e. less that 120 hrs
d <- subset(d, d$conc==8)
data <- subset(d, d$Cell=="SC01")
lmd <- function(l2){
  timePt <- unique(s$Time)
  coef(lm(l2~timePt))[2]
}
data_dmso <- subset(data, data$Time<=0) # Before the drug is added,
data_plx <- subset(data, data$Time>=0)  # After the drug is added, 
s <- subset(data_dmso)
dn <- data.frame()
for(j in unique(s$cID)){
    temp1 <- subset(s, s$cID==j)
    dc <- data.frame(cID = unique(temp1$cID) ,rates=coef(lm(l2~Time, data=temp1))[2])
    dn <- rbind(dn, dc)
  }
dmso_rates <- dn

s <- subset(data_plx, data_plx$Time>24 & data_plx$Time<120)
dn <- data.frame()

for(j in unique(s$cID)){
    temp1 <- subset(s, s$cID==j)
    dc <- data.frame(cID = unique(temp1$cID) ,rates=coef(lm(l2~Time, data=temp1))[2])
    dn <- rbind(dn, dc)
}
plx_rates <- dn

d_rates2 <- data.frame(cID=dmso_rates$cID, dmso_rates = dmso_rates$rates, plx_rates = plx_rates$rates)
d_rates2 <- d_rates2[!is.na(d_rates2[,3]),]

#=====================================================================================================
#
#=====================================================================================================
d_rates <- rbind(d_rates1, d_rates2)
#=====================================================================================================
# Plot scatter plot between pre and post growth rates
#=====================================================================================================
# scatterplot
pdf(paste0(output,"/","SKMEL5 pre vs post DIP rates (8uM PLX).pdf"), width=5, height=5)
par(ps = 12, cex = 1, cex.axis = 1)
d_pos <- subset(d_rates, d_rates$plx_rates>=0)
d_neg <- subset(d_rates, d_rates$plx_rates<0)

plot(dmso_rates~plx_rates, data=d_pos, xlab="", ylab="", 
     xlim=c(-0.1, 0.055), ylim=c(-0.055, 0.055), pch=1, col="black")
points(dmso_rates~plx_rates, data=d_neg,  pch=1, col="black")
abline(h=0, v=0, lty=2)
dev.off()

#=====================================================================================================

#=====================================================================================================

write.csv(d_rates, file="SKMEL5_cFP_assay_pre_vs_post_rates.csv")
