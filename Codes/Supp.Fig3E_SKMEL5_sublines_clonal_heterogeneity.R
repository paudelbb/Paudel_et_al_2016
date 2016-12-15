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
fig = "Supp.Fig3" # Which figure is this in paper?
output = paste0(output, "/", fig)

#=====================================================================================================
# Read experimental data
#=====================================================================================================

d <- read.csv("SKMEL5_Subclones_cFP_assay_PLX4720.csv", header=T, sep=",")
# Select 8uM PLX4720 and short-term response i.e. less that 120 hrs
d <- subset(d, d$conc==8)
d <- subset(d, d$Time >=0 & d$Time <=100)

#=====================================================================================================
# Plot clonal response for sublines
#=====================================================================================================
# SC01
data <- subset(d, d$Cell=="SC01")
length(unique(data$cID))
pdf(paste0(output,"/", unique(data$Cell), "_8uM PLX4720 cFP Assay.pdf"), width=3, height=4)
cl1 <- rainbow(length(unique(data$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~Time, data=data, type="n", ylim=c(-5.5,3), xlim=c(0,105), ylab="", xlab="")
for(w in 1:length(unique(data$cID))){
  temp <- subset(data, data$cID==w)
  lines(temp$Time, temp$nl2, col=cl1[w])
}
dev.off()

#SC07
data <- subset(d, d$Cell=="SC07")
length(unique(data$cID))
pdf(paste0(output,"/", unique(data$Cell), "_8uM PLX4720 cFP Assay.pdf"), width=3, height=4)
cl1 <- rainbow(length(unique(data$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~Time, data=data, type="n", ylim=c(-5.5,3), xlim=c(0,105), ylab="", xlab="")
for(w in 1:length(unique(data$cID))){
  temp <- subset(data, data$cID==w)
  lines(temp$Time, temp$nl2, col=cl1[w])
}
dev.off()

#SC10
data <- subset(d, d$Cell=="SC10")
length(unique(data$cID))
pdf(paste0(output,"/", unique(data$Cell), "_8uM PLX4720 cFP Assay.pdf"), width=3, height=4)
cl1 <- rainbow(length(unique(data$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~Time, data=data, type="n", ylim=c(-5.5,3), xlim=c(0,105), ylab="", xlab="")
for(w in 1:length(unique(data$cID))){
  temp <- subset(data, data$cID==w)
  lines(temp$Time, temp$nl2, col=cl1[w])
}
dev.off()

#=====================================================================================================