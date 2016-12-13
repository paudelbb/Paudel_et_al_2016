#=====================================================================================================
# Supp. Figure 1:  Population level response of SKMEL5 in DMSO vs 8uM PLX4720
#=====================================================================================================
#         Population level response of SKMEL5 to PLX4720
#         Drug-naive vs Post-idle cells treated with either DMSO or 8uM PLX4720
#=====================================================================================================

library(gplots)
require(ggplot2)
require(gdata)
setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures")
fig = "Supp.Fig1" # Which figure is this in paper?
output = paste0(output, "/", fig)

#=====================================================================================================
# Read the experimental data
#=====================================================================================================

data <- read.csv("SKMEL5_Parental_drug_removal.csv", header=T, sep=",")

#=====================================================================================================
# SKMEL5 drug-naive vs post-idle treated with DMSO
#=====================================================================================================

cnc <- 0
s1 <- subset(data, data$conc==cnc)
s1$conc <- as.character(s1$conc)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2)) +
  theme_bw()+ geom_smooth(aes(linetype=cond, fill=cond), span=.5, data=s1, method = "loess", 
  size=0.9, alpha=0.4, se=TRUE, col="black")+ 
  scale_colour_manual(values=c("blue"), name="a") + ylim(-.5, 4.)+ xlim(0,100)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") + 
  ggsave(paste0(unique(data$CellLine), " + ", "drug_removal_", unique(s1$conc),"μΜ.pdf"), 
         path=output, width=3, height=3)

#=====================================================================================================
# SKMEL5 drug-naive vs post-idle treated with 8uM PLX4720
#=====================================================================================================

cnc <- 8
s1 <- subset(data, data$conc==cnc)
s1$conc <- as.character(s1$conc)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2)) +
  theme_bw()+ geom_smooth(aes(linetype=cond, fill=cond), span=.5, data=s1, method = "loess", 
  size=0.9, alpha=0.4, se=TRUE, col="black")+ 
  scale_colour_manual(values=c("black","blue"), name="") + ylim(-0.5, 1.)+ xlim(0,100)+
  theme(legend.position=c(0.4, 0.8), legend.position="horizontal")+ theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") + 
  ggsave(paste0(unique(data$CellLine), " + ", "drug_removal_", unique(s1$conc),"μΜ.pdf"), 
         path=output, width=3, height=3)

#=====================================================================================================
