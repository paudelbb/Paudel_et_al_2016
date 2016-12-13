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

data <- read.csv("SKMEL5 Population Level DMSO vs 8uM PLX4720.csv", header=T, sep=",")

#=====================================================================================================
# SKMEL5 drug-naive vs post-idle treated with DMSO
#=====================================================================================================

cell <- "SKMEL5"
cnc <- c(0,8)
s1 <- subset(data, data$CellLine==cell & data$conc %in%  cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=drug))+ 
  theme_bw()+geom_smooth(span=.3, aes(fill=drug), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values=rainbow(2)) + ylim(-0.5,6.5)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=10)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc[2], "uM_vs_", "DMSO.pdf"), path=output, width=2, height=3.0)

