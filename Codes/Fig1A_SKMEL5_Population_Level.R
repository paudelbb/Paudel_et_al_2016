#=====================================================================================================
# Figure 1: Population Level Response for different cell lines
#=====================================================================================================
#         Population level response of BRAF-mutated melanoma cell lines to PLX4720
#         CellLines used: SKMEL5, SKMEL19, SKMEL28, WM88, WM164, WM793, A375
#                 Concentrations for idling phenotype in cell lines:
#         SKMEL5(8uM), SKMEL19(8uM), SKMEL28(16uM), WM88(8uM), WM164(32uM), WM793(16uM), A375(8uM)
#=====================================================================================================
#=====================================================================================================

library(gplots)
require(ggplot2)
setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Fig1" # Which figure is this in paper?
output = paste0(output, "/", fig)
#=====================================================================================================
# Read the experimental data
#=====================================================================================================

data <- read.csv("Population_Response_CellLines.csv", header=T, sep=",")
cell <- "SKMEL5"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

#=====================================================================================================
# Get the plot for the population level response
#=====================================================================================================

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-0.5,1.1)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)

#=====================================================================================================

#=====================================================================================================




