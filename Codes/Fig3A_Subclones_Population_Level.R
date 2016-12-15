#=====================================================================================================
# Figure 3: Population Level Response for SKMEL5 Single-cell-derived sublines
#=====================================================================================================
#         Population level response of Sublines to PLX4720
#         CellLines used: SC01, SC07 and SC10
#         Concentrations used: 8uM
#=====================================================================================================
#=====================================================================================================
library(gplots)
require(ggplot2)

setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Fig3" # Which figure is this in paper?
output = paste0(output, "/", fig)

#=====================================================================================================
# Read the experimental data
#=====================================================================================================

data <- read.csv("SKMEL5_Subclones_Response.csv", header=T, sep=",")

#=====================================================================================================
# Plot the response
#=====================================================================================================
ggplot(data = data, aes(x=data$Time, y=data$nl2, col=CellLine)) +
  theme_bw()+ 
  geom_smooth(span=.45, aes(fill=CellLine),data=data, method = "loess", size=.5, alpha=0.6)+ 
  scale_colour_manual(values=c("red", "green", "deepskyblue")) + ylim(-2, 2)+
  theme(legend.position="none") +
  ggtitle("")+ labs(x="", y="") +
  theme(axis.text=element_text(size=12))+theme(text = element_text(size=12)) +
  ggsave("SKMEL5 Subclones + 8uM PLX4720.pdf", path=output, width=3, height=3)

#=====================================================================================================