#=====================================================================================================
# Supp. Figure 3BC: Single-cell-derived SKMEL5 sublines 03,04,08 response and  
# fraction of control in Subclones, 01, 07 and 10 in idling state
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

data <- read.csv("SKMEL5_Subclones030408_Population_Level.csv", header=T, sep=",")
data <- data[(data$conc %in% c(8)),]
data <- subset(data, data$Time <= 250)

#=====================================================================================================
# Plot the responses of single-cell-derived subclones
#=====================================================================================================

ggplot(data = data, aes(x=data$Time, y=data$nl2, col=CellLine)) +
  theme_bw()+ 
  geom_smooth(span=.45, aes(fill=CellLine),data=data, method = "loess", size=.5, level=0.95, col="blue")+ 
  scale_colour_manual(values=c("black", "black", "black")) + ylim(-3, 3)+
  theme(legend.position="none") +
  ggtitle("")+ labs(x="", y="") +
  theme(axis.text=element_text(size=12))+theme(text = element_text(size=12)) +
  ggsave("SKMEL5 Subclones_030408 + 8uM PLX4720.pdf", path=output, width=3, height=3.5)

#=====================================================================================================
# Read experimental data for subclones idling state in SC01, SC07 and SC10
# That idling is not due to cell confluency
#=====================================================================================================

data <- read.csv("SKMEL5_Subclones_Fraction_of_control_idling_state.csv", header=T, sep=",")
pdf(paste0(output,"/","Fraction of control_idling_subclones.pdf"), width=3.5, height=3.5)
boxplot(fCont~CellLine, data=data, ylim=c(0,1))
dev.off()

#=====================================================================================================

#=====================================================================================================

data <- read.csv("SKMEL5_Subclones030408_Population_Level.csv", header=T, sep=",")
data <- data[(data$conc %in% c(0,8)),]
data <- subset(data, data$Time ==208.5317)
dn = data.frame()
for(i in unique(data$CellLine)){
  temp = subset(data, data$CellLine==i)
  for(j in unique(temp$rep)){
    t = subset(temp, temp$rep==j)
    t$rel = t$Cell.Nucleus/t$Cell.Nucleus[t$conc==0]
    dn = rbind(dn, t)
  }
}
dn = subset(dn, dn$conc==8)
pdf(paste0(output, "/SKMEL5_SC030408_idling_state_cell_fraction.pdf"), height=4, width=3)
par(ps = 8, cex = 1, cex.axis = 1)
boxplot(rel~CellLine, data=dn, ylim=c(0,1.), main="", xlab="")
dev.off()

