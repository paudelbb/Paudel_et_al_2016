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
fig = "Fig2" # Which figure is this in paper?
output = paste0(output, "/", fig)

#=====================================================================================================
# Read the experimental data
#=====================================================================================================

data <- read.csv("SKMEL5_cFP_Assay.csv", header=T, sep=",")
# Select 8uM PLX4720 and short-term response i.e. less that 120 hrs
data <- subset(data, data$conc==8)
data <- subset(data, data$Time >=0 & data$Time <=100)

#=====================================================================================================
# Plot individual clonal responses in the short-term
#=====================================================================================================

pdf(paste0(output,"/","SKMEL5 + 8uM PLX4720 cFP Assay.pdf"), width=3, height=4)
cl1 <- rainbow(length(unique(data$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~Time, data=data, type="n", ylim=c(-4,3), xlim=c(0,105), ylab="", xlab="")
for(w in 1:length(unique(data$cID))){
  temp <- subset(data, data$cID==w)
  lines(temp$Time, temp$nl2, col=cl1[w])
}
dev.off()

#=====================================================================================================
# Comparison of clonal composite vs population level response
#=====================================================================================================

d_pop <- read.csv("Population_Response_CellLines.csv", header=T, sep=",")
d_pop <- subset(d_pop,  d_pop$CellLine =="SKMEL5" & d_pop$conc ==8)
d_pop <- d_pop[(d_pop$Date %in% c("20150916", "20150831")),]
d_pop <- subset(d_pop, d_pop$Time<110)
d_pop <- d_pop[,c("CellLine", "rep", "Time", "l2", "nl2")]
d_pop$group <- "Population Response"

# Clonal composite
d_cfp <- read.csv("SKMEL5_cFP_Assay.csv", header=T, sep=",")
d_cfp <- d_cfp[(d_cfp$Date %in% c("20140808", "20140117")),]
d_cfp <- subset(d_cfp, d_cfp$Time>=0)
dn <- data.frame()
for(i in unique(d_cfp$Date)){
  s <- subset(d_cfp, d_cfp$Date==i)
  s <- aggregate(Count~Time, data=s, FUN='sum')
  s$l2 = log2(s$Count)
  s$nl2 = s$l2 - s$l2[s$Time==unique(s$Time[1])]
  s <- subset(s, s$Time<110)
  dn <- rbind(dn, s)
}
d_cfp <- dn
d_cfp$CellLine <- "SKMEL5"; d_cfp$rep <- 1
d_cfp <- d_cfp[,c("CellLine", "rep" , "Time", "l2", "nl2" )]
d_cfp$group <- "Clonal composite"

# Combine data frame into one
d <- rbind(d_pop, d_cfp)

#=====================================================================================================
# Plot for clonal composite vs population level response
#=====================================================================================================
ggplot(data = d, aes(x=d$Time, y=d$nl2)) +
  theme_bw()+ geom_smooth(aes(linetype=group, fill=group), span=.5, data=d, method = "loess", 
  size=0.9, alpha=0.4, se=TRUE, col="blue")+ 
  scale_colour_manual(values=c("black","blue"), name="") + ylim(-1., 1.)+ xlim(0,100)+
  theme(legend.position=c(0.6, 0.2), legend.position="horizontal")+ theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave("SKMEL5 Population Response vs cFP Aggregate.pdf", path=output, width=3, height=3)

#=====================================================================================================
