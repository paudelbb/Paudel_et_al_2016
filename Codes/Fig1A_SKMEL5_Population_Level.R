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
# Calculate DIP rates
#=====================================================================================================

data <- read.csv("Population_Response_CellLines.csv", header=T, sep=",")
data$drug = "plx"
cnc <- 8
s1 <- subset(data, data$conc==cnc)
s1$Date <- as.character(s1$Date)
#=====================================================================================================
# calculate DIP rates
#=====================================================================================================

dn = data.frame()
for(i in unique(s1$Date)){
  temp = subset(s1, s1$Date==i)
  for(j in unique(temp$CellLine)){
    temp1 = subset(temp, temp$CellLine==j)
    for(l in unique(temp1$conc)){
      temp3 = subset(temp1, temp1$conc==l)
      for(k in unique(temp3$Well)){
        temp2 = subset(temp3, temp3$Well==k)
        s = subset(temp2, temp2$Time>24 & temp2$Time<=100)
        df = data.frame(Date = unique(s$Date), Cell = unique(s$CellLine), Treatment = unique(s$drug),
                        conc = unique(s$conc), Well = unique(s$Well), 
                        rates = coef(lm(l2~Time, data=s))[2])
        dn = rbind(dn, df)
      }
    }
  }
}
calcRates = dn
df = summarySE(calcRates, measurevar = "rates", groupvars = c("Cell", "Treatment", "conc"))
#=====================================================================================================
df = df[(order(df$rates)),]
par(ps = 10, cex = 1, cex.axis = 1)
mp = barplot2(df$rates, beside = T, 
              col=rainbow(nrow(df)), ylab="", 
              main="", font.main=5, sub = "", col.sub = mybarcol, ylim=c(-0.04, 0.04), 
              plot.ci = T, ci.l = df$rates-df$se, ci.u = df$rates+df$se, 
              plot.grid = T)
group = c(as.character((df$Cell)))
group = noquote(group)
labels = paste0(group)
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1, 1.1), xpd = T, cex = 0.7)
#=====================================================================================================
