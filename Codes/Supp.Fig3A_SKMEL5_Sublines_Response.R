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

data <- read.csv("SKMEL5 Subclones_initial_screen.csv", header=T, sep=",")
data$Date <- "20150315"
s <- subset(data,  data$Time>=0 & data$Time <100)
dn <- data.frame() # getting the responses from 24-100 hours of treatment
for(i in unique(s$cID)){
  tem <- subset(s, s$cID==i)
  tem$nl2 = tem$nl2 - tem$nl2[tem$Time==unique(tem$Time)[2]]
  tem$Time = tem$Time - tem$Time[tem$Time==unique(tem$Time)[2]]
  dn <- rbind(dn, tem)
}
s <- subset(dn, dn$Time>=0)

#=====================================================================================================
#=====================================================================================================
dn = data.frame()
for(i in unique(s$Date)){
  temp = subset(s, s$Date==i)
  for(j in unique(temp$cID)){
  t2 = subset(temp, temp$cID==j)
  t2 = subset(t2, t2$Time<=60)
  dc = data.frame(Date = unique(t2$Date), cID = unique(t2$cID), 
                drug = unique(t2$drug),rate1 = coef(lm(nl2~Time, data=t2))[2])
    dn = rbind(dn, dc)
  }
}
ds = dn[(order(dn$rate1)),]
ds$se = 0.00001
#pdf(paste0(output,"/","SKMEL5 Sublines Short-term DIP rates.pdf"), width=2.5, height=3.2)
par(ps = 10, cex = 1, cex.axis = 1)
mp = barplot2(ds$rate1, beside = T,
              col=grey.colors(nrow(ds)), ylab="",
              main="", font.main=5, sub = "", col.sub = mybarcol, 
              plot.ci = T, ci.l = ds$rate1-ds$se, ci.u = ds$rate1+ds$se, ylim = c(-0.04, 0.04),
              plot.grid = T, horiz = F, cex.axis = 0.7)
# group = c(as.character((ds$cID)))
# group = noquote(group)
# labels = paste0("SC",group)
#text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1, 1.1), xpd = T, cex = 0.2)
#dev.off()



#=====================================================================================================
# Plot the responses of single-cell-derived subclones
#=====================================================================================================

pdf(paste0(output,"/","SKMEL5 + Sublines_Initial_Screen.pdf"), width=3, height=4)
cl1 <- grey.colors(length(unique(s$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~Time, data=s, type="n", ylim=c(-1.2,1.2), xlim=c(0,70), ylab="", xlab="")
for(w in (unique(s$cID))){
  temp <- subset(s, s$cID==w)
  lines(temp$Time, temp$nl2, col=cl1[w], lwd=2)
}
dev.off()


#=====================================================================================================
# Read the data for proliferation rates in growth media vs in drug for SKMEL5 subclones
# Selected clones are: SC01, SC07, SC10, SC03, SC04, SC08
#=====================================================================================================

data <- read.csv("SKMEL5 Subclones_Growth_rates_Pre_vs_Post_Treatment.csv", header=T, sep=",")
cor1 = round(cor(data$pre, data$post, method="pearson"), 3)
cor2 = round(cor(data$pre, data$post, method="spearman"), 3)

pdf(paste0(output,"/","Proliferation_rates_pre_vs_post_treatment_Sublines.pdf"), width=3.5, height=3.5)
par(ps = 12, cex = 1, cex.axis = 1)
plot(pre~post, data=data, xlab="", ylab="", xlim=c(-0.05, 0.05), ylim=c(-0.02, 0.06))
abline(v=0, lty=2)
abline(h=0, lty=2)
dev.off()

#=====================================================================================================
#=====================================================================================================


