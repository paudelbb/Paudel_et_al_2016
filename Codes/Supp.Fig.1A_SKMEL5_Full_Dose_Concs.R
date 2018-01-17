#=====================================================================================================
# Figure 1: Population Level Response for different cell lines
#=====================================================================================================
#         Population level response of BRAF-mutated melanoma cell lines to PLX4720
#         CellLines used: SKMEL5, SKMEL19, SKMEL28, WM88, WM164, WM793, A375
#                 Concentrations for idling phenotype in cell lines:
#         SKMEL5(8uM), SKMEL19(8uM), SKMEL28(16uM), WM88(8uM), WM164(32uM), WM793(16uM), A375(8uM)
#=====================================================================================================
#=====================================================================================================
library("Hmisc")
library(gplots)
require(ggplot2)
setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Supp.Fig1" # Which figure is this in paper?
output = paste0(output, "/", fig)
#=====================================================================================================
# Read the experimental data
#=====================================================================================================

data <- read.csv("Population_Response_CellLines.csv", header=T, sep=",")
cell <- "SKMEL5"
s1 <- subset(data, data$CellLine==cell & data$Date =="20150617")
s1$Date <- as.character(s1$Date)

#=====================================================================================================
#  Summarize the data and get mean+-sd
#=====================================================================================================

source("~/Paudel_et_al_2016/Codes/summarySE.R")
# d <- s1[s1$conc!=0.25, ]
d = s1
df <- summarySE(d, measurevar="nl2", groupvars=c("Time", "conc"))
d <- df

#=====================================================================================================
# Graph function
#=====================================================================================================

par(ps = 12, cex = 1, cex.axis = 1)
graph <- function(time, nl2, conc, title, sub)
{
  y_min = -1
  y_max = 6.5
  plot(  time, nl2, 
         xlim=c(0,(max(d$Time)+35)), ylim=c(y_min, y_max),
         xlab='',ylab='',
         main=title, sub=sub)
  with (
    data = d, expr = errbar(Time, nl2, nl2+se, nl2-se, add=T, pch=21, cap=.002))
  for(tc in unique(conc))
  {
    ccc <- conc[tc==conc]
    lines(time[tc==conc], nl2[tc==conc])
    text((max(d$Time)-5), nl2[tc==conc][length(ccc)], label=paste0(tc,'µM'), pos=4, cex=1.0)
  }
}

#=====================================================================================================
# Plot the full dose range in SKMEL5
#=====================================================================================================

pdf(paste0(output, "/SKMEL5_full_PLX4720_Range.pdf"), height=6, width=4)
graph(d$Time, d$nl2, d$conc, title="", sub="")
dev.off()

#=====================================================================================================
# Early rates
#=====================================================================================================
d <- s1[s1$conc!=0.25, ]
d <- d[d$conc!=0.125, ]
s1 = d

# Period 1: 15-50
dn = data.frame()
for(i in unique(s1$Date)){
  t = subset(s1, s1$Date==i)
  for(j in unique(t$conc)){
    t1 = subset(t, t$conc==j)
    for(k in unique(t1$Well)){
      t2 = subset(t1, t1$Well==k)
      s = subset(t2, t2$Time>15 & t2$Time <50)
      df = data.frame(Date = unique(s$Date), Cell = unique(s$CellLine), Treatment = unique(s$drug),
                      conc = unique(s$conc), Well = unique(s$Well),
                      rates = coef(lm(l2~Time, data=s))[2])
      dn = rbind(dn, df)
    }
  }
}
dn$period = 1
period1 = dn
#=====================================================================================================
#=====================================================================================================
# Period2 : 50-120
dn = data.frame()
for(i in unique(s1$Date)){
  t = subset(s1, s1$Date==i)
  for(j in unique(t$conc)){
    t1 = subset(t, t$conc==j)
    for(k in unique(t1$Well)){
      t2 = subset(t1, t1$Well==k)
      s = subset(t2, t2$Time>50 & t2$Time <120)
      df = data.frame(Date = unique(s$Date), Cell = unique(s$CellLine), Treatment = unique(s$drug),
                      conc = unique(s$conc), Well = unique(s$Well),
                      rates = coef(lm(l2~Time, data=s))[2])
      dn = rbind(dn, df)
    }
  }
}
dn$period = 2
period2 = dn
#=====================================================================================================
#=====================================================================================================
# Period3 : >170
dn = data.frame()
for(i in unique(s1$Date)){
  t = subset(s1, s1$Date==i)
  for(j in unique(t$conc)){
    t1 = subset(t, t$conc==j)
    for(k in unique(t1$Well)){
      t2 = subset(t1, t1$Well==k)
      s = subset(t2, t2$Time>160)
      df = data.frame(Date = unique(s$Date), Cell = unique(s$CellLine), Treatment = unique(s$drug),
                      conc = unique(s$conc), Well = unique(s$Well),
                      rates = coef(lm(l2~Time, data=s))[2])
      dn = rbind(dn, df)
    }
  }
}
dn$period = 3
period3 = dn

calcRates = rbind(period1, period2, period3)
#=====================================================================================================
#=====================================================================================================
df = summarySE(calcRates, measurevar = "rates", groupvars = c("Cell", "Treatment", "conc", "period"))

#df = df[(df$conc %in% c(0, 2, 4, 8, 16, 32)),]
df1 = subset(df, df$period==3)
df1$rates[df1$conc==0] <- df$rates[df$conc==0 & df$period==2]
df = df1
mybarcol = "gray20"
par(ps = 10, cex = 1, cex.axis = 1)
mp = barplot2(df$rates, beside = T,
              col=rainbow(nrow(df)), ylab="",
              main="", font.main=5, sub = "", col.sub = mybarcol, ylim = c(-0.025, 0.05),
              plot.ci = T, ci.l = df$rates-df$se, ci.u = df$rates+df$se,
              plot.grid = T)
group = c(as.character((df$conc)))
group = noquote(group)
labels = paste0(group)
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1, 1.1), xpd = T, cex = 1)


#=====================================================================================================
# Relative cell counts normalized to control
t = subset(s1, s1$Time<=200)
dn = data.frame()
for(i in unique(t$conc)){
  t1 = t[(t$conc %in% c(0, i)),]
  for(j in unique(t1$rep)){
    t3 = subset(t1, t1$rep==j)
    for(k in unique(t3$Time)){
      t2 = t3[(t3$Time %in% c(0, k)),]
      tc = t2$Cell.Nucleus[t2$Time == k & t2$conc == i] - t2$Cell.Nucleus[t2$Time == 0 & t2$conc == i]
      cc = t2$Cell.Nucleus[t2$Time == k & t2$conc == 0] - t2$Cell.Nucleus[t2$Time == 0 & t2$conc == 0]
      dc = data.frame(Date = unique(t2$Date), drug = unique(t2$drug), conc = max(t2$conc), 
                      Time = max(t2$Time), rel=tc/cc)
      dn = rbind(dn, dc)
    }
  }
}
dn$rel[dn$Time==0] <- 1
df = summarySE(dn, measurevar = "rel", groupvars = c("Date", "drug", "conc", "Time"))
df = subset(df, df$conc!=0)
#errbar(df$Time, df$rel, df$rel+df$se, df$rel-df$se, type="o", ylab = "", xlab = "", xlim=c(0,300))
graph <- function(time, nl2, conc, title, sub)
{
  y_min = -0.2
  y_max = 1
  plot(  time, nl2, 
         xlim=c(0,(max(df$Time)+15)), ylim=c(y_min, y_max),
         xlab='',ylab='',
         main=title, sub=sub)
  with (
    data = df, expr = errbar(Time, nl2, nl2+se, nl2-se, add=T, pch=21, cap=.002))
  for(tc in unique(conc))
  {
    ccc <- conc[tc==conc]
    lines(time[tc==conc], nl2[tc==conc])
    text((max(df$Time)-0), nl2[tc==conc][length(ccc)], label=paste0(tc,'µM'), pos=4, cex=1.0)
  }
}
graph(df$Time, df$rel, df$conc, title="", sub="")
#=====================================================================================================
t = s1[(s1$conc %in% c(0,8)),]
dn = data.frame()
for(i in unique(t$rep)){
  t1 = subset(t, t$rep==i)
  for(j in unique(t1$Time)){
    t2 = t1[(t1$Time %in% c(0, j)),]
    tc = t2$Cell.Nucleus[t2$Time == j & t2$conc == 8] - t2$Cell.Nucleus[t2$Time == 0 & t2$conc == 8]
    cc = t2$Cell.Nucleus[t2$Time == j & t2$conc == 0] - t2$Cell.Nucleus[t2$Time == 0 & t2$conc == 0]
    dc = data.frame(Date = unique(t2$Date), drug = unique(t2$drug), conc = max(t2$conc), 
                    Time = max(t2$Time), rel=tc/cc)
    dn = rbind(dn, dc)
  }
}
dn$rel[dn$Time==0] <- 1
df = summarySE(dn, measurevar = "rel", groupvars = c("Date", "drug", "conc", "Time"))

#pdf(paste0(output, "/SKMEL5_relative_to_DMSO.pdf"), height=4, width=3.5)
errbar(df$Time, df$rel, df$rel+df$se, df$rel-df$se, type="o", ylab = "", xlab = "", xlim=c(0,300))
#dev.off()

#=====================================================================================================



