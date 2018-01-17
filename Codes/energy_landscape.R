# ===========================================================================================================
# Figure 4A: SKMEL5 distribution of rate constants
# ===========================================================================================================
library(gplots)
require(ggplot2)
require(FME)
require(deSolve)
library(devtools)
library(ggbiplot)
library(preprocessCore)
setwd("~/Paudel_et_al_2016/SourceData/Model_Fitting") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData/Model_Fitting"))), "Figures") 
fig = "Supp.Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)
# ===========================================================================================================
# Get the MCMC estimated parameters
# ===========================================================================================================

set.seed(12345)
files <- grep(".RData", dir(), value=T)
chain1 <- load("MCMC_Subclones_Mix_150000_Run.RData")
dim(MCMC1$pars)
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
rates <- data.frame(MCMC1$pars[,c(1,2,3,4)])


# ===========================================================================================================
# PLot the energy landscape based on the model fitted transition rate constants
# ===========================================================================================================
# ===========================================================================================================
# Transition from state A to state B
k1 = sample(log10(1/rates$k_ab), 5)
x1 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k1)){
  dc = 1/2*k1[i]*x1^2
  dc = data.frame(x = x1, y=dc, k=k1[i])
  dn = rbind(dn, dc)
}
df1 = dn
plot(y~x, data=df1, ylim=c(-250, 350))
# Transition from state B to state A
k2 = sample(log10(1/rates$k_ba), 5)
x2 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k2)){
  dc = 1/2*k2[i]*x2^2
  dc = data.frame(x = x2, y=dc, k=k2[i])
  dn = rbind(dn, dc)
}
df2 = dn
plot(y~x, data=df2, ylim=c(-250, 350))
df2$x = df2$x+20
# ===========================================================================================================
# ===========================================================================================================
# Transition from state B to state C
k3 = sample(log10(1/rates$k_bc), 5)
x1 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k3)){
  dc = 1/2*k3[i]*x1^2
  dc = data.frame(x = x1, y=dc, k=k3[i])
  dn = rbind(dn, dc)
}
df3 = dn
df3$x = df3$x + 20+20
plot(y~x, data=df3, ylim=c(-250, 5))
# Transition from state C to state B
k4 = sample(log10(1/rates$k_cb), 1000)
x2 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k4)){
  dc = 1/2*k4[i]*x2^2
  dc = data.frame(x = x2, y=dc, k=k4[i])
  dn = rbind(dn, dc)
}
df4 = dn
df4$x = df4$x +30 + 10 + 20
plot(y~x, data=df4, ylim=c(-250, 350))

# ===========================================================================================================
# Quantile normalize the Well B
# ===========================================================================================================
# dty = data.frame(b1=df2$y, b2=df3$y)
# library(preprocessCore)
# dt = data.matrix(dty)
# nd = normalize.quantiles(dt,copy=TRUE)
# nd = data.frame(nd)
# df1
# df2$y = nd$X1
# df3$y = nd$X2
# df4
# ===========================================================================================================
dn = data.frame()
for(i in 1:length(unique(df1$k))){
  temp1 = subset(df1, df1$k==unique(df1$k)[i])
  for(j in 1:length(unique(df2$k))){
    temp2 = subset(df2, df2$k==unique(df2$k)[j])
    diff = max(temp2$y) - max(temp1$y)
    temp2$y = temp2$y-diff
    dc = rbind(temp1,temp2)
    dn = rbind(dn, dc)
  }
}
uid = seq(1:(length(dn$x)/42))
dn$id = rep(uid, each=42)
# 
dt1 = dn
df = dn
cl1 <- rainbow(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(y~x, data=df, type="n", ylab="", xlab="", ylim=c(-100, 150), xaxt="n", yaxt="n")
for(w in unique(df$id)){
  temp <- subset(df, df$id==w)
  lines(temp$x, temp$y, col=cl1[w])
}
# 
# 
# # ===========================================================================================================
# # ===========================================================================================================
dn = data.frame()
for(i in 1:length(unique(df3$k))){
  temp1 = subset(df3, df3$k==unique(df3$k)[i])
  for(j in 1:length(unique(df4$k))){
    temp2 = subset(df4, df4$k==unique(df4$k)[j])
    diff = max(temp2$y) - max(temp1$y)
    temp2$y = temp2$y - diff
    dc = rbind(temp1,temp2)
    dn = rbind(dn, dc)
  }
}
uid = seq(1:(length(dn$x)/42))
dn$id = rep(uid, each=42)

dt2 = dn
df = dn
cl1 <- rainbow(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
df = dn
cl1 <- rainbow(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(y~x, data=df, type="n", ylab="", xlab="", ylim=c(-100, 300))
for(w in unique(df$id)){
  temp <- subset(df, df$id==w)
  lines(temp$x, temp$y, col=cl1[w])
}


# ===========================================================================================================
# ===========================================================================================================
# 
# dn = data.frame()
# for(i in 1:length(unique(dt1$id))){
#   temp1 = subset(dt1, dt1$id==unique(dt1$id)[i])
#   for(j in 1:length(unique(dt2$id))){
#     temp2 = subset(dt2, dt2$id==unique(dt2$id)[j])
#     diff = max(temp2$y) - max(temp1$y)
#     temp2$y = temp2$y - diff
#     dc = rbind(temp1,temp2)
#     dn = rbind(dn, dc)
#   }
# }
# uid = seq(1:(length(dn$x)/42))
# dn$id = rep(uid, each=42)
# df = dn
# cl1 <- rainbow(length(unique(df$id)))
# par(ps = 12, cex = 1, cex.axis = 1)
# plot(y~x, data=df, type="n", ylab="", xlab="", ylim=c(-200, 150), xaxt="n", yaxt="n")
# for(w in unique(df$id)){
#   temp <- subset(df, df$id==w)
#   lines(temp$x, temp$y, col=cl1[w])
# }
# abline(h=0, lty=2)
# ===========================================================================================================
# a = max(df3$y)-min(df3$y)
# b = max(df2$y)-min(df2$y)
# c = max(df1$y)-min(df1$y)
# 
# 
# dn = data.frame()
# for(i in 1:length(unique(dt1$id))){
#   temp1 = subset(dt1, dt1$id==unique(dt1$id)[i])
#   for(j in 1:length(unique(dt2$id))){
#     temp2 = subset(dt2, dt2$id==unique(dt2$id)[j])
#     diff = abs(min(temp2$y)) - abs(min(temp1$y))
#     temp1$y = temp1$y - diff
#     temp2$y = temp2$y *(b/a)
#     #temp1 = subset(temp1, temp1$x <=10)
#     temp2$x = temp2$x-20
#     dc = rbind(temp1,temp2)
#     dn = rbind(dn, dc)
#   }
# }
# uid = seq(1:(length(dn$x)/42))
# dn$id = rep(uid, each=42)
# df = dn
# cl1 <- rainbow(length(unique(df$id)))
# par(ps = 12, cex = 1, cex.axis = 1)
# plot(y~x, data=df, type="n", ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(0,200))
# for(w in unique(df$id)){
#   temp <- subset(df, df$id==w)
#   lines(temp$x, temp$y, col=cl1[w])
# }
# a = subset(df, df$x <=10)
# abline(h=min(a$y), lty=2)
# ===========================================================================================================
# ===========================================================================================================
# find the activation energy for the wells
dn = data.frame()
for(i in unique(df1$k)){
  temp = subset(df1, df1$k==i)
  dc = data.frame(k = unique(temp$k), Ea = max(temp$y))
  dn = rbind(dn, dc)
}
sA = dn

dn = data.frame()
for(i in unique(df3$k)){
  temp = subset(df3, df3$k==i)
  dc = data.frame(k = unique(temp$k), Ea = max(temp$y))
  dn = rbind(dn, dc)
}
sC = dn


dn = data.frame()
for(i in unique(df2$k)){
  temp = subset(df2, df2$k==i)
  dc = data.frame(k = unique(temp$k), Ea = max(temp$y))
  dn = rbind(dn, dc)
}
sB = dn

dn = data.frame()
for(i in unique(df4$k)){
  temp = subset(df4, df4$k==i)
  dc = data.frame(k = unique(temp$k), Ea = max(temp$y))
  dn = rbind(dn, dc)
}
sD = dn

a1 = sample(log2(sA$Ea), 900)
a2 = sample(log2(sB$Ea), 900)
b1 = sample(log2(sC$Ea), 900)
b2 = sample(log2(sD$Ea), 900)

k1 = data.frame(ratio = a1/a2, state="AB")
k2 = data.frame(ratio = b1/b2, state="BC")
k = rbind(k1, k2)
boxplot(ratio~state, data=k, notch=T)
ks.test(k1, k2)
t.test(k1$ratio, k2$ratio)






