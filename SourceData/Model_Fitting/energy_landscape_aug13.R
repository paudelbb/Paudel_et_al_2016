# ===========================================================================================================
# Inferring energy landscape from the rate constants;
# ===========================================================================================================
# Different rate constants that we have: kab, kba, kbc & kcb
# How likely it is to transition to A from A; kab
# How likely it is to transition to B from A; kba
# How likely it is to transition to C from B; kcb
# How likely it is to transition to B from C; kbc
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
fig = "Supp.Fig7" # Which figure is this in paper?
output = paste0(output, "/", fig)
# ===========================================================================================================
# Calculating the energy of activation from the rate constants
# ===========================================================================================================
# Arrhenius Equation: k = A*exp(-Ea/RT)
# Solving this, we get: Ea = -K*ln(k); 
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
# ===========================================================================================================
# PLot the energy landscape based on the model fitted transition rate constants
# ===========================================================================================================
# Transition from state A to state B
n = 10
k1 = sample(log(1/rates$k_ab), n)
x1 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k1)){
  dc = 1/2*k1[i]*x1^2
  dc = data.frame(x = x1, y=dc, k=k1[i])
  dn = rbind(dn, dc)
}
df1 = dn
plot(y~x, data=df1)
# Transition from state B to state A
k2 = sample(log(1/rates$k_ba), n)
x2 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k2)){
  dc = 1/2*k2[i]*x2^2
  dc = data.frame(x = x2, y=dc, k=k2[i])
  dn = rbind(dn, dc)
}
df2 = dn
plot(y~x, data=df2)
# Transition from state B to state C
k3 = sample(log(1/rates$k_bc), n)
x1 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k3)){
  dc = 1/2*k3[i]*x1^2
  dc = data.frame(x = x1, y=dc, k=k3[i])
  dn = rbind(dn, dc)
}
df3 = dn
plot(y~x, data=df3)
# Transition from state C to state B
k4 = sample(log(1/rates$k_cb), n)
x2 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k4)){
  dc = 1/2*k4[i]*x2^2
  dc = data.frame(x = x2, y=dc, k=k4[i])
  dn = rbind(dn, dc)
}
df4 = dn
plot(y~x, data=df4)

df1
df2$x = df2$x + 21
df3$x = df3$x + 21 + 21
df4$x = df4$x + 21 + 21 + 21
# ===========================================================================================================
# ===========================================================================================================
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
plot(y~x, data=df, type="n", ylab="", xlab="", xaxt="n", yaxt="n")
for(w in unique(df$id)){
  temp <- subset(df, df$id==w)
  lines(temp$x, temp$y, col=cl1[w])
}
# 
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
plot(y~x, data=df, type="n", ylab="", xlab="")
for(w in unique(df$id)){
  temp <- subset(df, df$id==w)
  lines(temp$x, temp$y, col=cl1[w])
}

# ===========================================================================================================
dn = data.frame()
for(i in 1:length(unique(dt1$id))){
  temp1 = subset(dt1, dt1$id==unique(dt1$id)[i])
  temp2 = subset(dt2, dt2$id==unique(dt2$id)[i])
  diff = max(temp2$y) - max(temp1$y)
  temp2$y = temp2$y - diff
  dc = rbind(temp1,temp2)
  dn = rbind(dn, dc)
}

uid = seq(1:(length(dn$x)/84))
dn$id = rep(uid, each=84)
df = dn
cl1 <- rainbow(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(y~x, data=df, type="n", ylab="", xlab="", xaxt="n", yaxt="n")
for(w in unique(df$id)){
  temp <- subset(df, df$id==w)
  lines(temp$x, temp$y, col=cl1[w])
}
abline(h=0, lty=2)
plot(y~x, data=temp, type="l")
# ===========================================================================================================
dn = data.frame()
for(i in unique(df$id)){
  temp = subset(df, df$id == i)
  temp1 = subset(temp, temp$x <=10)
  temp2 = subset(temp, temp$x > 10 & temp$x <= 31)
  temp3 = subset(temp, temp$x > 31 & temp$x <= 52)
  temp4 = subset(temp, temp$x > 52)

  a = max(temp2$y) - min(temp2$y)
  b = max(temp3$y) - min(temp3$y)
  diff = b - a

  temp1$y = temp1$y - diff
  temp2$y = temp2$y - diff

  temp3a = subset(temp3, temp3$x <=34 & temp3$y <= (max(temp1$y)+10))
  temp3b = subset(temp3, temp3$x >34)
  temp3 = rbind(temp3a, temp3b)

  temp1$x = temp1$x + 23
  temp4$x = temp4$x - 1

  dc = rbind(temp1, temp3, temp4)
  dn = rbind(dn, dc)
}
plot(y~x, data=dc)

df = dn
df = df[(order(df$id)),]
cl1 <- rainbow(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(y~x, data=df, type="n", ylab="", xlab="", ylim=c(-500, 500), xaxt="n", yaxt="n")
for(w in 1:length(unique(df$id))){
  temp <- subset(df, df$id==w)
  lines(temp$x, temp$y, col = "grey")
}

id = c(1, 10, 21, 31, 41, 51, 61, 71, 81, 100)
db = df[(df$id %in% id),]
wd=1.5
plot(y~x, data=subset(db, db$id==1), ylim=c(-500, 500), xaxt="n", yaxt="n", xlab="", ylab = "", type="l", lwd=2)
points(y~x, data=subset(db, db$id==10), type="l", lty=2, lwd=wd)
points(y~x, data=subset(db, db$id==21), type="l", lty=3, lwd=wd)   
points(y~x, data=subset(db, db$id==31), type="l", lty=5, lwd=wd) 
points(y~x, data=subset(db, db$id==41), type="l", lty=6, lwd=wd)   
points(y~x, data=subset(db, db$id==51), type="l", lty=7, lwd=wd)   
points(y~x, data=subset(db, db$id==61), type="l", lty=8, lwd=wd)   
points(y~x, data=subset(db, db$id==71), type="l", lty=9, lwd=wd)   
points(y~x, data=subset(db, db$id==81), type="l", lty=9, lwd=wd)  
points(y~x, data=subset(db, db$id==100), type="l", lty=10, lwd=wd)   


# ===========================================================================================================
# ===========================================================================================================




# pdf(file = "new_energy_landscape.pdf")
# for(w in 1:length(unique(df$id))){
#   temp <- subset(df, df$id==w)
#   plot(temp$x, temp$y, col = "red", main = paste0(w), ylab="", xlab="", ylim=c(-500, 500), lwd=2,
#        xaxt="n", yaxt="n",type = "l")
# }
# dev.off()
# ===========================================================================================================
# find the activation energy for the wells by the height of the barrier
dn = data.frame()
for(i in unique(df1$k)){
  temp = subset(df1, df1$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = max(temp$y) - min(temp$y)
  dn = rbind(dn, dc)
}
sA = dn

dn = data.frame()
for(i in unique(df3$k)){
  temp = subset(df3, df3$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = max(temp$y) - min(temp$y)
  dn = rbind(dn, dc)
}
sC = dn


dn = data.frame()
for(i in unique(df2$k)){
  temp = subset(df2, df2$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = max(temp$y) - min(temp$y)
  dn = rbind(dn, dc)
}
sB = dn

dn = data.frame()
for(i in unique(df4$k)){
  temp = subset(df4, df4$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = max(temp$y) - min(temp$y)
  dn = rbind(dn, dc)
}
sD = dn
# ===========================================================================================================
# find the activation energy for the wells
dn = data.frame()
for(i in unique(df1$k)){
  temp = subset(df1, df1$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = -log(dc$k)
  dn = rbind(dn, dc)
}
sA = dn

dn = data.frame()
for(i in unique(df3$k)){
  temp = subset(df3, df3$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = -log(dc$k)
  dn = rbind(dn, dc)
}
sC = dn

dn = data.frame()
for(i in unique(df2$k)){
  temp = subset(df2, df2$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = -log(dc$k)
  dn = rbind(dn, dc)
}
sB = dn

dn = data.frame()
for(i in unique(df4$k)){
  temp = subset(df4, df4$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = -log(dc$k)
  dn = rbind(dn, dc)
}
sD = dn
# ===========================================================================================================
# ===========================================================================================================
# find the activation energy for the wells
dn = data.frame()
for(i in unique(df1$k)){
  temp = subset(df1, df1$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = -log(dc$k)
  dn = rbind(dn, dc)
}
sA = dn

dn = data.frame()
for(i in unique(df3$k)){
  temp = subset(df3, df3$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = -log(dc$k)
  dn = rbind(dn, dc)
}
sC = dn


dn = data.frame()
for(i in unique(df2$k)){
  temp = subset(df2, df2$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = -log(dc$k)
  dn = rbind(dn, dc)
}
sB = dn

dn = data.frame()
for(i in unique(df4$k)){
  temp = subset(df4, df4$k==i)
  dc = data.frame(k = unique(temp$k))
  dc$Ea = -log(dc$k)
  dn = rbind(dn, dc)
}
sD = dn


a1 = sample((sA$Ea), (n-1))
a2 = sample((sB$Ea), (n-1))
b1 = sample((sC$Ea), (n-1))
b2 = sample((sD$Ea), (n-1))

k1 = data.frame(ratio = a1/a2, state="AB")
k2 = data.frame(ratio = b1/b2, state="BC")
k = rbind(k1, k2)
boxplot(ratio~state, data=k, notch=T)
ks.test(k1$ratio, k2$ratio)
t.test(k1$ratio, k2$ratio)








