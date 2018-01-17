ftest<-function(model1,model2){
  ifelse(model1$ssr>model2$ssr,
         {
           model.large<-model2
           model.small<-model1},
         {
           model.large<-model1
           model.small<-model2
         })
  #print(model.large)
  delta1<-model.small$ssr-model.large$ssr
  npar<-length(model.large$par)-length(model.small$par)
  ms2<-model.large$ssr/model.large$df
  df2<-model.large$df
  p.value<-1-pf(delta1/npar/ms2,npar,df2)
  f.value<-delta1/npar/ms2
  df1<-npar
  df2<-df2
  list(p.value=p.value,f.value=f.value,df1=df1,df2=df2)
  
}

#test<-ftest(modelfit,modelfit2)
#test$p.value

AICc1 <- function(model1) {
  ssr<-model1$ssr
  npar<-length(model1$par)
  ndata<-length(model1$residuals)
  ndata*log(ssr/ndata)+2*npar+2*npar*(npar+1)/(ndata-npar-1)
} 


AICc <- function(...) {
  models <- list(...) 
  nmod <- length(models)
  aic <- c()
  for(i in 1:nmod) {
    aic <- c(aic,AICc1(models[[i]]))
  }
  
  minA <- min(aic)
  delta <- aic-minA
  weights <- c()
  for(i in 1:nmod) {
    weights <- c(weights,exp(-(aic[[i]]-minA)/2))
  }
  
  weights <- weights/sum(weights)
  
  list(AIC=aic,delta=delta,weights=weights)
  
}

#AICc(modelfit)

# or 

#AICc(modelfit,modelfit2)