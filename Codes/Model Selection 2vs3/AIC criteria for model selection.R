# ===========================================================================================================
# AIC criteria for model selection
# ===========================================================================================================
source('~/Paudel_et_al_2016/Codes/AIC model fit.R', encoding = 'UTF-8')
source("~/Paudel_et_al_2016/Codes/Model Selection 2vs3/Three-state_model.R") #Three-state one dimensional model
source("~/Paudel_et_al_2016/Codes/Model Selection 2vs3/Three-state_model_triangle.R")#Three-state triangle
source("~/Paudel_et_al_2016/Codes/Model Selection 2vs3/Two-state_model_RE.R") # Two-state modelRE
source("~/Paudel_et_al_2016/Codes/Model Selection 2vs3/Two-state_model_RS.R") # Two-state modelRS
source("~/Paudel_et_al_2016/Codes/Model Selection 2vs3/Two-state_model_SE.R") # Two-state modelSE
# ===========================================================================================================
# ===========================================================================================================
summary(model3st)
summary(modeltri)
summary(model2stRE)
summary(model2stRS)
summary(model2stSE)
# 
AICc(model3st)
# AICc(model2st)
# 
AICc(model2stRS,model2stSE, model2stRE, modeltri, model3st)
# 
# ===========================================================================================================
# ===========================================================================================================
test<-ftest(modeltri,model3st)
test$p.value
