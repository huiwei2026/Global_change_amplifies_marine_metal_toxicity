################################################
#                                              #
#                Figure S5B; Table S7          #
#                                              #
################################################
############################ Field ################################## 
############## overall #################### 
data <- subset(data, Research.setting == "Field")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb   ci.ub    
# -0.1345  0.1536  -0.8755  0.3813  -0.4357  0.1666   

############## Stressor #################### 
data <- subset(data, Research.setting == "Field")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
#             estimate      se     zval    pval    ci.lb    ci.ub    
# StressorA    0.1243  0.2234   0.5566  0.5778  -0.3136   0.5623    
# StressorW   -0.2121  0.1968  -1.0779  0.2811  -0.5977   0.1735    

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("A"=1,"H"=1,"W"=1), 
                                   type="Tukey"))), test=adjusted("none"))
# Estimate Std. Error z value Pr(>|z|)  
# H - A == 0  -1.2377     0.6089  -2.033   0.0421 *
# W - A == 0  -0.3364     0.2977  -1.130   0.2584  
# W - H == 0   0.9013     0.5996   1.503   0.1328  


############################ Laboratory ##################################
############## overall #################### 
data <- subset(data, Research.setting == "Laboratory")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb    ci.ub    
# -0.0681  0.0286  -2.3841  0.0171  -0.1240  -0.0121  * 
rm(list = ls())


############## Stressor #################### 
data <- subset(data, Research.setting == "Laboratory")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
#               estimate      se     zval    pval    ci.lb    ci.ub     
# StressorA     -0.0018  0.0348  -0.0513  0.9591  -0.0699   0.0663     
# StressorA-W   -0.2451  0.1172  -2.0923  0.0364  -0.4748  -0.0155   * 
# StressorH     -0.2342  0.1276  -1.8359  0.0664  -0.4842   0.0158   . 
# StressorW     -0.1451  0.0449  -3.2348  0.0012  -0.2331  -0.0572  ** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("A"=1,"A-W"=1,"H"=1,"W"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                Estimate Std. Error z value Pr(>|z|)   
# A-W - A == 0 -0.24336    0.11876  -2.049  0.04045 * 
# H - A == 0   -0.23241    0.13222  -1.758  0.07878 . 
# W - A == 0   -0.14336    0.05510  -2.602  0.00928 **
#   H - A-W == 0  0.01095    0.17321   0.063  0.94961   
# W - A-W == 0  0.10000    0.12119   0.825  0.40927   
# W - H == 0    0.08905    0.13523   0.659  0.51018     

rm(list = ls())





















