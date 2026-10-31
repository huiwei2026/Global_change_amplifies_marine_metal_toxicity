library(metafor)
library(dplyr)

################################################
#                                              #
#              Figure 3; Table S9              #
#                                              #
################################################
################# Figure 3A_All taxa #################
######### Overall #########
data$Latitude_abs <- abs(data$Latitude)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix

#                 estimate      se     zval    pval    ci.lb    ci.ub     
# intrcpt        -0.2465  0.0767  -3.2129  0.0013  -0.3968  -0.0961  ** 
# Latitude_abs    0.0050  0.0021   2.4309  0.0151   0.0010   0.0090   * 


######## Acidification ########
data <- subset(data, Stressor == "A")
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
data$Latitude_abs_sq <- data$Latitude_abs^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs + Latitude_abs_sq,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                   estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt            0.3305  0.1970   1.6777  0.0934  -0.0556  0.7166  . 
# Latitude_abs      -0.0180  0.0105  -1.7205  0.0853  -0.0386  0.0025  . 
# Latitude_abs_sq    0.0002  0.0001   1.6300  0.1031  -0.0000  0.0005    


######### Warming ########
data <- subset(data, Stressor == "W")
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                 estimate      se     zval    pval    ci.lb    ci.ub      
# intrcpt        -0.6119  0.0934  -6.5484  <.0001  -0.7951  -0.4288  *** 
# Latitude_abs    0.0139  0.0026   5.3073  <.0001   0.0088   0.0191  *** 

######### Acidification-Warming ########
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Stressor == "A-W")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#              estimate       se     zval    pval     ci.lb    ci.ub    
# intrcpt      20.0024  10.0271   1.9948  0.0461    0.3497  39.6551  * 
# logdist_sq    1.9947   0.9341   2.1354  0.0327    0.1638   3.8255  * 
# logdist     -12.8018   6.1547  -2.0800  0.0375  -24.8648  -0.7387  * 
  

######### Oxygen #########
data<-subset (data, Stressor == "H")
nrow(data)
data$Latitude_abs <- abs(data$Latitude)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt         0.6257  0.6616   0.9457  0.3443  -0.6710  1.9224    
# Latitude_abs   -0.0203  0.0146  -1.3936  0.1635  -0.0489  0.0083    

################# Figure 3B #################
################ Overall ################
# Algae+Protist
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Algae+Protist")
nrow(data)
data$Latitude_abs_sq <- data$Latitude_abs^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs + Latitude_abs_sq,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                   estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt            0.5385  0.3894   1.3831  0.1666  -0.2246  1.3017    
# Latitude_abs      -0.0244  0.0201  -1.2097  0.2264  -0.0638  0.0151    
# Latitude_abs_sq    0.0003  0.0002   1.5048  0.1324  -0.0001  0.0008    

# Radiata
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Radiata")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                 estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt        -0.5634  0.3010  -1.8719  0.0612  -1.1532  0.0265  . 
# Latitude_abs    0.0152  0.0140   1.0903  0.2756  -0.0122  0.0426     

# Lophotrochozoa
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Lophotrochozoa")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                 estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt        -0.2188  0.1614  -1.3555  0.1753  -0.5351  0.0976    
# Latitude_abs    0.0040  0.0040   0.9860  0.3241  -0.0039  0.0119    


# Ecdysozoa
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Ecdysozoa")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#              estimate      se     zval    pval     ci.lb    ci.ub     
# intrcpt      -6.2151  2.0876  -2.9771  0.0029  -10.3068  -2.1235  ** 
# logdist_sq   -0.5849  0.2040  -2.8673  0.0041   -0.9847  -0.1851  ** 
# logdist       3.7981  1.3145   2.8894  0.0039    1.2217   6.3744  ** 


# Deuterostomia
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Deuterostomia")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#              estimate      se     zval    pval     ci.lb    ci.ub    
# intrcpt      -7.1058  3.2565  -2.1820  0.0291  -13.4884  -0.7231  * 
# logdist_sq   -0.6435  0.3020  -2.1305  0.0331   -1.2355  -0.0515  * 
# logdist       4.2606  1.9942   2.1365  0.0326    0.3521   8.1691  * 


######### Acidification #######
# Algae+Protist+Latitude
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Algae+Protist")
data <- subset(data, Stressor == "A")
nrow(data)
data$Latitude_abs_sq <- data$Latitude_abs^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs + Latitude_abs_sq,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                   estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt            0.7028  0.4396   1.5987  0.1099  -0.1588  1.5644    
# Latitude_abs      -0.0295  0.0233  -1.2672  0.2051  -0.0752  0.0161    
# Latitude_abs_sq    0.0004  0.0003   1.3163  0.1881  -0.0002  0.0009   



# Radiata+Latitude
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Radiata")
data <- subset(data, Stressor == "A")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt        -0.2714  0.2677  -1.0138  0.3107  -0.7960  0.2533    
# Latitude_abs    0.0063  0.0110   0.5721  0.5673  -0.0153  0.0278    


# Lophotrochozoa+Latitude
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Lophotrochozoa")
data <- subset(data, Stressor == "A")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt        -0.1011  0.2169  -0.4660  0.6412  -0.5263  0.3241    
# Latitude_abs    0.0011  0.0056   0.2023  0.8397  -0.0098  0.0121    


# Ecdysozoa+Latitude
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Ecdysozoa")
data <- subset(data, Stressor == "A")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#               estimate      se     zval    pval    ci.lb    ci.ub    
# intrcpt      -5.0284  2.1070  -2.3866  0.0170  -9.1580  -0.8988  * 
# logdist_sq   -0.5166  0.2083  -2.4802  0.0131  -0.9248  -0.1084  * 
# logdist       3.2652  1.3268   2.4609  0.0139   0.6646   5.8657  * 
  

# Deuterostomia+Latitude
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Deuterostomia")
data <- subset(data, Stressor == "A")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#              estimate       se     zval    pval     ci.lb     ci.ub    
# intrcpt      59.0298  23.1693   2.5478  0.0108   13.6188  104.4408  * 
# logdist_sq    3.7344   1.6199   2.3053  0.0211    0.5594    6.9093  * 
# logdist     -29.8234  12.2564  -2.4333  0.0150  -53.8455   -5.8012  * 

######## Warming ########
# Algae+Protist+Latitude
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Algae+Protist")
data <- subset(data, Stressor == "W")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                estimate      se     zval    pval    ci.lb   ci.ub     
# intrcpt        -0.4151  0.2416  -1.7185  0.0857  -0.8886  0.0583   . 
# Latitude_abs    0.0154  0.0052   2.9815  0.0029   0.0053  0.0255  ** 
  

# Radiata+Latitude
data$Latitude_abs <- abs(data$Latitude)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
data<-subset (data, Taxa == "Radiata")
data <- subset(data, Stressor == "W")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#              estimate       se     zval    pval     ci.lb     ci.ub    
# intrcpt      34.0224  35.7903   0.9506  0.3418  -36.1252  104.1701    
# logdist_sq    3.8902   3.8477   1.0111  0.3120   -3.6511   11.4315    
# logdist     -23.1620  23.5012  -0.9856  0.3243  -69.2235   22.8996    


# Lophotrochozoa+Latitude
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Lophotrochozoa")
data <- subset(data, Stressor == "W")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt        -0.4890  0.2627  -1.8614  0.0627  -1.0039  0.0259  . 
# Latitude_abs    0.0111  0.0062   1.7832  0.0745  -0.0011  0.0232  .     


# Ecdysozoa+Latitude
data$Latitude_abs <- abs(data$Latitude)
data<-subset (data, Taxa == "Ecdysozoa")
data <- subset(data, Stressor == "W")
nrow(data)
data$logdist<-log(data$Latitude_abs)
data$logdist_sq <- data$logdist^2
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ logdist_sq + logdist,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#              estimate      se     zval    pval    ci.lb   ci.ub    
# intrcpt      -3.5669  2.7594  -1.2926  0.1961  -8.9752  1.8414    
# logdist_sq   -0.2654  0.2821  -0.9408  0.3468  -0.8184  0.2875    
# logdist       1.8903  1.7807   1.0615  0.2885  -1.5999  5.3805    



# Deuterostomia+Latitude
data$Latitude_abs <- abs(data$Latitude)
data$Latitude_abs_sq <- data$Latitude_abs^2
data<-subset (data, Taxa == "Deuterostomia")
data <- subset(data, Stressor == "W")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Latitude_abs + Latitude_abs_sq,
                   method = "REML",
                   random = list(~1 | Paper_id / Observation_id), 
                   data = data)
rmod.mix
#                   estimate      se     zval    pval    ci.lb    ci.ub      
# intrcpt           -1.7598  0.4541  -3.8751  0.0001  -2.6499  -0.8697  *** 
# Latitude_abs       0.1406  0.0412   3.4142  0.0006   0.0599   0.2213  *** 
# Latitude_abs_sq   -0.0027  0.0008  -3.3718  0.0007  -0.0043  -0.0011  *** 











