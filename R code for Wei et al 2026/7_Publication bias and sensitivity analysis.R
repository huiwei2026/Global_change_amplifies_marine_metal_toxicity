

if (any(data$vi <= 0)) {
  data$vi[data$vi <= 0] <- 1e-6
  summary(data$vi)
  data[data$vi < 0.001, ]
  data <- data[data$vi >= 0.001, ]
}
rmod<-rma(yi,vi,data=data,method= "REML", subset = vi>0.00001)
funnel(rmod)
rmod.tf<-trimfill(rmod, ma.common= FALSE, estimator="R0")
funnel(rmod.tf)
fsn(yi, vi, data=data, type="Rosenthal", alpha=.05)
# Fail-safe N Calculation Using the Rosenthal Approach
# 
# Observed Significance Level: <.0001
# Target Significance Level:   0.05
# 
# Fail-safe N: 208573


################################################
#                                              #
#    Sensitivity analyses: Tables S18-S50      #
#                                              #
################################################

### A. Removing one study that accounts for the largest number of experiment ###
summary(data$Paper_id)
# 133     113      77      21      89     124     126     138      58      62 
# 24      20      19      15      14      13      13      13      12      11 
# 121      39      80      41     118     142      25      61      67      69 
# 11      10      10       9       9       9       8       8       8       8 
# 140      37      38      43      55      59      66      81     102     108 
# 7       6       6       6       6       6       6       6       6       6 
# 117     132     146      12      13      86     100     122       1       5 
# 6       6       6       5       5       5       5       5       4       4 
# 14      17      19      22      24      28      40      42      53      60 
# 4       4       4       4       4       4       4       4       4       4 
# 83      84      88     105     106     107     127     133     135     139 
# 4       4       4       4       4       4       4       4       4       4 
# 144     147       6      29      30      31      33      35      45      47 
# 4       4       3       3       3       3       3       3       3       3 
# 51      71      73      76      78      79      90      95     109     112 
# 3       3       3       3       3       3       3       3       3       3 
# 128     141     145       3       4       8       9      10      15      16 
# 3       3       3       2       2       2       2       2       2       2 
# 18      20      26      32      34      36      44      46      49 (Other) 
# 2       2       2       2       2       2       2       2       2      69
############## overall #################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb    ci.ub     
# -0.0785  0.0279  -2.8162  0.0049  -0.1331  -0.0239  ** 

############## Stressor #################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
#              estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA     -0.0079  0.0346  -0.2291  0.8188  -0.0758   0.0599   
# StressorA-W   -0.2488  0.1200  -2.0733  0.0381  -0.4840  -0.0136  * 
# StressorH     -0.2921  0.1249  -2.3382  0.0194  -0.5370  -0.0473  * 
# StressorW     -0.1504  0.0431  -3.4921  0.0005  -0.2348  -0.0660  *** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("A"=1,"A-W"=1,"H"=1,"W"=1), 
                                   type="Tukey"))), test=adjusted("none"))
# Estimate Std. Error z value Pr(>|z|)   
# A-W - A == 0 -0.24088    0.12178  -1.978  0.04792  * 
# H - A == 0   -0.28422    0.12965  -2.192  0.02837  * 
# W - A == 0   -0.14250    0.05388  -2.645  0.00817  ** 
# H - A-W == 0 -0.04334    0.17324  -0.250  0.80246   
# W - A-W == 0  0.09839    0.12388   0.794  0.42709   
# W - H == 0    0.14172    0.13216   1.072  0.28357   



############## Taxa ####################
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Taxa))
colnames(rep_data) <- c("Taxa","rep")
data <- merge(data, rep_data, by = c("Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                       estimate      se     zval    pval    ci.lb    ci.ub     
# TaxaAlgae+Protist     0.1961  0.0620   3.1610  0.0016   0.0745   0.3176  ** 
# TaxaDeuterostomia    -0.2018  0.0818  -2.4664  0.0136  -0.3622  -0.0414  * 
# TaxaEcdysozoa        -0.1919  0.0604  -3.1756  0.0015  -0.3103  -0.0734  ** 
# TaxaLophotrochozoa   -0.0757  0.0401  -1.8876  0.0591  -0.1543  0.0029  . 
# TaxaRadiata          -0.2351  0.0816  -2.8833  0.0039  -0.3950  -0.0753  ** 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Algae+Protist"=1,"Deuterostomia"=1,"Ecdysozoa"=1,"Lophotrochozoa"=1,"Radiata"=1), 
                                   type="Tukey"))), test=adjusted("none"))
# Estimate Std. Error z value Pr(>|z|)    
# Deuterostomia - Algae+Protist == 0  -0.39790    0.10269  -3.875 0.000107  *** 
# Ecdysozoa - Algae+Protist == 0      -0.38793    0.08659  -4.480 7.46e-06  *** 
# Lophotrochozoa - Algae+Protist == 0 -0.27177    0.07386  -3.679 0.000234  *** 
# Radiata - Algae+Protist == 0        -0.43121    0.10246  -4.209 2.57e-05  *** 
# Ecdysozoa - Deuterostomia == 0       0.00997    0.10157   0.098 0.921771   
# Lophotrochozoa - Deuterostomia == 0  0.12613    0.09106   1.385 0.166032   
# Radiata - Deuterostomia == 0        -0.03330    0.11553  -0.288 0.773148   
# Lophotrochozoa - Ecdysozoa == 0      0.11615    0.07202   1.613 0.106806   
# Radiata - Ecdysozoa == 0            -0.04328    0.10149  -0.426 0.669812   
# Radiata - Lophotrochozoa == 0       -0.15943    0.09088  -1.754 0.079381  . 



############## Habitat ####################
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Habitat))
colnames(rep_data) <- c("Habitat","rep")
data <- merge(data, rep_data, by = c("Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                    estimate      se     zval    pval    ci.lb    ci.ub      
# HabitatFrigid       0.5028  0.1636   3.0743  0.0021   0.1823   0.8234  ** 
# HabitatTemperate   -0.0813  0.0368  -2.2062  0.0274  -0.1535   -0.0091  * 
# HabitatTropical    -0.1081  0.0397  -2.7186  0.0066  -0.1860  -0.0302  ** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Frigid"=1,"Temperate"=1,"Tropical"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                           Estimate Std. Error z value Pr(>|z|)   
# Temperate - Frigid == 0   -0.58410    0.16765  -3.484  0.00049  *** 
# Tropical - Frigid == 0    -0.61087    0.16831  -3.629  0.00028  *** 
# Tropical - Temperate == 0 -0.02677    0.05418  -0.494  0.62128   


############## Trait ####################
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Trait))
colnames(rep_data) <- c("Trait","rep")
data <- merge(data, rep_data, by = c("Trait"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trait - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                      estimate      se     zval    pval    ci.lb    ci.ub     
# TraitAbundance       0.0976  0.1161   0.8400  0.4009  -0.1301   0.3252   
# TraitAccumulation   -0.1129  0.0434  -2.6001  0.0093  -0.1980  -0.0278  ** 
# TraitBehaviour      -0.1394  0.1375  -1.0144  0.3104  -0.4088   0.1300   
# TraitGrowth         -0.0503  0.0578  -0.8700  0.3843  -0.1635   0.0630   
# TraitMetabolism     -0.0555  0.0441  -1.2581  0.2084  -0.1419   0.0309   
# TraitReproduction   -0.1507  0.1033  -1.4585  0.1447  -0.3532   0.0518   
# TraitSurvival       -0.1226  0.0892  -1.3740  0.1694  -0.2975   0.0523   

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Abundance"=1,"Accumulation"=1,"Behaviour"=1,
                                     "Growth"=1,"Metabolism"=1,"Reproduction"=1,"Survival"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                  Estimate Std. Error z value Pr(>|z|)  
# Accumulation - Abundance == 0    -0.210440   0.123048  -1.710   0.0872  . 
# Behaviour - Abundance == 0       -0.236989   0.179314  -1.322   0.1863   
# Growth - Abundance == 0          -0.147820   0.125924  -1.174   0.2404   
# Metabolism - Abundance == 0      -0.153027   0.121671  -1.258   0.2085   
# Reproduction - Abundance == 0    -0.248247   0.154355  -1.608   0.1078   
# Survival - Abundance == 0        -0.220170   0.144336  -1.525   0.1272   
# Behaviour - Accumulation == 0    -0.026548   0.142247  -0.187   0.8519   
# Growth - Accumulation == 0        0.062621   0.069394   0.902   0.3668   
# Metabolism - Accumulation == 0    0.057413   0.058186   0.987   0.3238   
# Reproduction - Accumulation == 0 -0.037806   0.109426  -0.345   0.7297   
# Survival - Accumulation == 0     -0.009730   0.096638  -0.101   0.9198   
# Growth - Behaviour == 0           0.089169   0.145428   0.613   0.5398   
# Metabolism - Behaviour == 0       0.083962   0.143285   0.586   0.5579   
# Reproduction - Behaviour == 0    -0.011258   0.170063  -0.066   0.9472   
# Survival - Behaviour == 0         0.016818   0.160505   0.105   0.9165   
# Metabolism - Growth == 0         -0.005207   0.068662  -0.076   0.9395   
# Reproduction - Growth == 0       -0.100427   0.110910  -0.905   0.3652   
# Survival - Growth == 0           -0.072351   0.102509  -0.706   0.4803   
# Reproduction - Metabolism == 0   -0.095220   0.110532  -0.861   0.3890   
# Survival - Metabolism == 0       -0.067143   0.097881  -0.686   0.4927   
# Survival - Reproduction == 0      0.028076   0.130497   0.215   0.8297   


############## Metal####################
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Metal))
colnames(rep_data) <- c("Metal","rep")
data <- merge(data, rep_data, by = c("Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Metal - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                  estimate      se     zval    pval    ci.lb    ci.ub     
# MetalArsenic     -0.1231  0.2160  -0.5698  0.5688  -0.5464   0.3003   
# MetalCadmium     -0.0425  0.0561  -0.7583  0.4483  -0.1524   0.0674   
# MetalCobalt      -0.0393  0.1219  -0.3223  0.7473  -0.2782   0.1996   
# MetalCopper      -0.1439  0.0456  -3.1569  0.0016  -0.2333  -0.0546  ** 
# MetalIron         0.2087  0.1195   1.7468  0.0807  -0.0255   0.4429  . 
# MetalLead        -0.3484  0.1089  -3.1980  0.0014  -0.5620  -0.1339  ** 
# MetalManganese   -0.0966  0.1428  -0.6764  0.4988  -0.3766   0.1833   
# MetalMercury     -0.0377  0.0748  -0.5034  0.6147  -0.1843   0.1090   
# MetalNickel      -0.0133  0.1743  -0.0768  0.9388  -0.3549   0.3282   
# MetalSilver      -0.1190  0.1207  -0.9863  0.3240  -0.3556   0.1175   
# MetalZinc         0.0355  0.0897   0.3960  0.6921  -0.1403   0.2113   
  
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Arsenic "=1,"Cadmium"=1,"Cobalt"=1,
                                     "Copper"=1,"Iron"=1,"Lead "=1,"Manganese"=1,
                                     "Mercury"=1,"Nickel"=1,"Silver"=1,"Zinc"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                            Estimate Std. Error z value Pr(>|z|)    
# Cadmium - Arsenic  == 0    0.080554   0.221336   0.364 0.715911   
# Cobalt - Arsenic  == 0     0.083782   0.247883   0.338 0.735371   
# Copper - Arsenic  == 0    -0.020843   0.220460  -0.095 0.924678   
# Iron - Arsenic  == 0       0.331756   0.246821   1.344 0.178911   
# Lead  - Arsenic  == 0     -0.225356   0.240829  -0.936 0.349401   
# Manganese - Arsenic  == 0  0.026455   0.258785   0.102 0.918575   
# Mercury - Arsenic  == 0    0.085400   0.228075   0.374 0.708078   
# Nickel - Arsenic  == 0     0.109686   0.277502   0.395 0.692650   
# Silver - Arsenic  == 0     0.004033   0.246936   0.016 0.986969   
# Zinc - Arsenic  == 0       0.158580   0.233091   0.680 0.496293   
# Cobalt - Cadmium == 0      0.003228   0.132972   0.024 0.980634   
# Copper - Cadmium == 0     -0.101397   0.070933  -1.429 0.152867   
# Iron - Cadmium == 0        0.251202   0.131964   1.904 0.056968  . 
# Lead  - Cadmium == 0      -0.305910   0.121763  -2.512 0.011993  * 
# Manganese - Cadmium == 0  -0.054099   0.152200  -0.355 0.722255   
# Mercury - Cadmium == 0     0.004846   0.093204   0.052 0.958535   
# Nickel - Cadmium == 0      0.029132   0.183043   0.159 0.873548   
# Silver - Cadmium == 0     -0.076521   0.126747  -0.604 0.546023   
# Zinc - Cadmium == 0        0.078025   0.104683   0.745 0.456058   
# Copper - Cobalt == 0      -0.104625   0.129947  -0.805 0.420742   
# Iron - Cobalt == 0         0.247974   0.170360   1.456 0.145507   
# Lead  - Cobalt == 0       -0.309138   0.163418  -1.892 0.058531  . 
# Manganese - Cobalt == 0   -0.057326   0.155888  -0.368 0.713066   
# Mercury - Cobalt == 0      0.001618   0.142397   0.011 0.990933   
# Nickel - Cobalt == 0       0.025904   0.212619   0.122 0.903031   
# Silver - Cobalt == 0      -0.079749   0.167638  -0.476 0.634273   
# Zinc - Cobalt == 0         0.074798   0.148259   0.505 0.613904   
# Iron - Copper == 0         0.352599   0.127864   2.758 0.005823  ** 
# Lead  - Copper == 0       -0.204513   0.116631  -1.754 0.079514  . 
# Manganese - Copper == 0    0.047298   0.149698   0.316 0.752036   
# Mercury - Copper == 0      0.106243   0.087507   1.214 0.224706   
# Nickel - Copper == 0       0.130529   0.180096   0.725 0.468591   
# Silver - Copper == 0       0.024876   0.128003   0.194 0.845911   
# Zinc - Copper == 0         0.179422   0.098322   1.825 0.068025  . 
# Lead  - Iron == 0         -0.557112   0.161686  -3.446 0.000570  *** 
# Manganese - Iron == 0     -0.305301   0.186088  -1.641 0.100876   
# Mercury - Iron == 0       -0.246356   0.140965  -1.748 0.080526  . 
# Nickel - Iron == 0        -0.222070   0.211277  -1.051 0.293219   
# Silver - Iron == 0        -0.327723   0.169788  -1.930 0.053584  . 
# Zinc - Iron == 0          -0.173176   0.149107  -1.161 0.245470   
# Manganese - Lead  == 0     0.251811   0.179550   1.402 0.160778   
# Mercury - Lead  == 0       0.310756   0.131919   2.356 0.018490  * 
# Nickel - Lead  == 0        0.335042   0.205502   1.630 0.103025   
# Silver - Lead  == 0        0.229389   0.162346   1.413 0.157667   
# Zinc - Lead  == 0          0.383935   0.140265   2.737 0.006196  ** 
# Mercury - Manganese == 0   0.058945   0.160292   0.368 0.713072   
# Nickel - Manganese == 0    0.083231   0.225256   0.369 0.711760   
# Silver - Manganese == 0   -0.022422   0.182393  -0.123 0.902159   
# Zinc - Manganese == 0      0.132124   0.164872   0.801 0.422914   
# Nickel - Mercury == 0      0.024286   0.189635   0.128 0.898096   
# Silver - Mercury == 0     -0.081367   0.141109  -0.577 0.564193   
# Zinc - Mercury == 0        0.073180   0.116043   0.631 0.528287   
# Silver - Nickel == 0      -0.105653   0.211920  -0.499 0.618096   
# Zinc - Nickel == 0         0.048894   0.195112   0.251 0.802129   
# Zinc - Silver == 0         0.154546   0.147142   1.050 0.293570   




############## Trophic.level #################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore          -0.1046  0.0363  -2.8802  0.0040  -0.1757  -0.0334  ** 
# Trophic.levelPredator           -0.1721  0.0505  -3.4053  0.0007  -0.2711  -0.0730  *** 
# Trophic.levelPrimary producer    0.1277  0.0600   2.1286  0.0333   0.0101   0.2452  * 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                    Estimate Std. Error z value Pr(>|z|)    
# Predator - Herbivore == 0         -0.06750    0.06174  -1.093 0.274231   
# Primary producer - Herbivore == 0  0.23223    0.07009   3.313 0.000922  *** 
# Primary producer - Predator == 0   0.29973    0.07711   3.887 1.01e-04  *** 








############## Stressor:Taxa #################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Taxa))
colnames(rep_data) <- c("Stressor","Taxa","rep")
data <- merge(data, rep_data, by = c("Stressor","Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TaxaAlgae+Protist       0.1857  0.0651   2.8533  0.0043   0.0581   0.3133  ** 
# StressorA-W:TaxaAlgae+Protist    -0.1206  0.2802  -0.4305  0.6669  -0.6697   0.4285   
# StressorW:TaxaAlgae+Protist       0.2450  0.1299   1.8865  0.0592  -0.0095   0.4996  . 
# StressorA:TaxaDeuterostomia      -0.1118  0.1014  -1.1017  0.2706  -0.3106   0.0871   
# StressorA-W:TaxaDeuterostomia    -0.1787  0.2075  -0.8614  0.3890  -0.5853   0.2279   
# StressorW:TaxaDeuterostomia      -0.3356  0.1196  -2.8068  0.0050  -0.5700  -0.1013  ** 
# StressorA:TaxaEcdysozoa           0.0582  0.0869   0.6700  0.5029  -0.1121   0.2285   
# StressorH:TaxaEcdysozoa          -0.4313  0.1873  -2.3021  0.0213  -0.7985  -0.0641  * 
# StressorW:TaxaEcdysozoa          -0.3316  0.0814  -4.0750  <.0001  -0.4911  -0.1721  *** 
# StressorA:TaxaLophotrochozoa     -0.0791  0.0475  -1.6666  0.0956  -0.1721   0.0139  . 
# StressorA-W:TaxaLophotrochozoa   -0.2260  0.1826  -1.2378  0.2158  -0.5839   0.1319   
# StressorH:TaxaLophotrochozoa     -0.2226  0.1690  -1.3169  0.1879  -0.5539  0.1087   
# StressorW:TaxaLophotrochozoa     -0.0366  0.0623  -0.5870  0.5572  -0.1588   0.0856   
# StressorA:TaxaRadiata            -0.1352  0.1195  -1.1312  0.2580  -0.3694   0.0990   
# StressorW:TaxaRadiata            -0.3133  0.1031  -3.0377  0.0024  -0.5154  -0.1111  ** 




############## Stressor:Habitat #################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Habitat))
colnames(rep_data) <- c("Stressor","Habitat","rep")
data <- merge(data, rep_data, by = c("Stressor","Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:HabitatPolar            0.4841   0.2388   2.0276   0.0426   0.0162   0.9521  * 
# StressorW:HabitatPolar            0.5173   0.1899   2.7246   0.0064   0.1452   0.8894  ** 
# StressorA:HabitatTemperate       -0.0545   0.0405  -1.3451   0.1786  -0.1339   0.0249   
# StressorA-W:HabitatTemperate     -0.0873   0.1486  -0.5875   0.5569  -0.3785   0.2039   
# StressorH:HabitatTemperate       -0.3515   0.1226  -2.8664   0.0042  -0.5918  -0.1112  ** 
# StressorW:HabitatTemperate       -0.0554   0.0584  -0.9483   0.3430  -0.1698   0.0591   
# StressorA:HabitatTropical         0.0480   0.0471   1.0183   0.3085  -0.0444   0.1403   
# StressorA-W:HabitatTropical      -0.4808   0.1839  -2.6140   0.0089  -0.8413  -0.1203  ** 
# StressorW:HabitatTropical        -0.2923   0.0536  -5.4571   <.0001  -0.3972  -0.1873  *** 

############## Stressor:Habitat:Taxa#################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Habitat,data$Taxa))
colnames(rep_data) <- c("Stressor","Habitat","Taxa","rep")
data <- merge(data, rep_data, by = c("Stressor","Habitat","Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Habitat:Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                     estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:HabitatPolar:TaxaAlgae+Protist             0.4838   0.2397   2.0178   0.0436   0.0139   0.9537  * 
# StressorW:HabitatPolar:TaxaAlgae+Protist             0.7032   0.2682   2.6216   0.0088   0.1775   1.2289  ** 
# StressorA:HabitatTemperate:TaxaAlgae+Protist         0.1174   0.0964   1.2186   0.2230  -0.0715   0.3063   
# StressorW:HabitatTemperate:TaxaAlgae+Protist         0.2724   0.2082   1.3083   0.1908  -0.1357   0.6804   
# StressorA:HabitatTropical:TaxaAlgae+Protist          0.1927   0.0909   2.1198   0.0340   0.0145   0.3708  * 
# StressorA-W:HabitatTropical:TaxaAlgae+Protist       -0.1807   0.2782  -0.6497   0.5159  -0.7260   0.3646   
# StressorW:HabitatTropical:TaxaAlgae+Protist         -0.0373   0.2006  -0.1861   0.8523  -0.4306   0.3559   
# StressorA:HabitatTemperate:TaxaDeuterostomia        -0.1357   0.0998  -1.3592   0.1741  -0.3313   0.0600   
# StressorA-W:HabitatTemperate:TaxaDeuterostomia      -0.1813   0.2049  -0.8846   0.3764  -0.5829   0.2204   
# StressorW:HabitatTemperate:TaxaDeuterostomia        -0.3418   0.1487  -2.2990   0.0215  -0.6332  -0.0504  * 
# StressorW:HabitatTropical:TaxaDeuterostomia         -0.3262   0.1915  -1.7031   0.0885  -0.7015   0.0492  . 
# StressorA:HabitatTemperate:TaxaEcdysozoa            -0.1845   0.2379  -0.7757   0.4380  -0.6508   0.2817   
# StressorH:HabitatTemperate:TaxaEcdysozoa            -0.4292   0.1837  -2.3360   0.0195  -0.7892  -0.0691  * 
# StressorW:HabitatTemperate:TaxaEcdysozoa            -0.2430   0.1931  -1.2586   0.2082  -0.6214   0.1354   
# StressorA:HabitatTropical:TaxaEcdysozoa              0.0952   0.0907   1.0494   0.2940  -0.0826   0.2731   
# StressorW:HabitatTropical:TaxaEcdysozoa             -0.3486   0.0872  -3.9983   <.0001  -0.5195  -0.1777  *** 
# StressorW:HabitatPolar:TaxaLophotrochozoa            0.3109   0.2793   1.1130   0.2657  -0.2366   0.8584   
# StressorA:HabitatTemperate:TaxaLophotrochozoa       -0.0789   0.0532  -1.4825   0.1382  -0.1831   0.0254   
# StressorA-W:HabitatTemperate:TaxaLophotrochozoa     -0.0814   0.2136  -0.3812   0.7031  -0.5002   0.3373   
# StressorH:HabitatTemperate:TaxaLophotrochozoa       -0.3601   0.1876  -1.9196   0.0549  -0.7277   0.0076  . 
# StressorW:HabitatTemperate:TaxaLophotrochozoa       -0.0124   0.0707  -0.1751   0.8610  -0.1509   0.1261   
# StressorA:HabitatTropical:TaxaLophotrochozoa        -0.0750   0.0945  -0.7930   0.4278  -0.2602   0.1103   
# StressorA-W:HabitatTropical:TaxaLophotrochozoa      -0.5888   0.3349  -1.7582   0.0787  -1.2451   0.0676  . 
# StressorW:HabitatTropical:TaxaLophotrochozoa        -0.1952   0.1365  -1.4296   0.1528  -0.4628   0.0724   
# StressorA:HabitatTropical:TaxaRadiata               -0.1351   0.1173  -1.1519   0.2494  -0.3649   0.0947   
# StressorW:HabitatTropical:TaxaRadiata               -0.3139   0.1013  -3.1000   0.0019  -0.5124  -0.1154  ** 


############## Stressor:Trophic.level #################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Trophic.level))
colnames(rep_data) <- c("Stressor","Trophic.level","rep")
data <- merge(data, rep_data, by = c("Stressor","Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trophic.level- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                           estimate      se     zval    pval    ci.lb    ci.ub     
# StressorA:Trophic.levelHerbivore          -0.0343  0.0485  -0.7076  0.4792  -0.1294   0.0608   
# StressorA-W:Trophic.levelHerbivore        -0.2654  0.1923  -1.3799  0.1676  -0.6423   0.1116   
# StressorH:Trophic.levelHerbivore          -0.1861  0.1648  -1.1289  0.2590  -0.5091   0.1370   
# StressorW:Trophic.levelHerbivore          -0.1785  0.0553  -3.2294  0.0012  -0.2869  -0.0702  ** 
# StressorA:Trophic.levelPredator           -0.0954  0.0669  -1.4266  0.1537  -0.2265   0.0357   
# StressorA-W:Trophic.levelPredator         -0.2875  0.1536  -1.8721  0.0612  -0.5885   0.0135  . 
# StressorH:Trophic.levelPredator           -0.4278  0.1872  -2.2848  0.0223  -0.7948  -0.0608  * 
# StressorW:Trophic.levelPredator           -0.2076  0.0736  -2.8206  0.0048  -0.3519  -0.0633  ** 
# StressorA:Trophic.levelPrimary producer    0.1273  0.0676   1.8842  0.0595  -0.0051   0.2598  . 
# StressorW:Trophic.levelPrimary producer    0.1289  0.1239   1.0405  0.2981  -0.1139   0.3718   


############## Habitat:Trophic.level #################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Habitat,data$Trophic.level))
colnames(rep_data) <- c("Habitat","Trophic.level","rep")
data <- merge(data, rep_data, by = c("Habitat","Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Habitat:Trophic.level- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                  estimate      se     zval    pval    ci.lb 
# HabitatFrigid:Trophic.levelHerbivore              0.3150  0.2986   1.0547  0.2916  -0.2704 
# HabitatTemperate:Trophic.levelHerbivore          -0.1138  0.0469  -2.4258  0.0153  -0.2057 
# HabitatTropical:Trophic.levelHerbivore           -0.1029  0.0542  -1.9008  0.0573  -0.2091 
# HabitatTemperate:Trophic.levelPredator           -0.1090  0.0668  -1.6305  0.1030  -0.2400 
# HabitatTropical:Trophic.levelPredator            -0.2479  0.0728  -3.4037  0.0007  -0.3906 
# HabitatFrigid:Trophic.levelPrimary producer       0.5792  0.1901   3.0465  0.0023  0.2066 
# HabitatTemperate:Trophic.levelPrimary producer    0.1148  0.0962   1.1930  0.2329  -0.0738 
# HabitatTropical:Trophic.levelPrimary producer     0.0508  0.0802   0.6342  0.5259  -0.1063 
#                                                    ci.ub      
# HabitatFrigid:Trophic.levelHerbivore             0.9003   
# HabitatTemperate:Trophic.levelHerbivore          -0.0219  * 
# HabitatTropical:Trophic.levelHerbivore          0.0032  . 
#   HabitatTemperate:Trophic.levelPredator           0.0220   
# HabitatTropical:Trophic.levelPredator           -0.1051  *** 
#   HabitatFrigid:Trophic.levelPrimary producer      0.9519  ** 
# HabitatTemperate:Trophic.levelPrimary producer   0.3035   
# HabitatTropical:Trophic.levelPrimary producer    0.2079   


############## Trophic.level:Metal #################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level,data$Metal))
colnames(rep_data) <- c("Trophic.level","Metal","rep")
data <- merge(data, rep_data, by = c("Trophic.level","Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level:Metal- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                               estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore:MetalArsenic          -0.0721  0.2633  -0.2736  0.7844  -0.5882   0.4441   
# Trophic.levelHerbivore:MetalCadmium          -0.1319  0.0659  -2.0009  0.0454  -0.2610  -0.0027  * 
# Trophic.levelPredator:MetalCadmium           -0.0436  0.1259  -0.3464  0.7290  -0.2903   0.2031   
# Trophic.levelPrimary producer:MetalCadmium    0.3788  0.1431   2.6464  0.0081   0.0983   0.6594  ** 
# Trophic.levelHerbivore:MetalCobalt           -0.0419  0.1987  -0.2111  0.8328  -0.4313   0.3474   
# Trophic.levelPredator:MetalCobalt            -0.0879  0.1617  -0.5437  0.5867  -0.4049   0.2290   
# Trophic.levelPrimary producer:MetalCobalt     0.1713  0.3161   0.5419  0.5879  -0.4483   0.7910   
# Trophic.levelHerbivore:MetalCopper           -0.1725  0.0598  -2.8858  0.0039  -0.2897  -0.0553  ** 
# Trophic.levelPredator:MetalCopper            -0.2147  0.0846  -2.5372  0.0112  -0.3805  -0.0488  * 
# Trophic.levelPrimary producer:MetalCopper    -0.0223  0.0894  -0.2493  0.8031  -0.1975   0.1529   
# Trophic.levelPrimary producer:MetalIron       0.2015  0.1133   1.7765  0.0756  -0.0208   0.4239  . 
# Trophic.levelHerbivore:MetalLead             -0.1659  0.1174  -1.4133  0.1576  -0.3959   0.0642   
# Trophic.levelPredator:MetalLead              -1.1299  0.2467  -4.5804  <.0001  -1.6133  -0.6464  *** 
# Trophic.levelHerbivore:MetalManganese        -0.0125  0.2923  -0.0426  0.9660  -0.5853   0.5604   
# Trophic.levelPredator:MetalManganese         -0.1427  0.1646  -0.8668  0.3860  -0.4652   0.1799   
# Trophic.levelHerbivore:MetalMercury           0.0439  0.0904   0.4862  0.6268  -0.1332   0.2211   
# Trophic.levelPredator:MetalMercury           -0.1571  0.1140  -1.3786  0.1680  -0.3805   0.0663   
# Trophic.levelHerbivore:MetalNickel           -0.0126  0.1733  -0.0728  0.9420  -0.3522   0.3270   
# Trophic.levelHerbivore:MetalSilver           -0.1530  0.1464  -1.0456  0.2957  -0.4399   0.1338   
# Trophic.levelPredator:MetalSilver            -0.1221  0.2057  -0.5935  0.5528  -0.5252   0.2811   
# Trophic.levelHerbivore:MetalZinc              0.0261  0.1466   0.1778  0.8589  -0.2614   0.3135   
# Trophic.levelPredator:MetalZinc              -0.0449  0.1403  -0.3197  0.7492  -0.3199   0.2301   
# Trophic.levelPrimary producer:MetalZinc       0.1883  0.1752   1.0743  0.2827  -0.1552   0.5317   




############## Stressor:Metal ####################
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Metal))
colnames(rep_data) <- c("Stressor","Metal","rep")
data <- merge(data, rep_data, by = c("Stressor","Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Metal - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                            estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:MetalArsenic     -0.1259  0.2127  -0.5917  0.5540  -0.5429   0.2911   
# StressorA:MetalCadmium      0.0223  0.0670   0.3331  0.7390  -0.1090   0.1537   
# StressorA-W:MetalCadmium   -0.2255  0.2590  -0.8705  0.3840  -0.7331   0.2821   
# StressorH:MetalCadmium     -0.3179  0.1505  -2.1130  0.0346  -0.6128  -0.0230  * 
# StressorW:MetalCadmium     -0.0307  0.1201  -0.2558  0.7981  -0.2662   0.2047   
# StressorA:MetalCobalt      -0.0744  0.1412  -0.5268  0.5983  -0.3511   0.2023   
# StressorW:MetalCobalt       0.2684  0.2009   1.3364  0.1814  -0.1252   0.6621   
# StressorA:MetalCopper      -0.0984  0.0577  -1.7049  0.0882  -0.2115   0.0147  . 
# StressorA-W:MetalCopper    -0.2862  0.2870  -0.9974  0.3186  -0.8487   0.2762   
# StressorW:MetalCopper      -0.2066  0.0694  -2.9761  0.0029  -0.3426  -0.0705  ** 
# StressorA:MetalIron         0.1376  0.1435   0.9590  0.3375  -0.1436   0.4189   
# StressorW:MetalIron         0.3352  0.1985   1.6893  0.0912  -0.0537   0.7242  . 
# StressorA:MetalLead        -0.1292  0.1810  -0.7135  0.4756  -0.4840   0.2256   
# StressorW:MetalLead        -0.4603  0.1327  -3.4687  0.0005  -0.7204  -0.2002  *** 
# StressorA:MetalManganese    0.1402  0.1853   0.7563  0.4495  -0.2231   0.5034   
# StressorW:MetalManganese    0.0852  0.2345   0.3634  0.7163  -0.3744   0.5448   
# StressorA:MetalMercury      0.1597  0.1008   1.5849  0.1130  -0.0378   0.3572   
# StressorA-W:MetalMercury   -0.1977  0.1835  -1.0775  0.2812  -0.5573   0.1619   
# StressorW:MetalMercury     -0.1993  0.0996  -2.0020  0.0453  -0.3945  -0.0042  * 
# StressorA:MetalNickel       0.2015  0.2496   0.8074  0.4195  -0.2877   0.6908   
# StressorW:MetalNickel      -0.1983  0.2288  -0.8666  0.3862  -0.6467   0.2502   
# StressorA:MetalSilver      -0.0946  0.1339  -0.7011  0.4832  -0.3589   0.1698   
# StressorW:MetalSilver      -0.1029  0.2587  -0.3979  0.6907  -0.6101   0.4042   
# StressorA:MetalZinc         0.0245  0.0975   0.2510  0.8018  -0.1667   0.2157   
# StressorW:MetalZinc         0.0098  0.2610   0.0374  0.9702  -0.5019   0.5214   




############## Stressor:Trait ####################
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Trait))
colnames(rep_data) <- c("Stressor","Trait","rep")
data <- merge(data, rep_data, by = c("Stressor","Trait"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trait - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                  estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TraitAbundance         0.2158  0.1724   1.2515  0.2107  -0.1221   0.5537   
# StressorW:TraitAbundance         0.0064  0.1543   0.0414  0.9670  -0.2960   0.3087   
# StressorA:TraitAccumulation     -0.0327  0.0499  -0.6547  0.5126  -0.1304   0.0651   
# StressorA-W:TraitAccumulation   -0.2136  0.1704  -1.2532  0.2101  -0.5477   0.1205   
# StressorH:TraitAccumulation     -0.1562  0.1636  -0.9546  0.3398  -0.4769   0.1645   
# StressorW:TraitAccumulation     -0.2570  0.0767  -3.3515  0.0008  -0.4073  -0.1067  *** 
# StressorA:TraitBehaviour         0.0078  0.2045   0.0380  0.9697  -0.3930   0.4086   
# StressorW:TraitBehaviour        -0.1661  0.2063  -0.8050  0.4208  -0.5704   0.2382   
# StressorA:TraitGrowth            0.0179  0.0674   0.2651  0.7909  -0.1142   0.1500   
# StressorW:TraitGrowth           -0.1410  0.1023  -1.3779  0.1682  -0.3415   0.0596   
# StressorA:TraitMetabolism       -0.0410  0.0572  -0.7175  0.4730  -0.1531   0.0710   
# StressorA-W:TraitMetabolism     -0.5583  0.2286  -2.4426  0.0146  -1.0063  -0.1103  * 
# StressorH:TraitMetabolism       -0.4126  0.2138  -1.9299  0.0536  -0.8317   0.0064  . 
# StressorW:TraitMetabolism       -0.0138  0.0677  -0.2041  0.8383  -0.1466   0.1190   
# StressorA:TraitReproduction      0.1616  0.1398   1.1566  0.2474  -0.1123   0.4355   
# StressorW:TraitReproduction     -0.4888  0.1486  -3.2891  0.0010  -0.7801  -0.1975  ** 
# StressorA:TraitSurvival          0.0623  0.1472   0.4229  0.6723  -0.2263   0.3508   
# StressorW:TraitSurvival         -0.2656  0.1117  -2.3783  0.0174  -0.4844  -0.0467  * 


############## Stressor:Trait:Habitat  #################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Trait,data$Habitat))
colnames(rep_data) <- c("Stressor","Trait","Habitat","rep")
data <- merge(data, rep_data, by = c("Stressor","Trait","Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trait:Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:TraitMetabolism:HabitatPolar              0.6650   0.2336   2.8472   0.0044   0.2072   1.1228  ** 
# StressorA:TraitAbundance:HabitatTemperate           0.0301   0.2211   0.1360   0.8918  -0.4032   0.4633   
# StressorA:TraitAccumulation:HabitatTemperate       -0.0579   0.0563  -1.0277   0.3041  -0.1682   0.0525   
# StressorA-W:TraitAccumulation:HabitatTemperate     -0.1299   0.1807  -0.7188   0.4723  -0.4841   0.2243   
# StressorH:TraitAccumulation:HabitatTemperate       -0.3229   0.1627  -1.9847   0.0472  -0.6417  -0.0040  * 
# StressorW:TraitAccumulation:HabitatTemperate       -0.0852   0.0919  -0.9266   0.3541  -0.2654   0.0950   
# StressorA:TraitBehaviour:HabitatTemperate          -0.0465   0.2669  -0.1742   0.8617  -0.5696   0.4766   
# StressorA:TraitGrowth:HabitatTemperate             -0.0195   0.1018  -0.1920   0.8477  -0.2190   0.1799   
# StressorW:TraitGrowth:HabitatTemperate             -0.3356   0.1478  -2.2704   0.0232  -0.6254  -0.0459  * 
# StressorA:TraitMetabolism:HabitatTemperate         -0.0588   0.0715  -0.8227   0.4107  -0.1989   0.0813   
# StressorH:TraitMetabolism:HabitatTemperate         -0.3131   0.2251  -1.3904   0.1644  -0.7543   0.1282   
# StressorW:TraitMetabolism:HabitatTemperate          0.1681   0.0995   1.6892   0.0912  -0.0270   0.3632  . 
# StressorA:TraitSurvival:HabitatTemperate           -0.2213   0.2721  -0.8131   0.4162  -0.7546   0.3121   
# StressorW:TraitSurvival:HabitatTemperate           -0.2294   0.1728  -1.3274   0.1844  -0.5681   0.1093   
# StressorA:TraitAbundance:HabitatTropical           -0.0963   0.2850  -0.3380   0.7354  -0.6549   0.4622   
# StressorW:TraitAbundance:HabitatTropical           -0.1278   0.1553  -0.8228   0.4106  -0.4322   0.1766   
# StressorA:TraitAccumulation:HabitatTropical         0.0464   0.0832   0.5576   0.5771  -0.1166   0.2094   
# StressorW:TraitAccumulation:HabitatTropical        -0.5510   0.1266  -4.3540   <.0001  -0.7990  -0.3030  *** 
# StressorA:TraitBehaviour:HabitatTropical            0.0747   0.2963   0.2520   0.8011  -0.5061   0.6554   
# StressorW:TraitBehaviour:HabitatTropical           -0.2404   0.2350  -1.0232   0.3062  -0.7009   0.2201   
# StressorA:TraitGrowth:HabitatTropical               0.0542   0.0853   0.6363   0.5246  -0.1128   0.2213   
# StressorW:TraitGrowth:HabitatTropical              -0.0391   0.1296  -0.3015   0.7630  -0.2931   0.2150   
# StressorA:TraitMetabolism:HabitatTropical          -0.0302   0.0856  -0.3528   0.7242  -0.1979   0.1375   
# StressorA-W:TraitMetabolism:HabitatTropical        -0.6381   0.2511  -2.5419   0.0110  -1.1302  -0.1461  * 
# StressorW:TraitMetabolism:HabitatTropical          -0.2558   0.0910  -2.8115   0.0049  -0.4342  -0.0775  ** 
# StressorA:TraitReproduction:HabitatTropical         0.2410   0.1409   1.7100   0.0873  -0.0352   0.5172  . 
# StressorW:TraitReproduction:HabitatTropical        -0.5169   0.1435  -3.6018   0.0003  -0.7982  -0.2356  *** 
# StressorA:TraitSurvival:HabitatTropical             0.1790   0.1662   1.0770   0.2815  -0.1468   0.5048   
# StressorW:TraitSurvival:HabitatTropical            -0.3068   0.1326  -2.3135   0.0207  -0.5668  -0.0469  * 


############## Trophic.level-Accumulation #################### 
data<-subset (data,data$Paper_id!=133)
nrow(data)
data <- subset(data, Trait == "Accumulation")
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub    
# Trophic.levelHerbivore          -0.0948  0.0637  -1.4875  0.1369  -0.2196   0.0301   
# Trophic.levelPredator           -0.1969  0.0821  -2.3984  0.0165  -0.3578  -0.0360  * 
# Trophic.levelPrimary producer    0.4938  0.2011   2.4547  0.0141   0.0995   0.8880  * 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                     Estimate Std. Error z value Pr(>|z|)   
# Predator - Herbivore == 0          -0.1021     0.1039  -0.983  0.32572   
# Primary producer - Herbivore == 0   0.5885     0.2110   2.789  0.00528  ** 
# Primary producer - Predator == 0    0.6907     0.2173   3.179  0.00148  ** 


############## Trophic.level-Abundance #################### 
data <- subset(data, Trait == "Abundance")
nrow(data)
data<-subset (data,data$Paper_id!=133)
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Trophic.levelHerbivore          -0.1820  0.4154  -0.4381  0.6613  -0.9961  0.6321   
# Trophic.levelPredator           -0.0873  0.3618  -0.2412  0.8094  -0.7964  0.6219   
# Trophic.levelPrimary producer    0.2546  0.2869   0.8875  0.3748  -0.3077  0.8169   

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                     Estimate Std. Error z value Pr(>|z|)
# Predator - Herbivore == 0          0.09472    0.55085   0.172    0.863   
# Primary producer - Herbivore == 0  0.43659    0.50479   0.865    0.387   
# Primary producer - Predator == 0   0.34187    0.46175   0.740    0.459   











### B. Removing one trait of experiment ###

############## overall #################### 
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb    ci.ub     
# -0.0817  0.0284  -2.8783  0.0040  -0.1373  -0.0261  ** 

rm(list = ls())



############## Stressor #################### 
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
#               estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA     -0.0114  0.0347  -0.3299  0.7414  -0.0794   0.0565   
# StressorA-W   -0.2556  0.1122  -2.2784  0.0227  -0.4754  -0.0357  * 
# StressorH     -0.2991  0.1203  -2.4860  0.0129  -0.5349  -0.0633  * 
# StressorW     -0.1527  0.0436  -3.5006  0.0005  -0.2382  -0.0672  *** 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("A"=1,"A-W"=1,"H"=1,"W"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#               Estimate Std. Error z value Pr(>|z|)   
# A-W - A == 0 -0.24412    0.11320  -2.157  0.03104  * 
# H - A == 0   -0.28767    0.12522  -2.297  0.02160  * 
# W - A == 0   -0.14125    0.05366  -2.632  0.00848  ** 
# H - A-W == 0 -0.04355    0.16449  -0.265  0.79120   
# W - A-W == 0  0.10287    0.11550   0.891  0.37309   
# W - H == 0    0.14642    0.12798   1.144  0.25259   






############## Taxa ####################
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Taxa))
colnames(rep_data) <- c("Taxa","rep")
data <- merge(data, rep_data, by = c("Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                      estimate      se     zval    pval    ci.lb    ci.ub      
# TaxaAlgae+Protist     0.1592  0.0588   2.7095  0.0067   0.0441   0.2744  ** 
# TaxaDeuterostomia    -0.2112  0.0805  -2.6226  0.0087  -0.3690  -0.0534  ** 
# TaxaEcdysozoa        -0.2014  0.0612  -3.2917  0.0010  -0.3213  -0.0815  *** 
# TaxaLophotrochozoa   -0.0609  0.0390  -1.5621  0.1183  -0.1373   0.0155   
# TaxaRadiata          -0.2496  0.0829  -3.0105  0.0026  -0.4122  -0.0871  ** 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Algae+Protist"=1,"Deuterostomia"=1,"Ecdysozoa"=1,"Lophotrochozoa"=1,"Radiata"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                      Estimate Std. Error z value Pr(>|z|)    
# Deuterostomia - Algae+Protist == 0  -0.370434   0.099684  -3.716  0.00020  *** 
# Ecdysozoa - Algae+Protist == 0      -0.360648   0.084791  -4.253 2.11e-05  *** 
# Lophotrochozoa - Algae+Protist == 0 -0.220137   0.068734  -3.203  0.00136  ** 
# Radiata - Algae+Protist == 0        -0.408886   0.101639  -4.023 5.75e-05  *** 
# Ecdysozoa - Deuterostomia == 0       0.009786   0.100866   0.097  0.92271   
# Lophotrochozoa - Deuterostomia == 0  0.150296   0.089346   1.682  0.09253  . 
# Radiata - Deuterostomia == 0        -0.038452   0.115588  -0.333  0.73939   
# Lophotrochozoa - Ecdysozoa == 0      0.140511   0.071851   1.956  0.05051  . 
# Radiata - Ecdysozoa == 0            -0.048238   0.103052  -0.468  0.63972   
# Radiata - Lophotrochozoa == 0       -0.188749   0.091628  -2.060  0.03940  * 

rm(list = ls())


############## Habitat ####################
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Habitat))
colnames(rep_data) <- c("Habitat","rep")
data <- merge(data, rep_data, by = c("Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                     estimate      se     zval    pval    ci.lb    ci.ub     
# HabitatFrigid       0.2999  0.1773   1.6913  0.0908  -0.0476   0.6474  . 
# HabitatTemperate   -0.0762  0.0380  -2.0054  0.0449  -0.1507   -0.0017  * 
# HabitatTropical    -0.1098  0.0425  -2.5831  0.0098  -0.1931  -0.0265  ** 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Frigid"=1,"Temperate"=1,"Tropical"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                           Estimate Std. Error z value Pr(>|z|)  
# Temperate - Frigid == 0   -0.37612    0.18133  -2.074   0.0381  * 
# Tropical - Frigid == 0    -0.40967    0.18234  -2.247   0.0247  * 
# Tropical - Temperate == 0 -0.03355    0.05698  -0.589   0.5560   

rm(list = ls())


############## Trait ####################
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Trait))
colnames(rep_data) <- c("Trait","rep")
data <- merge(data, rep_data, by = c("Trait"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trait - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                     estimate      se     zval    pval    ci.lb    ci.ub    
# TraitAccumulation   -0.0964  0.0411  -2.3438  0.0191  -0.1770  -0.0158  * 
# TraitBehaviour      -0.1333  0.1276  -1.0446  0.2962  -0.3835   0.1168   
# TraitGrowth         -0.0506  0.0529  -0.9577  0.3382  -0.1543   0.0530   
# TraitMetabolism     -0.0633  0.0414  -1.5275  0.1266  -0.1445   0.0179   
# TraitReproduction   -0.1510  0.0948  -1.5924  0.1113  -0.3368   0.0348   
# TraitSurvival       -0.1172  0.0838  -1.3981  0.1621  -0.2814   0.0471   

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Accumulation"=1,"Behaviour"=1,
                                     "Growth"=1,"Metabolism"=1,"Reproduction"=1,"Survival"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                   Estimate Std. Error z value Pr(>|z|)
# Behaviour - Accumulation == 0    -0.03691    0.13142  -0.281    0.779   
# Growth - Accumulation == 0        0.04577    0.06189   0.740    0.460   
# Metabolism - Accumulation == 0    0.03310    0.05257   0.630    0.529   
# Reproduction - Accumulation == 0 -0.05454    0.09960  -0.548    0.584   
# Survival - Accumulation == 0     -0.02074    0.08966  -0.231    0.817   
# Growth - Behaviour == 0           0.08268    0.13372   0.618    0.536   
# Metabolism - Behaviour == 0       0.07002    0.13254   0.528    0.597   
# Reproduction - Behaviour == 0    -0.01763    0.15623  -0.113    0.910   
# Survival - Behaviour == 0         0.01617    0.14805   0.109    0.913   
# Metabolism - Growth == 0         -0.01267    0.06123  -0.207    0.836   
# Reproduction - Growth == 0       -0.10031    0.09975  -1.006    0.315   
# Survival - Growth == 0           -0.06651    0.09409  -0.707    0.480   
# Reproduction - Metabolism == 0   -0.08765    0.10069  -0.870    0.384   
# Survival - Metabolism == 0       -0.05385    0.09089  -0.592    0.554   
# Survival - Reproduction == 0      0.03380    0.11896   0.284    0.776   



############## Metal####################
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Metal))
colnames(rep_data) <- c("Metal","rep")
data <- merge(data, rep_data, by = c("Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Metal - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                  estimate      se     zval    pval    ci.lb    ci.ub     
# MetalArsenic      0.0644  0.1568   0.4103  0.6816  -0.2431   0.3718   
# MetalCadmium     -0.0346  0.0556  -0.6225  0.5336  -0.1436   0.0744   
# MetalCobalt      -0.0597  0.1172  -0.5094  0.6105  -0.2895   0.1701   
# MetalCopper      -0.1458  0.0456  -3.1946  0.0014  -0.2353  -0.0564  ** 
# MetalIron         0.1005  0.1267   0.7928  0.4279  -0.1479   0.3488   
# MetalLead        -0.3733  0.1139  -3.2778  0.0010  -0.5966  -0.1501  ** 
# MetalManganese   -0.1129  0.1334  -0.8405  0.4006  -0.3763   0.1504   
# MetalMercury     -0.0450  0.0757  -0.5939  0.5526  -0.1933   0.1034   
# MetalNickel       -0.0004  0.1795   -0.0022  0.9983  -0.3521   0.3514   
# MetalSilver      -0.1264  0.1128  -1.1200  0.2627  -0.3475   0.0948   
# MetalZinc         0.0269  0.0850   0.3159  0.7520  -0.1398   0.1935   

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Arsenic "=1,"Cadmium"=1,"Cobalt"=1,
                                     "Copper"=1,"Iron"=1,"Lead "=1,"Manganese"=1,
                                     "Mercury"=1,"Nickel"=1,"Silver"=1,"Zinc"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                            Estimate Std. Error z value Pr(>|z|)   
# Cadmium - Arsenic  == 0   -0.098967   0.164363  -0.602  0.54709   
# Cobalt - Arsenic  == 0    -0.124081   0.195621  -0.634  0.52589   
# Copper - Arsenic  == 0    -0.210189   0.162995  -1.290  0.19721   
# Iron - Arsenic  == 0       0.036098   0.201628   0.179  0.85791   
# Lead  - Arsenic  == 0     -0.437684   0.192598  -2.273  0.02305  * 
# Manganese - Arsenic  == 0 -0.177306   0.206305  -0.859  0.39010   
# Mercury - Arsenic  == 0   -0.109319   0.173602  -0.630  0.52888   
# Nickel - Arsenic  == 0    -0.064746   0.238327  -0.272  0.78588   
# Silver - Arsenic  == 0    -0.190717   0.192586  -0.990  0.32203   
# Zinc - Arsenic  == 0      -0.037497   0.177606  -0.211  0.83279   
# Cobalt - Cadmium == 0     -0.025115   0.127814  -0.196  0.84422   
# Copper - Cadmium == 0     -0.111222   0.069871  -1.592  0.11142   
# Iron - Cadmium == 0        0.135065   0.138349   0.976  0.32894   
# Lead  - Cadmium == 0      -0.338718   0.125482  -2.699  0.00695  ** 
# Manganese - Cadmium == 0  -0.078339   0.143472  -0.546  0.58505   
# Mercury - Cadmium == 0     -0.010353   0.093444   -0.111  0.91178   
# Nickel - Cadmium == 0      0.034220   0.187848   0.182  0.85545   
# Silver - Cadmium == 0     -0.091750   0.117815  -0.779  0.43612   
# Zinc - Cadmium == 0        0.061469   0.099842   0.616  0.53812   
# Copper - Cobalt == 0      -0.086107   0.125409  -0.687  0.49233   
# Iron - Cobalt == 0         0.160179   0.172107   0.931  0.35201   
# Lead  - Cobalt == 0       -0.313603   0.163266  -1.921  0.05476  . 
# Manganese - Cobalt == 0   -0.053225   0.141022  -0.377  0.70586   
# Mercury - Cobalt == 0      0.014762   0.138565   0.107  0.91516   
# Nickel - Cobalt == 0       0.059335   0.214256   0.277  0.78183   
# Silver - Cobalt == 0      -0.066635   0.157275  -0.424  0.67179   
# Zinc - Cobalt == 0         0.086584   0.140335   0.617  0.53725   
# Iron - Copper == 0         0.246287   0.133661   1.829  0.06741  . 
# Lead  - Copper == 0       -0.227496   0.120850  -1.882  0.05977  . 
# Manganese - Copper == 0    0.032883   0.141461   0.232  0.81619   
# Mercury - Copper == 0      0.100869   0.088199   1.144  0.25276   
# Nickel - Copper == 0       0.145442   0.185120   0.786  0.43206   
# Silver - Copper == 0       0.019472   0.120056   0.162  0.87116   
# Zinc - Copper == 0         0.172691   0.093277   1.851  0.06411  . 
# Lead  - Iron == 0         -0.473783   0.170364  -2.781  0.00542  ** 
# Manganese - Iron == 0     -0.213304   0.184455  -1.157  0.24730   
# Mercury - Iron == 0       -0.145417   0.147592  -0.985  0.32449   
# Nickel - Iron == 0        -0.100844   0.219682  -0.459  0.64620   
# Silver - Iron == 0        -0.226815   0.169597  -1.337  0.18110   
# Zinc - Iron == 0          -0.073595   0.152139  -0.484  0.62857   
# Manganese - Lead  == 0     0.260378   0.175928   1.480  0.13887   
# Mercury - Lead  == 0       0.328365   0.136325   2.409  0.01601  * 
# Nickel - Lead  == 0        0.372938   0.212523   1.755  0.07929  . 
# Silver - Lead  == 0        0.246968   0.159823   1.545  0.12228   
# Zinc - Lead  == 0          0.400187   0.140751   2.843  0.00447  ** 
# Mercury - Manganese == 0   0.067987   0.152846   0.445  0.65646   
# Nickel - Manganese == 0    0.112560   0.224054   0.502  0.61540   
# Silver - Manganese == 0   -0.013311   0.169328  -0.079  0.93687   
# Zinc - Manganese == 0      0.139809   0.153707   0.910  0.36304   
# Nickel - Mercury == 0      0.044573   0.194759   0.229  0.81898   
# Silver - Mercury == 0     -0.081397   0.133570  -0.605  0.54526   
# Zinc - Mercury == 0        0.071822   0.112720   0.637  0.52401   
# Silver - Nickel == 0      -0.125970   0.211872  -0.595  0.55214   
# Zinc - Nickel == 0         0.027249   0.197136   0.138  0.89006   
# Zinc - Silver == 0         0.153219   0.136809   1.120  0.26273   





############## Trophic.level #################### 
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore          -0.0915  0.0361  -2.5351  0.0112  -0.1622  -0.0208  * 
# Trophic.levelPredator           -0.1878  0.0517  -3.6348  0.0003  -0.2891  -0.0866  *** 
# Trophic.levelPrimary producer    0.0881  0.0572   1.5385  0.1239  -0.0241   0.2003   
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                   Estimate Std. Error z value Pr(>|z|)    
# Predator - Herbivore == 0         -0.09634    0.06235  -1.545 0.122298   
# Primary producer - Herbivore == 0  0.17957    0.06611   2.716 0.006606  ** 
# Primary producer - Predator == 0   0.27591    0.07550   3.654 0.000258  *** 





############## Stressor:Taxa #################### 
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Taxa))
colnames(rep_data) <- c("Stressor","Taxa","rep")
data <- merge(data, rep_data, by = c("Stressor","Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TaxaAlgae+Protist       0.1508  0.0603   2.5012  0.0124   0.0326   0.2691  * 
# StressorA-W:TaxaAlgae+Protist    -0.1403  0.2937  -0.4778  0.6328  -0.7159   0.4353   
# StressorW:TaxaAlgae+Protist       0.2541  0.1413   1.7980  0.0722  -0.0229   0.5311  . 
# StressorA:TaxaDeuterostomia      -0.1211  0.0989  -1.2247  0.2207  -0.3148   0.0727   
# StressorA-W:TaxaDeuterostomia    -0.2209  0.1928  -1.1454  0.2520  -0.5988   0.1571   
# StressorW:TaxaDeuterostomia      -0.3446  0.1142  -3.0171  0.0026  -0.5685  -0.1208  ** 
# StressorA:TaxaEcdysozoa           0.0493  0.0874   0.5641  0.5727  -0.1220   0.2206   
# StressorH:TaxaEcdysozoa          -0.4562  0.1802  -2.5313  0.0114  -0.8095  -0.1030  * 
# StressorW:TaxaEcdysozoa          -0.3463  0.0845  -4.0980  <.0001  -0.5119  -0.1807  *** 
# StressorA:TaxaLophotrochozoa     -0.0608  0.0462  -1.3166  0.1880  -0.1513   0.0297   
# StressorA-W:TaxaLophotrochozoa   -0.2161  0.1674  -1.2908  0.1968  -0.5442   0.1120   
# StressorH:TaxaLophotrochozoa     -0.2193  0.1598  -1.3718  0.1701  -0.5326  0.0940   
# StressorW:TaxaLophotrochozoa     -0.0341  0.0602  -0.5669  0.5708  -0.1521   0.0839   
# StressorA:TaxaRadiata            -0.1477  0.1218  -1.2121  0.2255  -0.3865   0.0911   
# StressorW:TaxaRadiata            -0.3281  0.1050  -3.1236  0.0018  -0.5340  -0.1222  ** 



############## Stressor:Habitat #################### 
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Habitat))
colnames(rep_data) <- c("Stressor","Habitat","rep")
data <- merge(data, rep_data, by = c("Stressor","Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                  estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:HabitatPolar            0.0980   0.2619   0.3741   0.7083  -0.4154   0.6113   
# StressorW:HabitatPolar            0.4317   0.2085   2.0705   0.0384   0.0230   0.8403  * 
# StressorA:HabitatTemperate       -0.0402   0.0423  -0.9512   0.3415  -0.1231   0.0427   
# StressorA-W:HabitatTemperate     -0.1118   0.1375  -0.8136   0.4159  -0.3813   0.1576   
# StressorH:HabitatTemperate       -0.3643   0.1207  -3.0185   0.0025  -0.6009  -0.1278  ** 
# StressorW:HabitatTemperate       -0.0530   0.0594  -0.8916   0.3726  -0.1695   0.0635   
# StressorA:HabitatTropical         0.0400   0.0506   0.7906   0.4292  -0.0592   0.1392   
# StressorA-W:HabitatTropical      -0.4998   0.1827  -2.7361   0.0062  -0.8578  -0.1418  ** 
# StressorW:HabitatTropical        -0.3024   0.0586  -5.1573   <.0001  -0.4173  -0.1874  *** 


############## Stressor:Habitat:Taxa #################### 
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Habitat,data$Taxa))
colnames(rep_data) <- c("Stressor","Habitat","Taxa","rep")
data <- merge(data, rep_data, by = c("Stressor","Habitat","Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Habitat:Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                      estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:HabitatPolar:TaxaAlgae+Protist             0.0974   0.2585   0.3766   0.7065  -0.4094   0.6041   
# StressorW:HabitatPolar:TaxaAlgae+Protist             0.5856   0.2997   1.9542   0.0507  -0.0017   1.1729  . 
# StressorA:HabitatTemperate:TaxaAlgae+Protist         0.1116   0.0849   1.3149   0.1885  -0.0548   0.2780   
# StressorW:HabitatTemperate:TaxaAlgae+Protist         0.3146   0.2139   1.4711   0.1413  -0.1046   0.7338   
# StressorA:HabitatTropical:TaxaAlgae+Protist          0.1907   0.0920   2.0718   0.0383   0.0103   0.3711  * 
# StressorA-W:HabitatTropical:TaxaAlgae+Protist       -0.1919   0.2961  -0.6479   0.5170  -0.7723   0.3885   
# StressorW:HabitatTropical:TaxaAlgae+Protist         -0.0281   0.2443  -0.1150   0.9085  -0.5069   0.4507   
# StressorA:HabitatTemperate:TaxaDeuterostomia        -0.1433   0.0999  -1.4347   0.1514  -0.3392   0.0525   
# StressorA-W:HabitatTemperate:TaxaDeuterostomia      -0.2436   0.1927  -1.2642   0.2062  -0.6214   0.1331   
# StressorW:HabitatTemperate:TaxaDeuterostomia        -0.3721   0.1431  -2.5994   0.0093  -0.6526  -0.0915  ** 
# StressorW:HabitatTropical:TaxaDeuterostomia         -0.3075   0.1890  -1.6271   0.1037  -0.6779   0.0629   
# StressorA:HabitatTemperate:TaxaEcdysozoa            -0.1892   0.2328  -0.8127   0.4164  -0.6455   0.2671   
# StressorH:HabitatTemperate:TaxaEcdysozoa            -0.4621   0.1804  -2.5620   0.0104  -0.8156  -0.1086  * 
# StressorW:HabitatTemperate:TaxaEcdysozoa            -0.2408   0.1871  -1.2874   0.1980  -0.6075   0.1258   
# StressorA:HabitatTropical:TaxaEcdysozoa              0.0865   0.0955   0.9053   0.3653  -0.1008   0.2737   
# StressorW:HabitatTropical:TaxaEcdysozoa             -0.3703   0.0957  -3.8672   0.0001  -0.5579  -0.1826  *** 
# StressorW:HabitatPolar:TaxaLophotrochozoa            0.2898   0.2840   1.0204   0.3075  -0.2669   0.8466   
# StressorA:HabitatTemperate:TaxaLophotrochozoa       -0.0556   0.0529  -1.0512   0.2932  -0.1593   0.0481   
# StressorA-W:HabitatTemperate:TaxaLophotrochozoa     -0.0887   0.1942  -0.4568   0.6478  -0.4693   0.2919   
# StressorH:HabitatTemperate:TaxaLophotrochozoa       -0.3527   0.1791  -1.9687   0.0490  -0.7038  -0.0016  * 
# StressorW:HabitatTemperate:TaxaLophotrochozoa       -0.0056   0.0695  -0.0806   0.9358  -0.1418   0.1306   
# StressorA:HabitatTropical:TaxaLophotrochozoa        -0.0864   0.0992  -0.8712   0.3837  -0.2809   0.1080   
# StressorA-W:HabitatTropical:TaxaLophotrochozoa      -0.5646   0.3154  -1.7899   0.0735  -1.1828   0.0537  . 
# StressorW:HabitatTropical:TaxaLophotrochozoa        -0.1964   0.1370  -1.4335   0.1517  -0.4650   0.0721   
# StressorA:HabitatTropical:TaxaRadiata               -0.1479   0.1219  -1.2135   0.2249  -0.3867   0.0910   
# StressorW:HabitatTropical:TaxaRadiata               -0.3264   0.1052  -3.1029   0.0019  -0.5326  -0.1202  ** 






############## Stressor:Trophic.level #################### 
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Trophic.level))
colnames(rep_data) <- c("Stressor","Trophic.level","rep")
data <- merge(data, rep_data, by = c("Stressor","Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trophic.level- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                           estimate      se     zval    pval    ci.lb    ci.ub     
# StressorA:Trophic.levelHerbivore          -0.0211  0.0479  -0.4411  0.6591  -0.1150   0.0728   
# StressorA-W:Trophic.levelHerbivore        -0.2654  0.1868  -1.4208  0.1554  -0.6316   0.1007   
# StressorH:Trophic.levelHerbivore          -0.1808  0.1570  -1.1517  0.2494  -0.4886   0.1269   
# StressorW:Trophic.levelHerbivore          -0.1686  0.0555  -3.0366  0.0024  -0.2774  -0.0598  ** 
# StressorA:Trophic.levelPredator           -0.1061  0.0665  -1.5960  0.1105  -0.2364   0.0242   
# StressorA-W:Trophic.levelPredator         -0.3122  0.1422  -2.1946  0.0282  -0.5909  -0.0334  * 
# StressorH:Trophic.levelPredator           -0.4584  0.1831  -2.5033  0.0123  -0.8173  -0.0995  * 
# StressorW:Trophic.levelPredator           -0.2261  0.0744  -3.0396  0.0024  -0.3719  -0.0803  ** 
# StressorA:Trophic.levelPrimary producer    0.0921  0.0633   1.4546  0.1458  -0.0320   0.2162   
# StressorW:Trophic.levelPrimary producer    0.1058  0.1287   0.8216  0.4113  -0.1465   0.3580   

############## Habitat:Trophic.level #################### 
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Habitat,data$Trophic.level))
colnames(rep_data) <- c("Habitat","Trophic.level","rep")
data <- merge(data, rep_data, by = c("Habitat","Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Habitat:Trophic.level- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                  estimate      se     zval    pval 
# HabitatFrigid:Trophic.levelHerbivore              0.2959  0.3054   0.9688  0.3326 
# HabitatTemperate:Trophic.levelHerbivore          -0.0946  0.0465  -2.0364  0.0417 
# HabitatTropical:Trophic.levelHerbivore           -0.0979  0.0574  -1.7043  0.0883 
# HabitatTemperate:Trophic.levelPredator           -0.1250  0.0689  -1.8133  0.0698 
# HabitatTropical:Trophic.levelPredator            -0.2690  0.0772  -3.4820  0.0005 
# HabitatFrigid:Trophic.levelPrimary producer       0.3041  0.2085   1.4588  0.1446 
# HabitatTemperate:Trophic.levelPrimary producer    0.0932  0.0860   1.0836  0.2785 
# HabitatTropical:Trophic.levelPrimary producer     0.0444  0.0820   0.5421  0.5877 
#                                                  ci.lb    ci.ub      
# HabitatFrigid:Trophic.levelHerbivore            -0.3027   0.8946   
# HabitatTemperate:Trophic.levelHerbivore         -0.1857   -0.0035  * 
# HabitatTropical:Trophic.levelHerbivore          -0.2104  0.0147  . 
# HabitatTemperate:Trophic.levelPredator          -0.2600   0.0101  . 
# HabitatTropical:Trophic.levelPredator           -0.4204  -0.1176  *** 
# HabitatFrigid:Trophic.levelPrimary producer     -0.1045   0.7127   
# HabitatTemperate:Trophic.levelPrimary producer  -0.0754   0.2619   
# HabitatTropical:Trophic.levelPrimary producer   -0.1162   0.2051   



############## Trophic.level:Metal #################### 
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level,data$Metal))
colnames(rep_data) <- c("Trophic.level","Metal","rep")
data <- merge(data, rep_data, by = c("Trophic.level","Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level:Metal- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                               estimate      se     zval    pval 
# Trophic.levelHerbivore:MetalArsenic           0.1739  0.1786   0.9738  0.3302 
# Trophic.levelPrimary producer:MetalArsenic    0.1732  0.1914   0.9053  0.3653 
# Trophic.levelHerbivore:MetalCadmium          -0.1239  0.0637  -1.9448  0.0518 
# Trophic.levelPredator:MetalCadmium           -0.0391  0.1246  -0.3142  0.7534 
# Trophic.levelPrimary producer:MetalCadmium    0.4028  0.1421   2.8355  0.0046 
# Trophic.levelHerbivore:MetalCobalt           -0.0440  0.1825  -0.2410  0.8096 
# Trophic.levelPredator:MetalCobalt            -0.1330  0.1541  -0.8699  0.3844 
# Trophic.levelPrimary producer:MetalCobalt     0.1250  0.2907   0.4301  0.6671 
# Trophic.levelHerbivore:MetalCopper           -0.1756  0.0580  -3.0279  0.0025 
# Trophic.levelPredator:MetalCopper            -0.2333  0.0825  -2.8268  0.0047 
# Trophic.levelPrimary producer:MetalCopper    -0.0259  0.0865  -0.2997  0.7644 
# Trophic.levelPrimary producer:MetalIron       0.0992  0.1186   0.8363  0.4030 
# Trophic.levelHerbivore:MetalLead             -0.1331  0.1155  -1.1617  0.2454 
# Trophic.levelPredator:MetalLead              -2.8953  0.3896  -7.4312  <.0001 
# Trophic.levelHerbivore:MetalManganese        -0.0045  0.2624  -0.0172  0.9862 
# Trophic.levelPredator:MetalManganese         -0.1901  0.1556  -1.2217  0.2218 
# Trophic.levelHerbivore:MetalMercury           0.0339  0.0889   0.3809  0.7033 
# Trophic.levelPredator:MetalMercury           -0.1687  0.1148  -1.4693  0.1418 
# Trophic.levelHerbivore:MetalNickel            0.0106  0.1786   0.0592  0.9528 
# Trophic.levelHerbivore:MetalSilver           -0.1615  0.1335  -1.2005  0.2299 
# Trophic.levelPredator:MetalSilver            -0.1561  0.1881  -0.8298  0.4067 
# Trophic.levelHerbivore:MetalZinc              0.0100  0.1326   0.0755  0.9398 
# Trophic.levelPredator:MetalZinc              -0.0383  0.1312  -0.2921  0.7702 
# Trophic.levelPrimary producer:MetalZinc       0.1681  0.1666   1.0090  0.3130 
#                                              ci.lb    ci.ub      
# Trophic.levelHerbivore:MetalArsenic         -0.1761   0.5239   
# Trophic.levelPrimary producer:MetalArsenic  -0.2018   0.5483   
# Trophic.levelHerbivore:MetalCadmium         -0.2487  0.0010  . 
# Trophic.levelPredator:MetalCadmium          -0.2833   0.2050   
# Trophic.levelPrimary producer:MetalCadmium   0.1244   0.6812  ** 
# Trophic.levelHerbivore:MetalCobalt          -0.4018   0.3138   
# Trophic.levelPredator:MetalCobalt           -0.4360   0.1679   
# Trophic.levelPrimary producer:MetalCobalt   -0.4447   0.6947   
# Trophic.levelHerbivore:MetalCopper          -0.2893  -0.0619  ** 
# Trophic.levelPredator:MetalCopper           -0.3950  -0.0715  ** 
# Trophic.levelPrimary producer:MetalCopper   -0.1956   0.1437   
# Trophic.levelPrimary producer:MetalIron     -0.1333   0.3316   
# Trophic.levelHerbivore:MetalLead            -0.3604   0.0922   
# Trophic.levelPredator:MetalLead             -3.6589  -2.1317  *** 
# Trophic.levelHerbivore:MetalManganese       -0.5188   0.5098   
# Trophic.levelPredator:MetalManganese        -0.4950   0.1149   
# Trophic.levelHerbivore:MetalMercury         -0.1404   0.2082   
# Trophic.levelPredator:MetalMercury          -0.3937   0.0563   
# Trophic.levelHerbivore:MetalNickel          -0.3395   0.3606   
# Trophic.levelHerbivore:MetalSilver          -0.4252   0.1022   
# Trophic.levelPredator:MetalSilver           -0.5247   0.2126   
# Trophic.levelHerbivore:MetalZinc            -0.2499   0.2700   
# Trophic.levelPredator:MetalZinc             -0.2954   0.2188   
# Trophic.levelPrimary producer:MetalZinc     -0.1584   0.4946   




############## Stressor:Metal ####################
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Metal))
colnames(rep_data) <- c("Stressor","Metal","rep")
data <- merge(data, rep_data, by = c("Stressor","Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Metal - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                            estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:MetalArsenic      0.3543  0.2354   1.5049  0.1323  -0.1071   0.8157   
# StressorW:MetalArsenic     -0.1486  0.2014  -0.7380  0.4605  -0.5434   0.2461   
# StressorA:MetalCadmium      0.0357  0.0665   0.5364  0.5917  -0.0947   0.1661   
# StressorA-W:MetalCadmium   -0.1950  0.2376  -0.8207  0.4118  -0.6606   0.2707   
# StressorH:MetalCadmium     -0.3284  0.1465  -2.2422  0.0249  -0.6154  -0.0413  * 
# StressorW:MetalCadmium     -0.0294  0.1151  -0.2554  0.7984  -0.2549   0.1961   
# StressorA:MetalCobalt      -0.0843  0.1335  -0.6312  0.5279  -0.3460   0.1774   
# StressorW:MetalCobalt       0.2795  0.1881   1.4853  0.1375  -0.0893   0.6482   
# StressorA:MetalCopper      -0.1050  0.0568  -1.8472  0.0647  -0.2164   0.0064  . 
# StressorA-W:MetalCopper    -0.3256  0.2970  -1.0963  0.2729  -0.9076   0.2565   
# StressorW:MetalCopper      -0.2025  0.0708  -2.8593  0.0042  -0.3413  -0.0637  ** 
# StressorA:MetalIron         0.0131  0.1527   0.0858  0.9317  -0.2863   0.3124   
# StressorW:MetalIron         0.2633  0.2103   1.2518  0.2106  -0.1490   0.6756   
# StressorA:MetalLead        -0.1285  0.1726  -0.7445  0.4566  -0.4667   0.2098   
# StressorW:MetalLead        -0.5425  0.1463  -3.7094  0.0002  -0.8292  -0.2559  *** 
# StressorA:MetalManganese    0.1310  0.1706   0.7681  0.4424  -0.2033   0.4654   
# StressorW:MetalManganese    0.0923  0.2159   0.4273  0.6692  -0.3310   0.5155   
# StressorA:MetalMercury      0.1588  0.1011   1.5711  0.1162  -0.0393   0.3570   
# StressorA-W:MetalMercury   -0.1926  0.1708  -1.1279  0.2593  -0.5273   0.1421   
# StressorW:MetalMercury     -0.2049  0.0964  -2.1248  0.0336  -0.3939  -0.0159  * 
# StressorA:MetalNickel       0.1835  0.2449   0.7495  0.4536  -0.2964   0.6634   
# StressorW:MetalNickel      -0.1907  0.2482  -0.7683  0.4423  -0.6770   0.2957   
# StressorA:MetalSilver      -0.1062  0.1261  -0.8417  0.4000  -0.3534   0.1411   
# StressorW:MetalSilver      -0.0614  0.2389  -0.2571  0.7971  -0.5296   0.4068   
# StressorA:MetalZinc         0.0164  0.0925   0.1769  0.8596  -0.1649   0.1976   
# StressorW:MetalZinc        -0.0064  0.2431  -0.0262  0.9791  -0.4829   0.4701   



############## Stressor:Trait ####################
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Trait))
colnames(rep_data) <- c("Stressor","Trait","rep")
data <- merge(data, rep_data, by = c("Stressor","Trait"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trait - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TraitAccumulation     -0.0097  0.0468  -0.2065  0.8364  -0.1014   0.0821   
# StressorA-W:TraitAccumulation   -0.2250  0.1540  -1.4603  0.1442  -0.5269   0.0770   
# StressorH:TraitAccumulation     -0.1571  0.1575  -0.9976  0.3185  -0.4659   0.1516   
# StressorW:TraitAccumulation     -0.2611  0.0725  -3.6001  0.0003  -0.4033  -0.1190  *** 
# StressorA:TraitBehaviour         0.0092  0.1874   0.0491  0.9608  -0.3582   0.3766   
# StressorW:TraitBehaviour        -0.1621  0.1902  -0.8523  0.3940  -0.5349   0.2107   
# StressorA:TraitGrowth            0.0032  0.0607   0.0526  0.9580  -0.1157   0.1221   
# StressorW:TraitGrowth           -0.1102  0.0959  -1.1495  0.2503  -0.2981   0.0777   
# StressorA:TraitMetabolism       -0.0488  0.0527  -0.9274  0.3537  -0.1521   0.0544   
# StressorA-W:TraitMetabolism     -0.5465  0.2082  -2.6247  0.0087  -0.9546  -0.1384  ** 
# StressorH:TraitMetabolism       -0.4202  0.2010  -2.0902  0.0366  -0.8141  -0.0262  * 
# StressorW:TraitMetabolism       -0.0152  0.0644  -0.2370  0.8127  -0.1414   0.1109   
# StressorA:TraitReproduction      0.1369  0.1278   1.0718  0.2838  -0.1135   0.3874   
# StressorW:TraitReproduction     -0.4631  0.1363  -3.3985  0.0007  -0.7301  -0.1960  *** 
# StressorA:TraitSurvival          0.0482  0.1337   0.3579  0.7204  -0.2157   0.3121   
# StressorW:TraitSurvival         -0.2566  0.1064  -2.4111  0.0159  -0.4651  -0.0480  * 


############## Stressor:Trait:Habitat #################### 
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Trait,data$Habitat))
colnames(rep_data) <- c("Stressor","Trait","Habitat","rep")
data <- merge(data, rep_data, by = c("Stressor","Trait","Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trait:Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:TraitMetabolism:HabitatPolar              0.6482   0.2359   2.7480   0.0060   0.1859   1.1106  ** 
# StressorA:TraitAccumulation:HabitatTemperate       -0.0318   0.0552  -0.5766   0.5642  -0.1401   0.0764   
# StressorA-W:TraitAccumulation:HabitatTemperate     -0.1533   0.1651  -0.9287   0.3530  -0.4769   0.1703   
# StressorH:TraitAccumulation:HabitatTemperate       -0.3524   0.1632  -2.1592   0.0308  -0.6723  -0.0325  * 
# StressorW:TraitAccumulation:HabitatTemperate       -0.0795   0.0889  -0.8938   0.3714  -0.2537   0.0948   
# StressorA:TraitBehaviour:HabitatTemperate          -0.0714   0.2528  -0.2824   0.7777  -0.5669   0.4241   
# StressorA:TraitGrowth:HabitatTemperate             -0.0129   0.0874  -0.1478   0.8825  -0.1843   0.1584   
# StressorW:TraitGrowth:HabitatTemperate             -0.3342   0.1448  -2.3087   0.0210  -0.6179  -0.0505  * 
# StressorA:TraitMetabolism:HabitatTemperate         -0.0614   0.0669  -0.9172   0.3590  -0.1926   0.0698   
# StressorH:TraitMetabolism:HabitatTemperate         -0.2922   0.2163  -1.3508   0.1767  -0.7163   0.1318   
# StressorW:TraitMetabolism:HabitatTemperate          0.1650   0.0970   1.7009   0.0890  -0.0251   0.3550  . 
# StressorA:TraitSurvival:HabitatTemperate           -0.1740   0.2603  -0.6682   0.5040  -0.6842   0.3363   
# StressorW:TraitSurvival:HabitatTemperate           -0.2284   0.1934  -1.1813   0.2375  -0.6074   0.1506   
# StressorA:TraitAccumulation:HabitatTropical         0.0700   0.0826   0.8471   0.3969  -0.0919   0.2319   
# StressorW:TraitAccumulation:HabitatTropical        -0.5944   0.1258  -4.7236   <.0001  -0.8410  -0.3478  *** 
# StressorA:TraitBehaviour:HabitatTropical            0.1141   0.2690   0.4242   0.6714  -0.4131   0.6413   
# StressorW:TraitBehaviour:HabitatTropical           -0.2355   0.2209  -1.0662   0.2863  -0.6684   0.1974   
# StressorA:TraitGrowth:HabitatTropical               0.0195   0.0832   0.2343   0.8148  -0.1436   0.1826   
# StressorW:TraitGrowth:HabitatTropical               0.0011   0.1231   0.0086   0.9932  -0.2402   0.2423   
# StressorA:TraitMetabolism:HabitatTropical          -0.0374   0.0825  -0.4528   0.6507  -0.1990   0.1243   
# StressorA-W:TraitMetabolism:HabitatTropical        -0.6151   0.2319  -2.6524   0.0080  -1.0696  -0.1606  ** 
# StressorW:TraitMetabolism:HabitatTropical          -0.2437   0.0883  -2.7607   0.0058  -0.4167  -0.0707  ** 
# StressorA:TraitReproduction:HabitatTropical         0.2069   0.1325   1.5613   0.1185  -0.0528   0.4666   
# StressorW:TraitReproduction:HabitatTropical        -0.4957   0.1350  -3.6724   0.0002  -0.7603  -0.2312  *** 
# StressorA:TraitSurvival:HabitatTropical             0.1381   0.1531   0.9022   0.3670  -0.1619   0.4381   
# StressorW:TraitSurvival:HabitatTropical            -0.3063   0.1242  -2.4665   0.0136  -0.5498  -0.0629  * 




############## Trophic.level-Accumulation #################### 
data <- subset(data, Trait == "Accumulation")
nrow(data)
data<-subset(data,data$Trait!="Abundance")
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub 
# Trophic.levelHerbivore          -0.0637  0.0625  -1.0185  0.3084  -0.1862   0.0588   
# Trophic.levelPredator           -0.2017  0.0837  -2.4103  0.0159  -0.3657  -0.0377  * 
# Trophic.levelPrimary producer    0.3452  0.1517   2.2754  0.0229   0.0479   0.6426  * 
# 
# Trophic.levelHerbivore   
# Trophic.levelPredator  * 
# Trophic.levelPrimary producer  * 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                   Estimate Std. Error z value Pr(>|z|)   
# Predator - Herbivore == 0          -0.1380     0.1044  -1.321  0.18634   
# Primary producer - Herbivore == 0   0.4089     0.1559   2.623  0.00870  ** 
# Primary producer - Predator == 0    0.5469     0.1733   3.157  0.00160  ** 








### C. Removing extreme value ###

############## overall #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb    ci.ub     
# -0.0747  0.0276  -2.7036  0.0069  -0.1288  -0.0205  ** 


############## Stressor #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
#               estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA     -0.0061  0.0341  -0.1793  0.8577  -0.0729   0.0607   
# StressorA-W   -0.2457  0.1158  -2.1214  0.0339  -0.4728  -0.0187  * 
# StressorH     -0.2948  0.1225  -2.4059  0.0161  -0.5350  -0.0546  * 
# StressorW     -0.1428  0.0427  -3.3427  0.0008  -0.2265  -0.0591  *** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("A"=1,"A-W"=1,"H"=1,"W"=1), 
                                   type="Tukey"))), test=adjusted("none"))

#             Estimate Std. Error z value Pr(>|z|)   
# A-W - A == 0 -0.23963    0.11731  -2.043  0.04108  * 
# H - A == 0   -0.28872    0.12720  -2.270  0.02322  * 
# W - A == 0   -0.13668    0.05308  -2.575  0.01002  * 
# H - A-W == 0 -0.04909    0.16863  -0.291  0.77097   
# W - A-W == 0  0.10295    0.11941   0.862  0.38862   
# W - H == 0    0.15203    0.12978   1.171  0.24140   


############## Taxa ####################
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Taxa))
colnames(rep_data) <- c("Taxa","rep")
data <- merge(data, rep_data, by = c("Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                      estimate      se     zval    pval    ci.lb    ci.ub     
# TaxaAlgae+Protist     0.1653  0.0572   2.8870  0.0039   0.0531   0.2775  ** 
# TaxaDeuterostomia    -0.2044  0.0804  -2.5425  0.0110  -0.3619  -0.0468  * 
# TaxaEcdysozoa        -0.1931  0.0595  -3.2468  0.0012  -0.3096  -0.0765  ** 
# TaxaLophotrochozoa   -0.0614  0.0388  -1.5847  0.1130  -0.1374   0.0145   
# TaxaRadiata          -0.2136  0.0799  -2.6719  0.0075  -0.3702  -0.0569  ** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Algae+Protist"=1,"Deuterostomia"=1,"Ecdysozoa"=1,"Lophotrochozoa"=1,"Radiata"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                     Estimate Std. Error z value Pr(>|z|)    
# Deuterostomia - Algae+Protist == 0  -0.36966    0.09868  -3.746 0.000180  *** 
# Ecdysozoa - Algae+Protist == 0      -0.35835    0.08252  -4.343 1.41e-05  *** 
# Lophotrochozoa - Algae+Protist == 0 -0.22669    0.06772  -3.347 0.000816  *** 
# Radiata - Algae+Protist == 0        -0.37884    0.09832  -3.853 1.17e-04  *** 
# Ecdysozoa - Deuterostomia == 0       0.01131    0.09982   0.113 0.909809   
# Lophotrochozoa - Deuterostomia == 0  0.14296    0.08916   1.603 0.108847   
# Radiata - Deuterostomia == 0        -0.00918    0.11336  -0.081 0.935442   
# Lophotrochozoa - Ecdysozoa == 0      0.13165    0.07046   1.868 0.061698  . 
# Radiata - Ecdysozoa == 0            -0.02049    0.09962  -0.206 0.837051   
# Radiata - Lophotrochozoa == 0       -0.15214    0.08883  -1.713 0.086760  . 


############## Habitat ####################
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Habitat))
colnames(rep_data) <- c("Habitat","rep")
data <- merge(data, rep_data, by = c("Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                    estimate      se     zval    pval    ci.lb    ci.ub     
# HabitatPolar        0.5003  0.1625   3.0788  0.0021   0.1818   0.8188  ** 
# HabitatTemperate   -0.0700  0.0364  -1.9258  0.0541  -0.1413   0.0012  . 
# HabitatTropical    -0.1123  0.0398  -2.8175  0.0048  -0.1903  -0.0342  ** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Polar"=1,"Temperate"=1,"Tropical"=1), 
                                   type="Tukey"))), test=adjusted("none"))
# Estimate Std. Error z value Pr(>|z|)   
# Temperate - Polar == 0    -0.57031    0.16651  -3.425  0.00061  *** 
# Tropical - Polar == 0     -0.61254    0.16731  -3.661  0.00025  *** 
# Tropical - Temperate == 0 -0.04223    0.05391  -0.783  0.43339   



############## Trait ####################
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Trait))
colnames(rep_data) <- c("Trait","rep")
data <- merge(data, rep_data, by = c("Trait"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trait - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                     estimate      se     zval    pval    ci.lb    ci.ub    
# TraitAbundance       0.0879  0.1129   0.7787  0.4362  -0.1334   0.3092   
# TraitAccumulation   -0.0860  0.0415  -2.0719  0.0383  -0.1674  -0.0046  * 
# TraitBehaviour      -0.1384  0.1330  -1.0404  0.2982  -0.3991   0.1223   
# TraitGrowth         -0.0462  0.0543  -0.8513  0.3946  -0.1526   0.0602   
# TraitMetabolism     -0.0637  0.0422  -1.5095  0.1312  -0.1465   0.0190   
# TraitReproduction   -0.1521  0.0996  -1.5278  0.1266  -0.3473   0.0430   
# TraitSurvival       -0.1685  0.0872  -1.9327  0.0533  -0.3395   0.0024  . 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Abundance"=1,"Accumulation"=1,"Behaviour"=1,
                                     "Growth"=1,"Metabolism"=1,"Reproduction"=1,"Survival"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                  Estimate Std. Error z value Pr(>|z|)  
# Accumulation - Abundance == 0    -0.17394    0.11923  -1.459   0.1446   
# Behaviour - Abundance == 0       -0.22633    0.17379  -1.302   0.1928   
# Growth - Abundance == 0          -0.13313    0.12147  -1.104   0.2695   
# Metabolism - Abundance == 0      -0.15166    0.11783  -1.287   0.1981   
# Reproduction - Abundance == 0    -0.24005    0.14936  -1.607   0.1080   
# Survival - Abundance == 0        -0.25648    0.14048  -1.826   0.0679  . 
# Behaviour - Accumulation == 0    -0.05240    0.13729  -0.382   0.7027   
# Growth - Accumulation == 0        0.03980    0.06430   0.619   0.5359   
# Metabolism - Accumulation == 0    0.02228    0.05460   0.408   0.6832   
# Reproduction - Accumulation == 0 -0.06611    0.10495  -0.630   0.5287   
# Survival - Accumulation == 0     -0.08254    0.09373  -0.881   0.3785   
# Growth - Behaviour == 0           0.09220    0.13996   0.659   0.5101   
# Metabolism - Behaviour == 0       0.07468    0.13836   0.540   0.5894   
# Reproduction - Behaviour == 0    -0.01372    0.16410  -0.084   0.9334   
# Survival - Behaviour == 0        -0.03014    0.15537  -0.194   0.8462   
# Metabolism - Growth == 0         -0.01752    0.06399  -0.274   0.7842   
# Reproduction - Growth == 0       -0.10592    0.10590  -1.000   0.3172   
# Survival - Growth == 0           -0.12234    0.09878  -1.239   0.2155   
# Reproduction - Metabolism == 0   -0.08839    0.10609  -0.833   0.4047   
# Survival - Metabolism == 0       -0.10482    0.09501  -1.103   0.2699   
# Survival - Reproduction == 0     -0.01643    0.12594  -0.130   0.8962   


############## Metal####################
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Metal))
colnames(rep_data) <- c("Metal","rep")
data <- merge(data, rep_data, by = c("Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Metal - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                  estimate      se     zval    pval    ci.lb    ci.ub      
# MetalArsenic      0.1036  0.1511   0.6857  0.4929  -0.1925   0.3997   
# MetalCadmium     -0.0394  0.0545  -0.7229  0.4698  -0.1463   0.0675   
# MetalCobalt      -0.0431  0.1183  -0.3647  0.7154  -0.2750   0.1887   
# MetalCopper      -0.1577  0.0445  -3.5455  0.0004  -0.2448  -0.0705  *** 
# MetalIron         0.2105  0.1167   1.8031  0.0714  -0.0183   0.4393  . 
# MetalLead        -0.2938  0.1066  -2.7568  0.0058  -0.5028  -0.0849  ** 
# MetalManganese   -0.1009  0.1380  -0.7312  0.4646  -0.3715   0.1696   
# MetalMercury     -0.0381  0.0731  -0.5213  0.6022  -0.1815   0.1052   
# MetalNickel      -0.0100  0.1707  -0.0588  0.9531  -0.3447   0.3246   
# MetalSilver      -0.1223  0.1167  -1.0481  0.2946  -0.3509   0.1064   
# MetalZinc         0.0336  0.0867   0.3875  0.6984  -0.1364   0.2036   

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Arsenic "=1,"Cadmium"=1,"Cobalt"=1,
                                     "Copper"=1,"Iron"=1,"Lead "=1,"Manganese"=1,
                                     "Mercury"=1,"Nickel"=1,"Silver"=1,"Zinc"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                           Estimate Std. Error z value Pr(>|z|)   
# Cadmium - Arsenic  == 0   -0.143006   0.159284  -0.898  0.36929   
# Cobalt - Arsenic  == 0    -0.146727   0.191783  -0.765  0.44423   
# Copper - Arsenic  == 0    -0.261257   0.157271  -1.661  0.09668  . 
# Iron - Arsenic  == 0       0.106889   0.190917   0.560  0.57557   
# Lead  - Arsenic  == 0     -0.397433   0.184137  -2.158  0.03090  * 
# Manganese - Arsenic  == 0 -0.204528   0.204535  -1.000  0.31733   
# Mercury - Arsenic  == 0   -0.141713   0.167482  -0.846  0.39748   
# Nickel - Arsenic  == 0    -0.113622   0.227975  -0.498  0.61820   
# Silver - Arsenic  == 0    -0.225852   0.190530  -1.185  0.23586   
# Zinc - Arsenic  == 0      -0.069986   0.173652  -0.403  0.68693   
# Cobalt - Cadmium == 0      -0.003721   0.128983   -0.029  0.97698   
# Copper - Cadmium == 0     -0.118251   0.068977  -1.714  0.08646  . 
# Iron - Cadmium == 0        0.249894   0.128831   1.940  0.05242  . 
# Lead  - Cadmium == 0      -0.254427   0.118925  -2.139  0.03240  * 
# Manganese - Cadmium == 0  -0.061523   0.147121  -0.418  0.67582   
# Mercury - Cadmium == 0     0.001292   0.090915   0.014  0.98866   
# Nickel - Cadmium == 0      0.029384   0.179223   0.164  0.86977   
# Silver - Cadmium == 0     -0.082847   0.122380  -0.677  0.49843   
# Zinc - Cadmium == 0        0.073019   0.101303   0.721  0.47103   
# Copper - Cobalt == 0      -0.114530   0.126155  -0.908  0.36396   
# Iron - Cobalt == 0         0.253615   0.165859   1.529  0.12624   
# Lead  - Cobalt == 0       -0.250706   0.159144  -1.575  0.11518   
# Manganese - Cobalt == 0   -0.057802   0.149626  -0.386  0.69927   
# Mercury - Cobalt == 0      0.005013   0.138401   0.036  0.97110   
# Nickel - Cobalt == 0       0.033105   0.207663   0.159  0.87334   
# Silver - Cobalt == 0      -0.079126   0.162104  -0.488  0.62547   
# Zinc - Cobalt == 0         0.076740   0.143466   0.535  0.59272   
# Iron - Copper == 0         0.368145   0.124907   2.947  0.00320  ** 
# Lead  - Copper == 0       -0.136176   0.114055  -1.194  0.23250   
# Manganese - Copper == 0    0.056728   0.144777   0.392  0.69518   
# Mercury - Copper == 0      0.119543   0.085473   1.399  0.16193   
# Nickel - Copper == 0       0.147635   0.176408   0.837  0.40265   
# Silver - Copper == 0       0.035404   0.123772   0.286  0.77485   
# Zinc - Copper == 0         0.191270   0.095099   2.011  0.04430  * 
# Lead  - Iron == 0         -0.504321   0.158071  -3.190  0.00142  ** 
# Manganese - Iron == 0     -0.311417   0.180653  -1.724  0.08474  . 
# Mercury - Iron == 0       -0.248602   0.137746  -1.805  0.07111  . 
# Nickel - Iron == 0        -0.220510   0.206828  -1.066  0.28635   
# Silver - Iron == 0        -0.332741   0.164998  -2.017  0.04373  * 
# Zinc - Iron == 0          -0.176875   0.145133  -1.219  0.22296   
# Manganese - Lead  == 0     0.192904   0.174303   1.107  0.26841   
# Mercury - Lead  == 0       0.255719   0.128999   1.982  0.04744  * 
# Nickel - Lead  == 0        0.283811   0.201267   1.410  0.15850   
# Silver - Lead  == 0        0.171580   0.157755   1.088  0.27676   
# Zinc - Lead  == 0          0.327446   0.136522   2.398  0.01646  * 
# Mercury - Manganese == 0   0.062815   0.155229   0.405  0.68573   
# Nickel - Manganese == 0    0.090907   0.219497   0.414  0.67876   
# Silver - Manganese == 0   -0.021324   0.176015  -0.121  0.90357   
# Zinc - Manganese == 0      0.133542   0.159109   0.846  0.39778   
# Nickel - Mercury == 0      0.028092   0.185735   0.151  0.87978   
# Silver - Mercury == 0     -0.084139   0.136752  -0.615  0.53838   
# Zinc - Mercury == 0        0.071727   0.112654   0.637  0.52432   
# Silver - Nickel == 0      -0.112231   0.206732  -0.543  0.58721   
# Zinc - Nickel == 0         0.043635   0.190571   0.229  0.81889   
# Zinc - Silver == 0         0.155866   0.142028   1.097  0.27245   



############## Trophic.level #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore          -0.1001  0.0351  -2.8533  0.0043  -0.1689  -0.0313  ** 
# Trophic.levelPredator           -0.1632  0.0495  -3.2977  0.0010  -0.2601  -0.0662  *** 
# Trophic.levelPrimary producer    0.1160  0.0552   2.1013  0.0356   0.0078   0.2242  * 

summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                   Estimate Std. Error z value Pr(>|z|)    
# Predator - Herbivore == 0         -0.06307    0.06011  -1.049 0.294119   
# Primary producer - Herbivore == 0  0.21610    0.06427   3.363 0.000772  *** 
# Primary producer - Predator == 0   0.27917    0.07285   3.832 1.27e-04  *** 


############## Stressor:Taxa #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Taxa))
colnames(rep_data) <- c("Stressor","Taxa","rep")
data <- merge(data, rep_data, by = c("Stressor","Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                  estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TaxaAlgae+Protist       0.1567  0.0596   2.6281  0.0086   0.0398   0.2735  ** 
# StressorA-W:TaxaAlgae+Protist    -0.1101  0.2726  -0.4038  0.6863  -0.6443   0.4242   
# StressorW:TaxaAlgae+Protist       0.2450  0.1280   1.9135  0.0557  -0.0059   0.4959  . 
# StressorA:TaxaDeuterostomia      -0.1152  0.0999  -1.1536  0.2487  -0.3110   0.0805   
# StressorA-W:TaxaDeuterostomia    -0.1929  0.2014  -0.9577  0.3382  -0.5876   0.2018   
# StressorW:TaxaDeuterostomia      -0.3382  0.1170  -2.8907  0.0038  -0.5675  -0.1089  ** 
# StressorA:TaxaEcdysozoa           0.0555  0.0864   0.6419  0.5209  -0.1139   0.2249   
# StressorH:TaxaEcdysozoa          -0.4393  0.1837  -2.3919  0.0168  -0.7993  -0.0793  * 
# StressorW:TaxaEcdysozoa          -0.3305  0.0807  -4.0961  <.0001  -0.4886  -0.1724  *** 
# StressorA:TaxaLophotrochozoa     -0.0595  0.0459  -1.2968  0.1947  -0.1495   0.0304   
# StressorA-W:TaxaLophotrochozoa   -0.2184  0.1765  -1.2371  0.2160  -0.5644   0.1276   
# StressorH:TaxaLophotrochozoa     -0.2216  0.1648  -1.3448  0.1787  -0.5446  0.1014   
# StressorW:TaxaLophotrochozoa     -0.0347  0.0612  -0.5669  0.5708  -0.1546   0.0852   
# StressorA:TaxaRadiata            -0.1354  0.1168  -1.1588  0.2466  -0.3644   0.0936   
# StressorW:TaxaRadiata            -0.2738  0.1020  -2.6844  0.0073  -0.4738  -0.0739  ** 


############## Stressor:Habitat #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Habitat))
colnames(rep_data) <- c("Stressor","Habitat","rep")
data <- merge(data, rep_data, by = c("Stressor","Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:HabitatPolar            0.4836   0.2391   2.0224   0.0431   0.0149   0.9522  * 
# StressorW:HabitatPolar            0.5122   0.1954   2.6209   0.0088   0.1292   0.8953  ** 
# StressorA:HabitatTemperate       -0.0347   0.0411  -0.8442   0.3986  -0.1152   0.0458   
# StressorA-W:HabitatTemperate     -0.0985   0.1443  -0.6826   0.4949  -0.3812   0.1843   
# StressorH:HabitatTemperate       -0.3582   0.1223  -2.9279   0.0034  -0.5980  -0.1184  ** 
# StressorW:HabitatTemperate       -0.0533   0.0594  -0.8979   0.3692  -0.1697   0.0631   
# StressorA:HabitatTropical         0.0282   0.0490   0.5750   0.5653  -0.0679   0.1243   
# StressorA-W:HabitatTropical      -0.4722   0.1805  -2.6167   0.0089  -0.8259  -0.1185  ** 
# StressorW:HabitatTropical        -0.2775   0.0552  -5.0229   <.0001  -0.3858  -0.1692  *** 


############## Stressor:Habitat:Taxa #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Habitat,data$Taxa))
colnames(rep_data) <- c("Stressor","Habitat","Taxa","rep")
data <- merge(data, rep_data, by = c("Stressor","Habitat","Taxa"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Habitat:Taxa - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                    estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:HabitatPolar:TaxaAlgae+Protist             0.4835   0.2363   2.0463   0.0407   0.0204   0.9466  * 
# StressorW:HabitatPolar:TaxaAlgae+Protist             0.7029   0.2680   2.6226   0.0087   0.1776   1.2282  ** 
# StressorA:HabitatTemperate:TaxaAlgae+Protist         0.1209   0.0822   1.4711   0.1413  -0.0402   0.2819   
# StressorW:HabitatTemperate:TaxaAlgae+Protist         0.2862   0.2090   1.3697   0.1708  -0.1233   0.6958   
# StressorA:HabitatTropical:TaxaAlgae+Protist          0.1427   0.0910   1.5691   0.1166  -0.0356   0.3210   
# StressorA-W:HabitatTropical:TaxaAlgae+Protist       -0.1792   0.2723  -0.6580   0.5105  -0.7130   0.3546   
# StressorW:HabitatTropical:TaxaAlgae+Protist         -0.0428   0.1982  -0.2158   0.8292  -0.4311   0.3456   
# StressorA:HabitatTemperate:TaxaDeuterostomia        -0.1386   0.0993  -1.3959   0.1627  -0.3332   0.0560   
# StressorA-W:HabitatTemperate:TaxaDeuterostomia      -0.2047   0.1997  -1.0249   0.3054  -0.5961   0.1867   
# StressorW:HabitatTemperate:TaxaDeuterostomia        -0.3526   0.1460  -2.4162   0.0157  -0.6387  -0.0666  * 
# StressorW:HabitatTropical:TaxaDeuterostomia         -0.3193   0.1894  -1.6861   0.0918  -0.6905   0.0519  . 
# StressorA:HabitatTemperate:TaxaEcdysozoa            -0.1855   0.2347  -0.7907   0.4291  -0.6455   0.2744   
# StressorH:HabitatTemperate:TaxaEcdysozoa            -0.4408   0.1815  -2.4290   0.0151  -0.7964  -0.0851  * 
# StressorW:HabitatTemperate:TaxaEcdysozoa            -0.2416   0.1904  -1.2692   0.2044  -0.6148   0.1315   
# StressorA:HabitatTropical:TaxaEcdysozoa              0.0921   0.0920   1.0013   0.3167  -0.0882   0.2723   
# StressorW:HabitatTropical:TaxaEcdysozoa             -0.3477   0.0880  -3.9517   <.0001  -0.5202  -0.1753  *** 
# StressorW:HabitatPolar:TaxaLophotrochozoa            0.3030   0.2791   1.0853   0.2778  -0.2442   0.8501   
# StressorA:HabitatTemperate:TaxaLophotrochozoa       -0.0550   0.0518  -1.0624   0.2881  -0.1565   0.0465   
# StressorA-W:HabitatTemperate:TaxaLophotrochozoa     -0.0797   0.2058  -0.3873   0.6985  -0.4830   0.3236   
# StressorH:HabitatTemperate:TaxaLophotrochozoa       -0.3577   0.1834  -1.9507   0.0511  -0.7171   0.0017  . 
# StressorW:HabitatTemperate:TaxaLophotrochozoa       -0.0090   0.0698  -0.1290   0.8973  -0.1459   0.1279   
# StressorA:HabitatTropical:TaxaLophotrochozoa        -0.0796   0.0957  -0.8310   0.4060  -0.2672   0.1081   
# StressorA-W:HabitatTropical:TaxaLophotrochozoa      -0.5793   0.3268  -1.7724   0.0763  -1.2198   0.0613  . 
# StressorW:HabitatTropical:TaxaLophotrochozoa        -0.1952   0.1359  -1.4359   0.1510  -0.4616   0.0712   
# StressorA:HabitatTropical:TaxaRadiata               -0.1354   0.1154  -1.1736   0.2406  -0.3615   0.0907   
# StressorW:HabitatTropical:TaxaRadiata               -0.2742   0.1010  -2.7159   0.0066  -0.4721  -0.0763  ** 


############## Stressor:Trophic.level #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Trophic.level))
colnames(rep_data) <- c("Stressor","Trophic.level","rep")
data <- merge(data, rep_data, by = c("Stressor","Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trophic.level- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                           estimate      se     zval    pval    ci.lb    ci.ub     
# StressorA:Trophic.levelHerbivore          -0.0309  0.0466  -0.6629  0.5074  -0.1221   0.0604   
# StressorA-W:Trophic.levelHerbivore        -0.2605  0.1860  -1.4007  0.1613  -0.6251   0.1040   
# StressorH:Trophic.levelHerbivore          -0.1849  0.1600  -1.1558  0.2477  -0.4985   0.1286   
# StressorW:Trophic.levelHerbivore          -0.1763  0.0541  -3.2608  0.0011  -0.2823  -0.0703  ** 
# StressorA:Trophic.levelPredator           -0.0945  0.0653  -1.4461  0.1482  -0.2226   0.0336   
# StressorA-W:Trophic.levelPredator         -0.2847  0.1484  -1.9190  0.0550  -0.5755   0.0061  . 
# StressorH:Trophic.levelPredator           -0.4343  0.1827  -2.3777  0.0174  -0.7924  -0.0763  * 
# StressorW:Trophic.levelPredator           -0.1861  0.0718  -2.5900  0.0096  -0.3268  -0.0453  ** 
# StressorA:Trophic.levelPrimary producer    0.1226  0.0614   1.9980  0.0457   0.0023   0.2428  * 
# StressorW:Trophic.levelPrimary producer    0.1166  0.1210   0.9641  0.3350  -0.1205   0.3537   

############## Habitat:Trophic.level #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Habitat,data$Trophic.level))
colnames(rep_data) <- c("Habitat","Trophic.level","rep")
data <- merge(data, rep_data, by = c("Habitat","Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Habitat:Trophic.level- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                estimate      se     zval    pval    ci.lb    ci.ub 
# HabitatPolar:Trophic.levelHerbivore              0.3078  0.2944   1.0457  0.2957  -0.2691   0.8848   
# HabitatTemperate:Trophic.levelHerbivore          -0.0947  0.0451  -2.0993  0.0358  -0.1830   -0.0063  * 
# HabitatTropical:Trophic.levelHerbivore           -0.1202  0.0537  -2.2358  0.0254  -0.2255  -0.0148  * 
# HabitatTemperate:Trophic.levelPredator           -0.1126  0.0660  -1.7044  0.0883  -0.2420   0.0169  . 
# HabitatTropical:Trophic.levelPredator            -0.2245  0.0715  -3.1381  0.0017  -0.3647  -0.0843  ** 
# HabitatPolar:Trophic.levelPrimary producer       0.5783  0.1865   3.1018  0.0019   0.2129   0.9438  ** 
# HabitatTemperate:Trophic.levelPrimary producer    0.1058  0.0822   1.2881  0.1977  -0.0552   0.2669   
# HabitatTropical:Trophic.levelPrimary producer     0.0438  0.0784   0.5590  0.5762  -0.1098   0.1974   
# 
# HabitatPolar:Trophic.levelHerbivore   
# HabitatTemperate:Trophic.levelHerbivore  * 
# HabitatTropical:Trophic.levelHerbivore  * 
# HabitatTemperate:Trophic.levelPredator  . 
# HabitatTropical:Trophic.levelPredator  ** 
# HabitatPolar:Trophic.levelPrimary producer  ** 
# HabitatTemperate:Trophic.levelPrimary producer   
# HabitatTropical:Trophic.levelPrimary producer   



############## Trophic.level:Metal #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level,data$Metal))
colnames(rep_data) <- c("Trophic.level","Metal","rep")
data <- merge(data, rep_data, by = c("Trophic.level","Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level:Metal- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                               estimate      se     zval    pval    ci.lb 
# Trophic.levelHerbivore:MetalArsenic           0.1829  0.1796   1.0185  0.3084  -0.1691 
# Trophic.levelPrimary producer:MetalArsenic    0.2075  0.1863   1.1137  0.2654  -0.1577 
# Trophic.levelHerbivore:MetalCadmium          -0.1309  0.0646  -2.0275  0.0426  -0.2574 
# Trophic.levelPredator:MetalCadmium           -0.0421  0.1227  -0.3433  0.7314  -0.2826 
# Trophic.levelPrimary producer:MetalCadmium    0.3751  0.1392   2.6937  0.0071   0.1022 
# Trophic.levelHerbivore:MetalCobalt           -0.0419  0.1926  -0.2174  0.8279  -0.4194 
# Trophic.levelPredator:MetalCobalt            -0.0979  0.1580  -0.6199  0.5353  -0.4075 
# Trophic.levelPrimary producer:MetalCobalt     0.1610  0.3065   0.5252  0.5995  -0.4398 
# Trophic.levelHerbivore:MetalCopper           -0.1931  0.0587  -3.2908  0.0010  -0.3081 
# Trophic.levelPredator:MetalCopper            -0.2190  0.0824  -2.6556  0.0079  -0.3805 
# Trophic.levelPrimary producer:MetalCopper    -0.0242  0.0871  -0.2776  0.7813  -0.1949 
# Trophic.levelPrimary producer:MetalIron       0.2047  0.1116   1.8348  0.0665  -0.0140 
# Trophic.levelHerbivore:MetalLead             -0.1625  0.1142  -1.4223  0.1549  -0.3864 
# Trophic.levelPredator:MetalLead              -0.8733  0.2522  -3.4623  0.0005  -1.3677 
# Trophic.levelHerbivore:MetalManganese        -0.0108  0.2821  -0.0383  0.9695  -0.5637 
# Trophic.levelPredator:MetalManganese         -0.1539  0.1605  -0.9592  0.3375  -0.4685 
# Trophic.levelHerbivore:MetalMercury           0.0422  0.0888   0.4759  0.6342  -0.1317 
# Trophic.levelPredator:MetalMercury           -0.1595  0.1126  -1.4158  0.1568  -0.3802 
# Trophic.levelHerbivore:MetalNickel           -0.0069  0.1718  -0.0403  0.9678  -0.3436 
# Trophic.levelHerbivore:MetalSilver           -0.1578  0.1420  -1.1112  0.2665  -0.4362 
# Trophic.levelPredator:MetalSilver            -0.1297  0.1993  -0.6510  0.5150  -0.5203 
# Trophic.levelHerbivore:MetalZinc              0.0195  0.1416   0.1377  0.8905  -0.2581 
# Trophic.levelPredator:MetalZinc              -0.0434  0.1364  -0.3186  0.7501  -0.3108 
# Trophic.levelPrimary producer:MetalZinc       0.1841  0.1709   1.0771  0.2814  -0.1509 
#                                               ci.ub      
# Trophic.levelHerbivore:MetalArsenic          0.5348   
# Trophic.levelPrimary producer:MetalArsenic   0.5727   
# Trophic.levelHerbivore:MetalCadmium         -0.0044  * 
# Trophic.levelPredator:MetalCadmium           0.1984   
# Trophic.levelPrimary producer:MetalCadmium   0.6480  ** 
# Trophic.levelHerbivore:MetalCobalt           0.3356   
# Trophic.levelPredator:MetalCobalt            0.2117   
# Trophic.levelPrimary producer:MetalCobalt    0.7618   
# Trophic.levelHerbivore:MetalCopper          -0.0781  *** 
# Trophic.levelPredator:MetalCopper           -0.0574  ** 
# Trophic.levelPrimary producer:MetalCopper    0.1466   
# Trophic.levelPrimary producer:MetalIron      0.4234  . 
# Trophic.levelHerbivore:MetalLead             0.0614   
# Trophic.levelPredator:MetalLead             -0.3790  *** 
# Trophic.levelHerbivore:MetalManganese        0.5421   
# Trophic.levelPredator:MetalManganese         0.1606   
# Trophic.levelHerbivore:MetalMercury          0.2162   
# Trophic.levelPredator:MetalMercury           0.0613   
# Trophic.levelHerbivore:MetalNickel           0.3298   
# Trophic.levelHerbivore:MetalSilver           0.1206   
# Trophic.levelPredator:MetalSilver            0.2608   
# Trophic.levelHerbivore:MetalZinc             0.2971   
# Trophic.levelPredator:MetalZinc              0.2239   
# Trophic.levelPrimary producer:MetalZinc      0.5190   




############## Stressor:Metal #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Metal))
colnames(rep_data) <- c("Stressor","Metal","rep")
data <- merge(data, rep_data, by = c("Stressor","Metal"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Metal - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                             estimate      se     zval    pval    ci.lb    ci.ub     
# StressorA:MetalArsenic      0.3554  0.2079   1.7099  0.0873  -0.0520   0.7628  . 
# StressorW:MetalArsenic     -0.1272  0.2048  -0.6211  0.5345  -0.5287   0.2743   
# StressorA:MetalCadmium      0.0234  0.0649   0.3599  0.7190  -0.1039   0.1506   
# StressorA-W:MetalCadmium   -0.2206  0.2491  -0.8855  0.3759  -0.7088   0.2676   
# StressorH:MetalCadmium     -0.3198  0.1457  -2.1951  0.0282  -0.6053  -0.0343  * 
# StressorW:MetalCadmium     -0.0299  0.1161  -0.2578  0.7965  -0.2574   0.1976   
# StressorA:MetalCobalt      -0.0763  0.1362  -0.5597  0.5757  -0.3433   0.1908   
# StressorW:MetalCobalt       0.2690  0.1934   1.3907  0.1643  -0.1101   0.6481   
# StressorA:MetalCopper      -0.1184  0.0561  -2.1113  0.0347  -0.2282  -0.0085  * 
# StressorA-W:MetalCopper    -0.2887  0.2765  -1.0440  0.2965  -0.8307   0.2533   
# StressorW:MetalCopper      -0.2113  0.0673  -3.1369  0.0017  -0.3433  -0.0793  ** 
# StressorA:MetalIron         0.1407  0.1396   1.0083  0.3133  -0.1328   0.4142   
# StressorW:MetalIron         0.3343  0.1932   1.7304  0.0836  -0.0443   0.7129  . 
# StressorA:MetalLead        -0.1282  0.1749  -0.7331  0.4635  -0.4710   0.2146   
# StressorW:MetalLead        -0.3803  0.1300  -2.9253  0.0034  -0.6351  -0.1255  ** 
# StressorA:MetalManganese    0.1363  0.1781   0.7649  0.4443  -0.2129   0.4854   
# StressorW:MetalManganese    0.0856  0.2253   0.3798  0.7041  -0.3560   0.5272   
# StressorA:MetalMercury      0.1597  0.0980   1.6300  0.1031  -0.0323   0.3517   
# StressorA-W:MetalMercury   -0.1961  0.1768  -1.1092  0.2673  -0.5426   0.1504   
# StressorW:MetalMercury     -0.1994  0.0963  -2.0717  0.0383  -0.3881  -0.0108  * 
# StressorA:MetalNickel       0.1992  0.2414   0.8249  0.4095  -0.2741   0.6724   
# StressorW:MetalNickel      -0.1955  0.2246  -0.8704  0.3841  -0.6358   0.2448   
# StressorA:MetalSilver      -0.0988  0.1301  -0.7589  0.4479  -0.3538   0.1563   
# StressorW:MetalSilver      -0.1010  0.2487  -0.4061  0.6847  -0.5885   0.3865   
# StressorA:MetalZinc         0.0208  0.0940   0.2207  0.8254  -0.1636   0.2051   
# StressorW:MetalZinc         0.0093  0.2515   0.0372  0.9704  -0.4835   0.5022   



############## Stressor:Trait #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Trait))
colnames(rep_data) <- c("Stressor","Trait","rep")
data <- merge(data, rep_data, by = c("Stressor","Trait"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trait - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:TraitAbundance         0.2314  0.1683   1.3748  0.1692  -0.0985   0.5613   
# StressorW:TraitAbundance        -0.0145  0.1505  -0.0967  0.9230  -0.3095   0.2804   
# StressorA:TraitAccumulation     -0.0086  0.0475  -0.1820  0.8555  -0.1018   0.0845   
# StressorA-W:TraitAccumulation   -0.2096  0.1637  -1.2804  0.2004  -0.5303   0.1112   
# StressorH:TraitAccumulation     -0.1567  0.1608  -0.9747  0.3297  -0.4718   0.1584   
# StressorW:TraitAccumulation     -0.2306  0.0752  -3.0671  0.0022  -0.3779  -0.0832  ** 
# StressorA:TraitBehaviour         0.0084  0.1975   0.0424  0.9662  -0.3787   0.3954   
# StressorW:TraitBehaviour        -0.1628  0.1994  -0.8163  0.4143  -0.5537   0.2281   
# StressorA:TraitGrowth            0.0176  0.0627   0.2800  0.7795  -0.1054   0.1405   
# StressorW:TraitGrowth           -0.1284  0.0995  -1.2897  0.1972  -0.3235   0.0667   
# StressorA:TraitMetabolism       -0.0480  0.0541  -0.8868  0.3752  -0.1541   0.0581   
# StressorA-W:TraitMetabolism     -0.5514  0.2200  -2.5062  0.0122  -0.9826  -0.1202  * 
# StressorH:TraitMetabolism       -0.4154  0.2081  -1.9958  0.0460  -0.8233  -0.0075  * 
# StressorW:TraitMetabolism       -0.0182  0.0662  -0.2754  0.7830  -0.1479   0.1115   
# StressorA:TraitReproduction      0.1458  0.1338   1.0817  0.2794  -0.1184   0.4099   
# StressorW:TraitReproduction     -0.4770  0.1435  -3.3251  0.0009  -0.7582  -0.1959  *** 
# StressorA:TraitSurvival         -0.0756  0.1450  -0.5211  0.6023  -0.3598   0.2087   
# StressorW:TraitSurvival         -0.2609  0.1095  -2.3831  0.0172  -0.4754  -0.0463  * 
  
############## Stressor:Trait:Habitat #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Stressor,data$Trait,data$Habitat))
colnames(rep_data) <- c("Stressor","Trait","Habitat","rep")
data <- merge(data, rep_data, by = c("Stressor","Trait","Habitat"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor:Trait:Habitat - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorW:TraitMetabolism:HabitatPolar              0.6582   0.2332   2.8230   0.0048   0.2012   1.1152  ** 
# StressorA:TraitAbundance:HabitatTemperate           0.0544   0.2191   0.2481   0.8040  -0.3750   0.4837   
# StressorA:TraitAccumulation:HabitatTemperate       -0.0234   0.0541  -0.4329   0.6651  -0.1295   0.0827   
# StressorA-W:TraitAccumulation:HabitatTemperate     -0.1365   0.1739  -0.7847   0.4326  -0.4774   0.2044   
# StressorH:TraitAccumulation:HabitatTemperate       -0.3349   0.1622  -2.0641   0.0390  -0.6528  -0.0169  * 
# StressorW:TraitAccumulation:HabitatTemperate       -0.0809   0.0905  -0.8937   0.3715  -0.2583   0.0965   
# StressorA:TraitBehaviour:HabitatTemperate          -0.0540   0.2608  -0.2072   0.8358  -0.5652   0.4571   
# StressorA:TraitGrowth:HabitatTemperate             -0.0028   0.0895  -0.0317   0.9747  -0.1783   0.1726   
# StressorW:TraitGrowth:HabitatTemperate             -0.3348   0.1460  -2.2926   0.0219  -0.6210  -0.0486  * 
# StressorA:TraitMetabolism:HabitatTemperate         -0.0630   0.0671  -0.9391   0.3477  -0.1946   0.0685   
# StressorH:TraitMetabolism:HabitatTemperate         -0.3053   0.2204  -1.3851   0.1660  -0.7372   0.1267   
# StressorW:TraitMetabolism:HabitatTemperate          0.1671   0.0980   1.7041   0.0884  -0.0251   0.3592  . 
# StressorA:TraitSurvival:HabitatTemperate           -0.1964   0.2671  -0.7351   0.4623  -0.7200   0.3272   
# StressorW:TraitSurvival:HabitatTemperate           -0.2287   0.1809  -1.2638   0.2063  -0.5833   0.1260   
# StressorA:TraitAbundance:HabitatTropical           -0.0872   0.2814  -0.3100   0.7566  -0.6388   0.4643   
# StressorW:TraitAbundance:HabitatTropical           -0.1547   0.1532  -1.0096   0.3127  -0.4549   0.1456   
# StressorA:TraitAccumulation:HabitatTropical         0.0501   0.0828   0.6049   0.5452  -0.1122   0.2124   
# StressorW:TraitAccumulation:HabitatTropical        -0.4922   0.1274  -3.8623   0.0001  -0.7419  -0.2424  *** 
# StressorA:TraitBehaviour:HabitatTropical            0.0844   0.2843   0.2969   0.7665  -0.4727   0.6416   
# StressorW:TraitBehaviour:HabitatTropical           -0.2349   0.2282  -1.0294   0.3033  -0.6823   0.2124   
# StressorA:TraitGrowth:HabitatTropical               0.0375   0.0842   0.4457   0.6558  -0.1275   0.2026   
# StressorW:TraitGrowth:HabitatTropical              -0.0197   0.1265  -0.1560   0.8760  -0.2677   0.2282   
# StressorA:TraitMetabolism:HabitatTropical          -0.0350   0.0839  -0.4168   0.6768  -0.1995   0.1295   
# StressorA-W:TraitMetabolism:HabitatTropical        -0.6299   0.2424  -2.5984   0.0094  -1.1051  -0.1548  ** 
# StressorW:TraitMetabolism:HabitatTropical          -0.2555   0.0894  -2.8568   0.0043  -0.4308  -0.0802  ** 
# StressorA:TraitReproduction:HabitatTropical         0.2125   0.1372   1.5488   0.1214  -0.0564   0.4815   
# StressorW:TraitReproduction:HabitatTropical        -0.5033   0.1396  -3.6046   0.0003  -0.7769  -0.2296  *** 
# StressorA:TraitSurvival:HabitatTropical            -0.0147   0.1649  -0.0890   0.9291  -0.3378   0.3084   
# StressorW:TraitSurvival:HabitatTropical            -0.3022   0.1287  -2.3490   0.0188  -0.5544  -0.0501  * 


############## Trophic.level-Accumulation #################### 
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
data <- subset(data, Trait == "Accumulation")
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                 estimate      se     zval    pval    ci.lb    ci.ub 
# Trophic.levelHerbivore          -0.0648  0.0598  -1.0822  0.2791  -0.1821   0.0525   
# Trophic.levelPredator           -0.1778  0.0797  -2.2298  0.0258  -0.3340  -0.0215  * 
# Trophic.levelPrimary producer    0.3518  0.1471   2.3919  0.0168   0.0635   0.6400  * 
# 
# Trophic.levelHerbivore   
# Trophic.levelPredator  * 
# Trophic.levelPrimary producer  * 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                    Estimate Std. Error z value Pr(>|z|)   
# Predator - Herbivore == 0          -0.1130     0.0997  -1.133  0.25699   
# Primary producer - Herbivore == 0   0.4165     0.1515   2.750  0.00596  ** 
# Primary producer - Predator == 0    0.5295     0.1673   3.166  0.00155  ** 

############## Trophic.level-Abundance #################### 
data <- subset(data, Trait == "Abundance")
nrow(data)
data<-subset (data,data$Observation_id!=215)
data<-subset (data,data$Observation_id!=355)
nrow(data)
rep_data <- as.data.frame(table(data$Trophic.level))
colnames(rep_data) <- c("Trophic.level","rep")
data <- merge(data, rep_data, by = c("Trophic.level"))
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Trophic.level - 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   subset = rep >= 3,
                   data = data)
rmod.mix
#                                  estimate      se     zval    pval    ci.lb   ci.ub 
# Trophic.levelHerbivore          -0.1820  0.4154  -0.4381  0.6613  -0.9961  0.6321   
# Trophic.levelPredator           -0.0873  0.3618  -0.2412  0.8094  -0.7964  0.6219   
# Trophic.levelPrimary producer    0.2546  0.2869   0.8875  0.3748  -0.3077  0.8169   
# 
# Trophic.levelHerbivore   
# Trophic.levelPredator   
# Trophic.levelPrimary producer   
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                    Estimate Std. Error z value Pr(>|z|)
# Predator - Herbivore == 0          0.09472    0.55085   0.172    0.863   
# Primary producer - Herbivore == 0  0.43659    0.50479   0.865    0.387   
# Primary producer - Predator == 0   0.34187    0.46175   0.740    0.459   












