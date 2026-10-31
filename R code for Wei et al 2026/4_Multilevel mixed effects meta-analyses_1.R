
library(metafor)    
library(multcomp)  
library(ape)        
library(phytools)   
library(dplyr)     
library(Matrix)        
library(parallel) 




# Calculate Blomberg's K
data$Species <- tolower(gsub(" ", "_", data$Species))
tree<-read.tree("phylogeny.tre")
vcov <- ape::vcv(tree)
Species_avg_effect <- aggregate(yi ~ Species, data = data, FUN = mean)
effect_sizes <- Species_avg_effect$yi
names(effect_sizes) <- Species_avg_effect$Species

result_k <- phylosig(tree, effect_sizes, method = "K", test = TRUE, nsim = 10000)
print(result_k )
# Phylogenetic signal K : 0.042012 
# P-value (based on 10000 randomizations) : 0.1866 



################################################
#                                              #
#        Figure 2D; Table S3                   #
#                                              #
################################################
############## overall #################### 
rmod.mix <- rma.mv(yi, vi, 
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
# estimate      se     zval    pval    ci.lb    ci.ub     
# -0.0730  0.0283  -2.5829  0.0098  -0.1285  -0.0176  ** 



############## Stressor #################### 
rmod.mix <- rma.mv(yi, vi, 
                   mods = ~ Stressor- 1,
                   method = "REML",
                   random=list (~1| Paper_id/Observation_id),
                   data = data)
rmod.mix
#              estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA     -0.0002  0.0348  -0.0053  0.9958  -0.0683   0.0680      
# StressorA-W   -0.2459  0.1179  -2.0859  0.0370  -0.4770  -0.0148    * 
# StressorH     -0.2948  0.1250  -2.3596  0.0183  -0.5397  -0.0499    * 
# StressorW     -0.1480  0.0436  -3.3969  0.0007  -0.2334  -0.0626  *** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("A"=1,"A-W"=1,"H"=1,"W"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#               Estimate Std. Error z value Pr(>|z|)   
# A-W - A == 0 -0.24575    0.11938  -2.059  0.03953 * 
# H - A == 0   -0.29466    0.12970  -2.272  0.02310 * 
# W - A == 0   -0.14780    0.05411  -2.731  0.00631 **
# H - A-W == 0 -0.04890    0.17180  -0.285  0.77590   
# W - A-W == 0  0.09796    0.12152   0.806  0.42020   
# W - H == 0    0.14686    0.13233   1.110  0.26708  



############## Taxa ####################
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
# TaxaAlgae+Protist     0.1829  0.0580   3.1518  0.0016   0.0692   0.2966  ** 
# TaxaDeuterostomia    -0.2040  0.0817  -2.4982  0.0125  -0.3641  -0.0440   * 
# TaxaEcdysozoa        -0.1931  0.0604  -3.1964  0.0014  -0.3115  -0.0747  ** 
# TaxaLophotrochozoa   -0.0610  0.0394  -1.5499  0.1212  -0.1381   0.0161     
# TaxaRadiata          -0.2343  0.0810  -2.8912  0.0038  -0.3931  -0.0755  ** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Algae+Protist"=1,"Deuterostomia"=1,"Ecdysozoa"=1,"Lophotrochozoa"=1,"Radiata"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                      Estimate Std. Error z value Pr(>|z|)    
# Deuterostomia - Algae+Protist == 0  -0.38692    0.10018  -3.862 0.000112 ***
# Ecdysozoa - Algae+Protist == 0      -0.37600    0.08374  -4.490 7.12e-06 ***
# Lophotrochozoa - Algae+Protist == 0 -0.24390    0.06870  -3.550 0.000385 ***
# Radiata - Algae+Protist == 0        -0.41718    0.09967  -4.186 2.84e-05 ***
# Ecdysozoa - Deuterostomia == 0       0.01092    0.10141   0.108 0.914222    
# Lophotrochozoa - Deuterostomia == 0  0.14302    0.09058   1.579 0.114343    
# Radiata - Deuterostomia == 0        -0.03026    0.11504  -0.263 0.792542    
# Lophotrochozoa - Ecdysozoa == 0      0.13209    0.07158   1.845 0.064978 .  
# Radiata - Ecdysozoa == 0            -0.04118    0.10107  -0.407 0.683683    
# Radiata - Lophotrochozoa == 0       -0.17327    0.09008  -1.923 0.054421 .  



############## Habitat ####################
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
# HabitatPolar        0.5009  0.1665   3.0092  0.0026   0.1746   0.8271  ** 
# HabitatTemperate   -0.0706  0.0373  -1.8911  0.0586  -0.1437   0.0026   . 
# HabitatTropical    -0.1085  0.0408  -2.6575  0.0079  -0.1885  -0.0285  ** 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Polar"=1,"Temperate"=1,"Tropical"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                            Estimate Std. Error z value Pr(>|z|)    
# Temperate - Polar == 0    -0.57144    0.17058  -3.350 0.000808 ***
# Tropical - Polar == 0     -0.60935    0.17138  -3.555 0.000377 ***
# Tropical - Temperate == 0 -0.03791    0.05527  -0.686 0.492799    


############## Trait ####################
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
# TraitAbundance       0.1000  0.1152   0.8685  0.3851  -0.1257   0.3258    
# TraitAccumulation   -0.0937  0.0425  -2.2037  0.0275  -0.1770  -0.0104  * 
# TraitBehaviour      -0.1347  0.1356  -0.9933  0.3206  -0.4005   0.1311    
# TraitGrowth         -0.0486  0.0555  -0.8758  0.3811  -0.1573   0.0601    
# TraitMetabolism     -0.0620  0.0432  -1.4351  0.1512  -0.1466   0.0227    
# TraitReproduction   -0.1486  0.1015  -1.4640  0.1432  -0.3476   0.0504    
# TraitSurvival       -0.1169  0.0884  -1.3216  0.1863  -0.2902   0.0565     
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Abundance"=1,"Accumulation"=1,"Behaviour"=1,
                                     "Growth"=1,"Metabolism"=1,"Reproduction"=1,"Survival"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                  Estimate Std. Error z value Pr(>|z|)
# Accumulation - Abundance == 0    -0.19372    0.12155  -1.594    0.111
# Behaviour - Abundance == 0       -0.23476    0.17716  -1.325    0.185
# Growth - Abundance == 0          -0.14863    0.12381  -1.200    0.230
# Metabolism - Abundance == 0      -0.16200    0.12011  -1.349    0.177
# Reproduction - Abundance == 0    -0.24866    0.15223  -1.633    0.102
# Survival - Abundance == 0        -0.21691    0.14281  -1.519    0.129
# Behaviour - Accumulation == 0    -0.04104    0.13992  -0.293    0.769
# Growth - Accumulation == 0        0.04509    0.06559   0.687    0.492
# Metabolism - Accumulation == 0    0.03172    0.05566   0.570    0.569
# Reproduction - Accumulation == 0 -0.05494    0.10694  -0.514    0.607
# Survival - Accumulation == 0     -0.02320    0.09509  -0.244    0.807
# Growth - Behaviour == 0           0.08613    0.14260   0.604    0.546
# Metabolism - Behaviour == 0       0.07276    0.14102   0.516    0.606
# Reproduction - Behaviour == 0    -0.01390    0.16719  -0.083    0.934
# Survival - Behaviour == 0         0.01784    0.15800   0.113    0.910
# Metabolism - Growth == 0         -0.01337    0.06522  -0.205    0.838
# Reproduction - Growth == 0       -0.10003    0.10779  -0.928    0.353
# Survival - Growth == 0           -0.06828    0.10015  -0.682    0.495
# Reproduction - Metabolism == 0   -0.08666    0.10811  -0.802    0.423
# Survival - Metabolism == 0       -0.05491    0.09639  -0.570    0.569
# Survival - Reproduction == 0      0.03175    0.12792   0.248    0.804  


############## Metal ####################
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
#                   estimate      se     zval    pval    ci.lb    ci.ub     
# MetalArsenic      0.0979  0.1562   0.6270  0.5307  -0.2082   0.4041     
# MetalCadmium     -0.0390  0.0560  -0.6962  0.4863  -0.1487   0.0708     
# MetalCobalt      -0.0448  0.1212  -0.3694  0.7118  -0.2823   0.1928     
# MetalCopper      -0.1429  0.0456  -3.1315  0.0017  -0.2324  -0.0535  ** 
# MetalIron         0.2125  0.1204   1.7656  0.0775  -0.0234   0.4484   . 
# MetalLead        -0.3504  0.1084  -3.2323  0.0012  -0.5628  -0.1379  ** 
# MetalManganese   -0.1011  0.1411  -0.7164  0.4738  -0.3777   0.1755     
# MetalMercury     -0.0394  0.0755  -0.5223  0.6014  -0.1874   0.1085     
# MetalNickel      -0.0092  0.1767  -0.0523  0.9583  -0.3557   0.3372     
# MetalSilver      -0.1205  0.1190  -1.0123  0.3114  -0.3537   0.1128     
# MetalZinc         0.0340  0.0888   0.3828  0.7019  -0.1400   0.2080        
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Arsenic "=1,"Cadmium"=1,"Cobalt"=1,
                                     "Copper"=1,"Iron"=1,"Lead "=1,"Manganese"=1,
                                     "Mercury"=1,"Nickel"=1,"Silver"=1,"Zinc"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                           Estimate Std. Error z value Pr(>|z|)    
# Cadmium - Arsenic  == 0   -0.1369181  0.1644852  -0.832 0.405181    
# Cobalt - Arsenic  == 0    -0.1426994  0.1975888  -0.722 0.470169    
# Copper - Arsenic  == 0    -0.2408674  0.1624945  -1.482 0.138258    
# Iron - Arsenic  == 0       0.1145606  0.1971839   0.581 0.561252    
# Lead  - Arsenic  == 0     -0.4482916  0.1893136  -2.368 0.017885 *  
# Manganese - Arsenic  == 0 -0.1990443  0.2103944  -0.946 0.344121    
# Mercury - Arsenic  == 0   -0.1373575  0.1730774  -0.794 0.427417    
# Nickel - Arsenic  == 0    -0.1071690  0.2358659  -0.454 0.649567    
# Silver - Arsenic  == 0    -0.2183948  0.1959799  -1.114 0.265119    
# Zinc - Arsenic  == 0      -0.0639487  0.1790728  -0.357 0.721009    
# Cobalt - Cadmium == 0     -0.0057813  0.1321034  -0.044 0.965093    
# Copper - Cadmium == 0     -0.1039494  0.0707379  -1.470 0.141697    
# Iron - Cadmium == 0        0.2514786  0.1327324   1.895 0.058142 .  
# Lead  - Cadmium == 0      -0.3113736  0.1211433  -2.570 0.010161 *  
# Manganese - Cadmium == 0  -0.0621262  0.1504208  -0.413 0.679595    
# Mercury - Cadmium == 0    -0.0004395  0.0936433  -0.005 0.996255    
# Nickel - Cadmium == 0      0.0297490  0.1853903   0.160 0.872513    
# Silver - Cadmium == 0     -0.0814768  0.1247312  -0.653 0.513615    
# Zinc - Cadmium == 0        0.0729694  0.1037220   0.704 0.481739    
# Copper - Cobalt == 0      -0.0981681  0.1292488  -0.760 0.447537    
# Iron - Cobalt == 0         0.2572599  0.1704214   1.510 0.131158    
# Lead  - Cobalt == 0       -0.3055923  0.1624921  -1.881 0.060018 .  
# Manganese - Cobalt == 0   -0.0563449  0.1522960  -0.370 0.711405    
# Mercury - Cobalt == 0      0.0053418  0.1420317   0.038 0.969999    
# Nickel - Cobalt == 0       0.0355303  0.2142433   0.166 0.868282    
# Silver - Cobalt == 0      -0.0756955  0.1654836  -0.457 0.647369    
# Zinc - Cobalt == 0         0.0787507  0.1467403   0.537 0.591498    
# Iron - Copper == 0         0.3554280  0.1287052   2.762 0.005752 ** 
# Lead  - Copper == 0       -0.2074242  0.1159621  -1.789 0.073659 .  
# Manganese - Copper == 0    0.0418231  0.1480563   0.282 0.777575    
# Mercury - Copper == 0      0.1035099  0.0880709   1.175 0.239874    
# Nickel - Copper == 0       0.1336984  0.1825122   0.733 0.463836    
# Silver - Copper == 0       0.0224726  0.1262884   0.178 0.858765    
# Zinc - Copper == 0         0.1769187  0.0972966   1.818 0.069012 .  
# Lead  - Iron == 0         -0.5628522  0.1619645  -3.475 0.000511 ***
# Manganese - Iron == 0     -0.3136049  0.1853400  -1.692 0.090636 .  
# Mercury - Iron == 0       -0.2519181  0.1420545  -1.773 0.076164 .  
# Nickel - Iron == 0        -0.2217296  0.2138282  -1.037 0.299758    
# Silver - Iron == 0        -0.3329554  0.1692128  -1.968 0.049106 *  
# Zinc - Iron == 0          -0.1785093  0.1492335  -1.196 0.231628    
# Manganese - Lead  == 0     0.2492473  0.1778466   1.401 0.161072    
# Mercury - Lead  == 0       0.3109341  0.1317954   2.359 0.018314 *  
# Nickel - Lead  == 0        0.3411226  0.2073224   1.645 0.099893 .  
# Silver - Lead  == 0        0.2298968  0.1606722   1.431 0.152475    
# Zinc - Lead  == 0          0.3843429  0.1391466   2.762 0.005742 ** 
# Mercury - Manganese == 0   0.0616868  0.1589808   0.388 0.698006    
# Nickel - Manganese == 0    0.0918753  0.2261106   0.406 0.684501    
# Silver - Manganese == 0   -0.0193505  0.1795441  -0.108 0.914173    
# Zinc - Manganese == 0      0.1350956  0.1625137   0.831 0.405811    
# Nickel - Mercury == 0      0.0301885  0.1921762   0.157 0.875176    
# Silver - Mercury == 0     -0.0810373  0.1399135  -0.579 0.562457    
# Zinc - Mercury == 0        0.0734088  0.1156722   0.635 0.525671    
# Silver - Nickel == 0      -0.1112258  0.2130107  -0.522 0.601559    
# Zinc - Nickel == 0         0.0432203  0.1967689   0.220 0.826144    
# Zinc - Silver == 0         0.1544461  0.1448900   1.066 0.286444       


################################################
#                                              #
#                   Figure 4                   #
#                                              #
################################################
############## Stressor:Taxa_Table S10 #################### 
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
# StressorA:TaxaAlgae+Protist       0.1775  0.0603   2.9433  0.0032   0.0593   0.2957   ** 
# StressorA-W:TaxaAlgae+Protist    -0.1087  0.2767  -0.3927  0.6945  -0.6511   0.4337      
# StressorW:TaxaAlgae+Protist       0.2462  0.1299   1.8956  0.0580  -0.0084   0.5008    . 
# StressorA:TaxaDeuterostomia      -0.1142  0.1013  -1.1271  0.2597  -0.3129   0.0844      
# StressorA-W:TaxaDeuterostomia    -0.1911  0.2047  -0.9335  0.3505  -0.5923   0.2101      
# StressorW:TaxaDeuterostomia      -0.3381  0.1188  -2.8467  0.0044  -0.5710  -0.1053   ** 
# StressorA:TaxaEcdysozoa           0.0555  0.0876   0.6338  0.5262  -0.1162   0.2272      
# StressorH:TaxaEcdysozoa          -0.4380  0.1865  -2.3485  0.0188  -0.8035  -0.0725    * 
# StressorW:TaxaEcdysozoa          -0.3308  0.0818  -4.0425  <.0001  -0.4911  -0.1704  *** 
# StressorA:TaxaLophotrochozoa     -0.0530  0.0465  -1.1389  0.2547  -0.1441   0.0382      
# StressorA-W:TaxaLophotrochozoa   -0.2166  0.1795  -1.2069  0.2275  -0.5683   0.1351      
# StressorH:TaxaLophotrochozoa     -0.2218  0.1673  -1.3253  0.1851  -0.5498   0.1062      
# StressorW:TaxaLophotrochozoa     -0.0344  0.0621  -0.5539  0.5797  -0.1561   0.0873      
# StressorA:TaxaRadiata            -0.1731  0.1191  -1.4541  0.1459  -0.4065   0.0602      
# StressorW:TaxaRadiata            -0.3117  0.1030  -3.0276  0.0025  -0.5135  -0.1099   ** 




############## Stressor:Habitat_Table S11 #################### 
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
# StressorA:HabitatPolar          0.4837  0.2432   1.9891  0.0467   0.0071   0.9603    * 
# StressorW:HabitatPolar          0.5136  0.1987   2.5842  0.0098   0.1241   0.9031   ** 
# StressorA:HabitatTemperate     -0.0347  0.0417  -0.8303  0.4064  -0.1165   0.0472      
# StressorA-W:HabitatTemperate   -0.0982  0.1469  -0.6687  0.5037  -0.3862   0.1897      
# StressorH:HabitatTemperate     -0.3578  0.1245  -2.8743  0.0040  -0.6018  -0.1138   ** 
# StressorW:HabitatTemperate     -0.0534  0.0604  -0.8831  0.3772  -0.1718   0.0651      
# StressorA:HabitatTropical       0.0425  0.0498   0.8534  0.3934  -0.0551   0.1400      
# StressorA-W:HabitatTropical    -0.4728  0.1834  -2.5783  0.0099  -0.8322  -0.1134   ** 
# StressorW:HabitatTropical      -0.2887  0.0561  -5.1460  <.0001  -0.3987  -0.1787  *** 

############## Stressor:Habitat:Taxa_Table S12 #################### 
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
#                                                   estimate      se     zval    pval    ci.lb    ci.ub      
# StressorA:HabitatPolar:TaxaAlgae+Protist           0.4836  0.2400   2.0147  0.0439   0.0131   0.9541    * 
# StressorW:HabitatPolar:TaxaAlgae+Protist           0.7030  0.2722   2.5824  0.0098   0.1694   1.2366   ** 
# StressorA:HabitatTemperate:TaxaAlgae+Protist       0.1214  0.0834   1.4553  0.1456  -0.0421   0.2849      
# StressorW:HabitatTemperate:TaxaAlgae+Protist       0.2860  0.2122   1.3479  0.1777  -0.1299   0.7020      
# StressorA:HabitatTropical:TaxaAlgae+Protist        0.1901  0.0918   2.0710  0.0384   0.0102   0.3700    * 
# StressorA-W:HabitatTropical:TaxaAlgae+Protist     -0.1697  0.2767  -0.6135  0.5396  -0.7120   0.3726      
# StressorW:HabitatTropical:TaxaAlgae+Protist       -0.0356  0.2013  -0.1766  0.8598  -0.4301   0.3590      
# StressorA:HabitatTemperate:TaxaDeuterostomia      -0.1382  0.1009  -1.3703  0.1706  -0.3359   0.0595      
# StressorA-W:HabitatTemperate:TaxaDeuterostomia    -0.2031  0.2032  -0.9993  0.3177  -0.6013   0.1952      
# StressorW:HabitatTemperate:TaxaDeuterostomia      -0.3517  0.1484  -2.3697  0.0178  -0.6425  -0.0608    * 
# StressorW:HabitatTropical:TaxaDeuterostomia       -0.3212  0.1924  -1.6699  0.0949  -0.6982   0.0558    . 
# StressorA:HabitatTemperate:TaxaEcdysozoa          -0.1869  0.2381  -0.7850  0.4324  -0.6535   0.2797      
# StressorH:HabitatTemperate:TaxaEcdysozoa          -0.4396  0.1845  -2.3827  0.0172  -0.8011  -0.0780    * 
# StressorW:HabitatTemperate:TaxaEcdysozoa          -0.2408  0.1936  -1.2437  0.2136  -0.6202   0.1386      
# StressorA:HabitatTropical:TaxaEcdysozoa            0.0922  0.0933   0.9885  0.3229  -0.0907   0.2751      
# StressorW:HabitatTropical:TaxaEcdysozoa           -0.3482  0.0893  -3.8984  <.0001  -0.5233  -0.1731  *** 
# StressorW:HabitatPolar:TaxaLophotrochozoa          0.3063  0.2832   1.0816  0.2794  -0.2487   0.8614      
# StressorA:HabitatTemperate:TaxaLophotrochozoa     -0.0550  0.0525  -1.0474  0.2949  -0.1580   0.0479      
# StressorA-W:HabitatTemperate:TaxaLophotrochozoa   -0.0791  0.2096  -0.3775  0.7058  -0.4900   0.3317      
# StressorH:HabitatTemperate:TaxaLophotrochozoa     -0.3577  0.1864  -1.9186  0.0550  -0.7231   0.0077    . 
# StressorW:HabitatTemperate:TaxaLophotrochozoa     -0.0093  0.0710  -0.1316  0.8953  -0.1485   0.1298      
# StressorA:HabitatTropical:TaxaLophotrochozoa      -0.0790  0.0971  -0.8133  0.4160  -0.2693   0.1113      
# StressorA-W:HabitatTropical:TaxaLophotrochozoa    -0.5796  0.3316  -1.7480  0.0805  -1.2295   0.0703    . 
# StressorW:HabitatTropical:TaxaLophotrochozoa      -0.1959  0.1380  -1.4193  0.1558  -0.4665   0.0746      
# StressorA:HabitatTropical:TaxaRadiata             -0.1354  0.1174  -1.1538  0.2486  -0.3654   0.0946      
# StressorW:HabitatTropical:TaxaRadiata             -0.3114  0.1020  -3.0526  0.0023  -0.5114  -0.1115   ** 
  

################################################
#                                              #
#                    Figure 5                  #
#                                              #
################################################
####### Trophic.level #########
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
# Trophic.levelHerbivore          -0.0936  0.0360  -2.6015  0.0093  -0.1641  -0.0231   ** 
# Trophic.levelPredator           -0.1731  0.0507  -3.4173  0.0006  -0.2725  -0.0738  *** 
# Trophic.levelPrimary producer    0.1186  0.0566   2.0959  0.0361   0.0077   0.2294    * 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                     Estimate Std. Error z value Pr(>|z|)    
# Predator - Herbivore == 0         -0.07953    0.06157  -1.292  0.19640    
# Primary producer - Herbivore == 0  0.21217    0.06584   3.223  0.00127 ** 
# Primary producer - Predator == 0   0.29171    0.07453   3.914 9.09e-05 ***

######### Stressor:Trophic.level_Figure 5C_Table S14 #######
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
#                                            estimate      se     zval    pval    ci.lb    ci.ub     
# StressorA:Trophic.levelHerbivore          -0.0191  0.0476  -0.4016  0.6880  -0.1125   0.0742     
# StressorA-W:Trophic.levelHerbivore        -0.2585  0.1893  -1.3659  0.1720  -0.6294   0.1124     
# StressorH:Trophic.levelHerbivore          -0.1848  0.1632  -1.1321  0.2576  -0.5047   0.1351     
# StressorW:Trophic.levelHerbivore          -0.1766  0.0553  -3.1908  0.0014  -0.2851  -0.0681  ** 
# StressorA:Trophic.levelPredator           -0.0970  0.0668  -1.4507  0.1469  -0.2280   0.0340     
# StressorA-W:Trophic.levelPredator         -0.2925  0.1511  -1.9361  0.0529  -0.5887   0.0036   . 
# StressorH:Trophic.levelPredator           -0.4357  0.1868  -2.3323  0.0197  -0.8019  -0.0696   * 
# StressorW:Trophic.levelPredator           -0.2066  0.0731  -2.8243  0.0047  -0.3499  -0.0632  ** 
# StressorA:Trophic.levelPrimary producer    0.1232  0.0628   1.9625  0.0497   0.0002   0.2463   * 
# StressorW:Trophic.levelPrimary producer    0.1280  0.1236   1.0355  0.3005  -0.1143   0.3702     



######### Trophic.level:Metal_Figure 5B_Table S15 #########
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
#                                              estimate      se     zval    pval    ci.lb    ci.ub      
# Trophic.levelHerbivore:MetalArsenic           0.1822  0.1842   0.9893  0.3225  -0.1788   0.5432      
# Trophic.levelPrimary producer:MetalArsenic    0.2045  0.1919   1.0659  0.2865  -0.1715   0.5806      
# Trophic.levelHerbivore:MetalCadmium          -0.1302  0.0660  -1.9708  0.0487  -0.2596  -0.0007    * 
# Trophic.levelPredator:MetalCadmium           -0.0411  0.1256  -0.3271  0.7436  -0.2872   0.2050      
# Trophic.levelPrimary producer:MetalCadmium    0.3744  0.1424   2.6291  0.0086   0.0953   0.6535   ** 
# Trophic.levelHerbivore:MetalCobalt           -0.0434  0.1960  -0.2215  0.8247  -0.4276   0.3408      
# Trophic.levelPredator:MetalCobalt            -0.1000  0.1617  -0.6185  0.5362  -0.4169   0.2169      
# Trophic.levelPrimary producer:MetalCobalt     0.1623  0.3123   0.5195  0.6034  -0.4499   0.7744      
# Trophic.levelHerbivore:MetalCopper           -0.1729  0.0599  -2.8850  0.0039  -0.2904  -0.0554   ** 
# Trophic.levelPredator:MetalCopper            -0.2126  0.0841  -2.5267  0.0115  -0.3775  -0.0477    * 
# Trophic.levelPrimary producer:MetalCopper    -0.0217  0.0890  -0.2437  0.8075  -0.1961   0.1527      
# Trophic.levelPrimary producer:MetalIron       0.2065  0.1145   1.8045  0.0712  -0.0178   0.4308    . 
# Trophic.levelHerbivore:MetalLead             -0.1694  0.1166  -1.4522  0.1464  -0.3980   0.0592      
# Trophic.levelPredator:MetalLead              -1.1298  0.2441  -4.6277  <.0001  -1.6084  -0.6513  *** 
# Trophic.levelHerbivore:MetalManganese        -0.0118  0.2868  -0.0410  0.9673  -0.5738   0.5503      
# Trophic.levelPredator:MetalManganese         -0.1550  0.1642  -0.9439  0.3452  -0.4767   0.1668      
# Trophic.levelHerbivore:MetalMercury           0.0414  0.0911   0.4540  0.6498  -0.1372   0.2199      
# Trophic.levelPredator:MetalMercury           -0.1607  0.1157  -1.3886  0.1650  -0.3875   0.0661      
# Trophic.levelHerbivore:MetalNickel           -0.0061  0.1771  -0.0343  0.9727  -0.3531   0.3410      
# Trophic.levelHerbivore:MetalSilver           -0.1550  0.1447  -1.0715  0.2840  -0.4386   0.1286      
# Trophic.levelPredator:MetalSilver            -0.1314  0.2029  -0.6475  0.5173  -0.5291   0.2663      
# Trophic.levelHerbivore:MetalZinc              0.0217  0.1444   0.1501  0.8807  -0.2613   0.3046      
# Trophic.levelPredator:MetalZinc              -0.0425  0.1394  -0.3046  0.7607  -0.3156   0.2307      
# Trophic.levelPrimary producer:MetalZinc       0.1837  0.1747   1.0512  0.2932  -0.1588   0.5262   


############## Habitat:Trophic.level_Figure 5D_Table S16 #################### 
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
#                                                   estimate      se     zval    pval    ci.lb    ci.ub      
# HabitatPolar:Trophic.levelHerbivore               0.3108  0.3011   1.0323  0.3019  -0.2793   0.9009      
# HabitatTemperate:Trophic.levelHerbivore          -0.0951  0.0461  -2.0598  0.0394  -0.1855  -0.0046    * 
# HabitatTropical:Trophic.levelHerbivore           -0.1039  0.0550  -1.8902  0.0587  -0.2117   0.0038    . 
# HabitatTemperate:Trophic.levelPredator           -0.1131  0.0677  -1.6716  0.0946  -0.2457   0.0195    . 
# HabitatTropical:Trophic.levelPredator            -0.2449  0.0729  -3.3593  0.0008  -0.3878  -0.1020  *** 
# HabitatPolar:Trophic.levelPrimary producer        0.5780  0.1907   3.0317  0.0024   0.2043   0.9517   ** 
# HabitatTemperate:Trophic.levelPrimary producer    0.1057  0.0842   1.2560  0.2091  -0.0593   0.2707      
# HabitatTropical:Trophic.levelPrimary producer     0.0476  0.0801   0.5949  0.5519  -0.1093   0.2046    

################################################
#                                              #
#           Figure S6C; Table S17              #
#                                              #
################################################
############## Stressor:Metal #################### 
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
# StressorA:MetalArsenic      0.3556  0.2182   1.6296  0.1032  -0.0721   0.7834      
# StressorW:MetalArsenic     -0.1311  0.2103  -0.6231  0.5332  -0.5433   0.2812      
# StressorA:MetalCadmium      0.0251  0.0668   0.3754  0.7073  -0.1058   0.1560      
# StressorA-W:MetalCadmium   -0.2177  0.2539  -0.8574  0.3912  -0.7154   0.2800      
# StressorH:MetalCadmium     -0.3207  0.1497  -2.1424  0.0322  -0.6141  -0.0273    * 
# StressorW:MetalCadmium     -0.0304  0.1191  -0.2554  0.7984  -0.2638   0.2030      
# StressorA:MetalCobalt      -0.0769  0.1395  -0.5512  0.5815  -0.3504   0.1966      
# StressorW:MetalCobalt       0.2716  0.1980   1.3713  0.1703  -0.1166   0.6598      
# StressorA:MetalCopper      -0.1005  0.0576  -1.7449  0.0810  -0.2133   0.0124    . 
# StressorA-W:MetalCopper    -0.2855  0.2819  -1.0128  0.3111  -0.8379   0.2669      
# StressorW:MetalCopper      -0.2018  0.0693  -2.9103  0.0036  -0.3377  -0.0659   ** 
# StressorA:MetalIron         0.1439  0.1441   0.9984  0.3181  -0.1386   0.4263      
# StressorW:MetalIron         0.3365  0.2002   1.6807  0.0928  -0.0559   0.7289    . 
# StressorA:MetalLead        -0.1293  0.1793  -0.7212  0.4708  -0.4806   0.2221      
# StressorW:MetalLead        -0.4656  0.1320  -3.5279  0.0004  -0.7242  -0.2069  *** 
# StressorA:MetalManganese    0.1378  0.1819   0.7573  0.4489  -0.2188   0.4944      
# StressorW:MetalManganese    0.0873  0.2303   0.3793  0.7045  -0.3640   0.5387      
# StressorA:MetalMercury      0.1594  0.1013   1.5739  0.1155  -0.0391   0.3579      
# StressorA-W:MetalMercury   -0.1962  0.1806  -1.0863  0.2773  -0.5501   0.1578      
# StressorW:MetalMercury     -0.2004  0.0990  -2.0253  0.0428  -0.3944  -0.0065    * 
# StressorA:MetalNickel       0.1972  0.2491   0.7916  0.4286  -0.2910   0.6853      
# StressorW:MetalNickel      -0.1965  0.2348  -0.8369  0.4027  -0.6568   0.2637      
# StressorA:MetalSilver      -0.0979  0.1329  -0.7361  0.4617  -0.3584   0.1627      
# StressorW:MetalSilver      -0.0911  0.2541  -0.3583  0.7201  -0.5892   0.4070      
# StressorA:MetalZinc         0.0223  0.0965   0.2307  0.8175  -0.1669   0.2114      
# StressorW:MetalZinc         0.0067  0.2569   0.0261  0.9791  -0.4967   0.5101     



################################################
#                                              #
#         Figure S4; Tables S5-S6              #
#                                              #
################################################
######### Stressor:Trait_Table S5 #########
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
# StressorA:TraitAbundance         0.2253  0.1715   1.3134  0.1890  -0.1109   0.5614      
# StressorW:TraitAbundance         0.0081  0.1534   0.0529  0.9578  -0.2925   0.3087      
# StressorA:TraitAccumulation     -0.0072  0.0486  -0.1479  0.8824  -0.1024   0.0880      
# StressorA-W:TraitAccumulation   -0.2167  0.1666  -1.3004  0.1935  -0.5433   0.1099      
# StressorH:TraitAccumulation     -0.1562  0.1643  -0.9506  0.3418  -0.4782   0.1659      
# StressorW:TraitAccumulation     -0.2589  0.0764  -3.3889  0.0007  -0.4086  -0.1092  *** 
# StressorA:TraitBehaviour         0.0114  0.2008   0.0570  0.9546  -0.3821   0.4050      
# StressorW:TraitBehaviour        -0.1632  0.2032  -0.8032  0.4219  -0.5615   0.2351      
# StressorA:TraitGrowth            0.0164  0.0640   0.2557  0.7982  -0.1091   0.1419      
# StressorW:TraitGrowth           -0.1280  0.1016  -1.2600  0.2077  -0.3272   0.0711      
# StressorA:TraitMetabolism       -0.0472  0.0552  -0.8553  0.3924  -0.1555   0.0610      
# StressorA-W:TraitMetabolism     -0.5498  0.2236  -2.4587  0.0139  -0.9881  -0.1115    * 
# StressorH:TraitMetabolism       -0.4156  0.2120  -1.9603  0.0500  -0.8311  -0.0001    * 
# StressorW:TraitMetabolism       -0.0154  0.0675  -0.2281  0.8195  -0.1477   0.1169      
# StressorA:TraitReproduction      0.1551  0.1373   1.1292  0.2588  -0.1141   0.4242      
# StressorW:TraitReproduction     -0.4781  0.1462  -3.2707  0.0011  -0.7647  -0.1916   ** 
# StressorA:TraitSurvival          0.0630  0.1446   0.4353  0.6634  -0.2205   0.3464      
# StressorW:TraitSurvival         -0.2619  0.1118  -2.3427  0.0191  -0.4811  -0.0428    * 

############## Stressor:Trait:Habitat_Table S6 ########## 
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
# StressorW:TraitMetabolism:HabitatPolar            0.6612  0.2380   2.7777  0.0055   0.1946   1.1277   ** 
# StressorA:TraitAbundance:HabitatTemperate         0.0579  0.2235   0.2593  0.7954  -0.3801   0.4960      
# StressorA:TraitAccumulation:HabitatTemperate     -0.0241  0.0553  -0.4357  0.6630  -0.1325   0.0843      
# StressorA-W:TraitAccumulation:HabitatTemperate   -0.1377  0.1774  -0.7763  0.4375  -0.4853   0.2099      
# StressorH:TraitAccumulation:HabitatTemperate     -0.3353  0.1657  -2.0234  0.0430  -0.6601  -0.0105    * 
# StressorW:TraitAccumulation:HabitatTemperate     -0.0811  0.0923  -0.8779  0.3800  -0.2620   0.0999      
# StressorA:TraitBehaviour:HabitatTemperate        -0.0553  0.2652  -0.2087  0.8347  -0.5750   0.4644      
# StressorA:TraitGrowth:HabitatTemperate           -0.0036  0.0914  -0.0398  0.9683  -0.1827   0.1755      
# StressorW:TraitGrowth:HabitatTemperate           -0.3348  0.1492  -2.2444  0.0248  -0.6272  -0.0424    * 
# StressorA:TraitMetabolism:HabitatTemperate       -0.0626  0.0685  -0.9145  0.3604  -0.1968   0.0716      
# StressorH:TraitMetabolism:HabitatTemperate       -0.3043  0.2245  -1.3553  0.1753  -0.7443   0.1357      
# StressorW:TraitMetabolism:HabitatTemperate        0.1671  0.1001   1.6703  0.0949  -0.0290   0.3632    . 
# StressorA:TraitSurvival:HabitatTemperate         -0.2006  0.2724  -0.7365  0.4614  -0.7345   0.3332      
# StressorW:TraitSurvival:HabitatTemperate         -0.2297  0.1855  -1.2386  0.2155  -0.5933   0.1338      
# StressorA:TraitAbundance:HabitatTropical         -0.1066  0.2864  -0.3722  0.7098  -0.6680   0.4548      
# StressorW:TraitAbundance:HabitatTropical         -0.1361  0.1561  -0.8717  0.3834  -0.4420   0.1699      
# StressorA:TraitAccumulation:HabitatTropical       0.0577  0.0845   0.6820  0.4952  -0.1080   0.2233      
# StressorW:TraitAccumulation:HabitatTropical      -0.5726  0.1284  -4.4596  <.0001  -0.8243  -0.3210  *** 
# StressorA:TraitBehaviour:HabitatTropical          0.0966  0.2894   0.3336  0.7387  -0.4707   0.6639      
# StressorW:TraitBehaviour:HabitatTropical         -0.2376  0.2330  -1.0198  0.3078  -0.6943   0.2190      
# StressorA:TraitGrowth:HabitatTropical             0.0399  0.0860   0.4645  0.6423  -0.1286   0.2084      
# StressorW:TraitGrowth:HabitatTropical            -0.0228  0.1290  -0.1763  0.8600  -0.2757   0.2301      
# StressorA:TraitMetabolism:HabitatTropical        -0.0328  0.0856  -0.3837  0.7012  -0.2006   0.1349      
# StressorA-W:TraitMetabolism:HabitatTropical      -0.6272  0.2466  -2.5440  0.0110  -1.1105  -0.1440    * 
# StressorW:TraitMetabolism:HabitatTropical        -0.2490  0.0911  -2.7335  0.0063  -0.4276  -0.0705   ** 
# StressorA:TraitReproduction:HabitatTropical       0.2270  0.1399   1.6225  0.1047  -0.0472   0.5013      
# StressorW:TraitReproduction:HabitatTropical      -0.5085  0.1424  -3.5720  0.0004  -0.7875  -0.2295  *** 
# StressorA:TraitSurvival:HabitatTropical           0.1655  0.1634   1.0129  0.3111  -0.1547   0.4857      
# StressorW:TraitSurvival:HabitatTropical          -0.3077  0.1312  -2.3462  0.0190  -0.5648  -0.0507    * 


################################################
#                                              #
#              Figure S7; Table S13            #
#                                              #
################################################
###### Trophic.level-Accumulation ######
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
# Trophic.levelHerbivore          -0.0637  0.0625  -1.0185  0.3084  -0.1862   0.0588    
# Trophic.levelPredator           -0.2017  0.0837  -2.4103  0.0159  -0.3657  -0.0377  * 
# Trophic.levelPrimary producer    0.3452  0.1517   2.2754  0.0229   0.0479   0.6426  * 
summary(glht(rmod.mix, 
             linfct=cbind(contrMat(c("Herbivore"=1,"Predator"=1,"Primary producer"=1), 
                                   type="Tukey"))), test=adjusted("none"))
#                                    Estimate Std. Error z value Pr(>|z|)   
# Predator - Herbivore == 0          -0.1380     0.1044  -1.321   0.1863   
# Primary producer - Herbivore == 0   0.4089     0.1559   2.623   0.0087 **
# Primary producer - Predator == 0    0.5469     0.1733   3.157   0.0016 **

############## Trophic.level-Abundance #################### 
data <- subset(data, Trait == "Abundance")
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
#                                    Estimate Std. Error z value Pr(>|z|)
# Predator - Herbivore == 0          0.09472    0.55085   0.172    0.863
# Primary producer - Herbivore == 0  0.43659    0.50479   0.865    0.387
# Primary producer - Predator == 0   0.34187    0.46175   0.740    0.459






