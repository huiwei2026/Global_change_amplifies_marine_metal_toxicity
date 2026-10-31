################################################
#                                              #
#Figure S2; Table S4: Likelihood ratio tests   #
#                                              #
################################################
library(metafor)


m1<-rma.mv(yi, vi, mods=~Stressor:Taxa,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#          df      AIC       BIC     AICc    logLik     LRT   pval         QE 
# Full    19 953.3854 1037.1787 954.6779 -457.6927                38382.6813 
# Reduced  3 972.4227  985.6532 972.4624 -483.2113 51.0373 <.0001 41656.8496 

# ΔAIC  -19.0373       



m1<-rma.mv(yi, vi, mods=~Stressor:Habitat,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC       BIC     AICc    logLik     LRT   pval         QE 
# Full    12 948.2625 1001.1846 948.7868 -462.1312                39567.8931 
# Reduced  3 972.4227  985.6532 972.4624 -483.2113 42.1602 <.0001 41656.8496 

# ΔAIC  -24.1602      



m1<-rma.mv(yi, vi, mods=~Stressor:Metal,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#        df      AIC       BIC     AICc    logLik     LRT   pval         QE 
# Full    33 976.0945 1121.6303 980.0039 -455.0472                35451.8399 
# Reduced  3 972.4227  985.6532 972.4624 -483.2113 56.3282 0.0025 41656.8496 

# ΔAIC  3.6718      



m1<-rma.mv(yi, vi, mods=~Stressor:Trait,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC       BIC     AICc    logLik     LRT   pval         QE 
# Full    26 978.7644 1093.4290 981.1810 -463.3822                25884.7231 
# Reduced  3 972.4227  985.6532 972.4624 -483.2113 39.6583 0.0168 41656.8496 

# ΔAIC 6.3417   








m1<-rma.mv(yi, vi, mods=~Stressor:Taxa:Trait,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC       BIC      AICc    logLik      LRT   pval         QE 
# Full    76 980.4702 1315.6435 1002.5117 -414.2351                 15274.0300 
# Reduced  3 972.4227  985.6532  972.4624 -483.2113 137.9524 <.0001 41656.8496 
# ΔAIC  8.0475   




m1<-rma.mv(yi, vi, mods=~Stressor:Taxa:Metal,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC       BIC     AICc    logLik      LRT   pval         QE 
# Full    72 969.4890 1287.0216 989.1376 -412.7445                 31084.3009 
# Reduced  3 972.4227  985.6532 972.4624 -483.2113 140.9336 <.0001 41656.8496 
# ΔAIC  -2.9337  



m1<-rma.mv(yi, vi, mods=~Stressor:Taxa:Habitat,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC       BIC     AICc    logLik     LRT   pval         QE 
# Full    32 952.0340 1093.1596 955.7071 -444.0170                37367.7562 
# Reduced  3 972.4227  985.6532 972.4624 -483.2113 78.3887 <.0001 41656.8496

# ΔAIC  -20.3887  



m1<-rma.mv(yi, vi, mods=~Stressor:Trait:Metal,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC       BIC      AICc    logLik      LRT   pval         QE 
# Full    103 984.1145 1438.3626 1026.6225 -389.0573                 14910.4611 
# Reduced   3 972.4227  985.6532  972.4624 -483.2113 188.3081 <.0001 41656.8496 

# ΔAIC  11.6918   



m1<-rma.mv(yi, vi, mods=~Stressor:Trait:Habitat,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC       BIC     AICc    logLik      LRT   pval         QE 
# Full    49 957.2696 1173.3682 966.0510 -429.6348                 17583.3985 
# Reduced  3 972.4227  985.6532 972.4624 -483.2113 107.1531 <.0001 41656.8496 
# ΔAIC  -15.1531   



m1<-rma.mv(yi, vi, mods=~Stressor:Metal:Habitat,random=list(~1|Paper_id/Observation_id),method="ML",data=data)
m2<-rma.mv(yi, vi, mods=~1,  random=list(~1|Paper_id/Observation_id),method="ML",data=data)
anova(m1,m2)
#         df      AIC       BIC     AICc    logLik     LRT   pval         QE 
# Full    50 976.5137 1197.0224 985.6699 -438.2568                33769.7939 
# Reduced  3 972.4227  985.6532 972.4624 -483.2113 89.9090 0.0002 41656.8496 


# ΔAIC   4.091 
