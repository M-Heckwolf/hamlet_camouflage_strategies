# Melanie J Heckwolf, 2024
# Analysis of activity in experimental tank (open field trial)

# Manuscript:
# Differences in colour pattern, behaviour and gene expression in the brain suggest divergent camouflage strategies in sympatric reef fish species
# MJ Heckwolf, J Gismann, M González-Santoro, F Coulmance, J Fuß, WO McMillan, O Puebla


#load libraries:
library(lme4)
library(car)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(ggsignif)
theme_set(theme_cowplot())

# set your working directory:
setwd("PATH/TO/YOUR/DATA")

# load data:
df <- read.csv2("activity_assessment_data.csv",header=T)



### ------------------ analyses ---------------------------------------------------


#-------------- Distance moved in cm (Suppl Fig 2A):
M1=lm(df$Distance_moved_cm~df$species)
summary(M1)
anova(M1)

# model validation visual:
par(mfrow = c(2, 2))
plot(M1)

# model validation tests (both should be >0.05):
leveneTest(residuals(M1)~data$species) #0.3334
shapiro.test(residuals(M1)) # <0.05
### -> residuals deviate from normal distribution.

# try the following transformations:
# x**2     # skews data more in the wrong direction
# 1/x      # skews data more in the wrong direction
# log()    # corrects in the right direction


# rerun model with log transformation:
M1=lm(log(df$Distance_moved_cm)~df$species)
summary(M1)
anova(M1)

# model validation visual:
par(mfrow = c(2, 2))
plot(M1)

# model validation tests (both should be >0.05):
leveneTest(residuals(M1)~data$species) #0.006
shapiro.test(residuals(M1)) #0.002
# tests are robust to slight deviations from normality and 
# both, original data and log() transformated data, are non significant, so our result is robust.





#-------------- Distance moved in body length (Suppl Fig 2B):
M2=lm(df$Distance_per_body_length~df$species)
summary(M2)
anova(M2)

# model validation visual:
par(mfrow = c(2, 2))
plot(M2)

# model validation tests (both should be >0.05):
leveneTest(residuals(M2)~data$species) #0.2904
shapiro.test(residuals(M2)) #<0.05


# try the following transformations:
# x**2     # skews data more in the wrong direction
# 1/x      # skews data more in the wrong direction
# log()    # corrects in the right direction


# rerun model with log transformation:
M2=lm(log(df$Distance_per_body_length)~df$species)
summary(M2)
anova(M2)

# model validation visual:
par(mfrow = c(2, 2))
plot(M2)

# model validation tests (both should be >0.05):
leveneTest(residuals(M2)~data$species) #0.006
shapiro.test(residuals(M2)) #0.002
# tests are robust to slight deviations from normality and 
# both, original data and log() transformated data, are non significant, so our result is robust.







#-------------- Velocity in cm per second (Suppl Fig 2C):
M3=lm(df$Velocity_cm.s~df$species)
summary(M3)
anova(M3)

# model validation visual:
par(mfrow = c(2, 2))
plot(M3)

# model validation tests (both should be >0.05):
leveneTest(residuals(M3)~data$species) #0.3321
shapiro.test(residuals(M3)) #<0.05



# try the following transformations:
# x**2     # skews data more in the wrong direction
# 1/x      # skews data more in the wrong direction
# log()    # corrects in the right direction


# rerun model with log transformation:
M3=lm(log(df$Velocity_cm.s)~df$species)
summary(M3)
anova(M3)

# model validation visual:
par(mfrow = c(2, 2))
plot(M3)

# model validation tests (both should be >0.05):
leveneTest(residuals(M3)~data$species) #0.006
shapiro.test(residuals(M3)) #0.002
# tests are robust to slight deviations from normality and 
# both, original data and log() transformated data, are non significant, so our result is robust.






#-------------- Velocity in body length per second (Suppl Fig 2D):
M4=lm(df$Velocity_body_length.s~df$species)
summary(M4)
anova(M4)

# model validation visual:
par(mfrow = c(2, 2))
plot(M4)

# model validation tests (both should be >0.05):
leveneTest(residuals(M4)~data$species) #0.2892
shapiro.test(residuals(M4)) #<0.05



# try the following transformations:
# x**2     # skews data more in the wrong direction
# 1/x      # skews data more in the wrong direction
# log()    # corrects in the right direction


# rerun model with log transformation:
M4=lm(log(df$Velocity_body_length.s)~df$species)
summary(M4)
anova(M4)

# model validation visual:
par(mfrow = c(2, 2))
plot(M4)

# model validation tests (both should be >0.05):
leveneTest(residuals(M4)~data$species) #0.006
shapiro.test(residuals(M4)) #0.003
# tests are robust to slight deviations from normality and 
# both, original data and log() transformated data, are non significant, so our result is robust.





### ------------------ plot data ---------------------------------------------------


### Supplementary Figure S2

anova(M1) # take values for annotations
A=ggplot(df,aes(x=species,y=Distance_moved_cm,fill=species,col=species))+
  geom_violin(alpha=0.25, width=.9,color=NA)+
  geom_boxplot(width=0.25)+
  scale_fill_manual(values = c("burlywood","grey28"))+
  scale_color_manual(values = c("black","black"))+
  geom_signif(comparisons = list(c("black","barred")), annotations="F(1,105)=0.615, P=0.435",
              y_position = 12500, tip_length = 0, vjust=0.1)+
  scale_x_discrete(labels = c("barred","black"),name="Species")+
  ylab(label = "distance moved (cm)")+
  theme_classic()+
  theme(legend.position = "none")


anova(M2) # take values for annotations
B=ggplot(df,aes(x=species,y=Distance_per_body_length,fill=species,col=species))+
  geom_violin(alpha=0.25, width=.9,color=NA)+
  geom_boxplot(width=0.25)+
  scale_fill_manual(values = c("burlywood","grey28"))+
  scale_color_manual(values = c("black","black"))+
  geom_signif(comparisons = list(c("black","barred")), annotations="F(1,105)=0.598, P=0.441",
              y_position = 2000, tip_length = 0, vjust=0.1)+
  scale_x_discrete(labels = c("barred","black"),name="Species")+
  ylab(label = "distance moved / body length")+
  theme_classic()


anova(M3) # take values for annotations
C=ggplot(df,aes(x=species,y=Velocity_cm.s,fill=species,col=species))+
  geom_violin(alpha=0.25, width=.9,color=NA)+
  geom_boxplot(width=0.25)+
  scale_fill_manual(values = c("burlywood","grey28"))+
  scale_color_manual(values = c("black","black"))+
  geom_signif(comparisons = list(c("black","barred")), annotations="F(1,105)=0.613, P=0.435",
              y_position = 23, tip_length = 0, vjust=0.1)+
  scale_x_discrete(labels = c("barred","black"),name="Species")+
  ylab(label = "velocity (cm/s)")+
  theme_classic()+
  theme(legend.position = "none")


anova(M4) # take values for annotations
D=ggplot(df,aes(x=species,y=Velocity_body_length.s,fill=species,col=species))+
  geom_violin(alpha=0.25, width=.9,color=NA)+
  geom_boxplot(width=0.25)+
  scale_fill_manual(values = c("burlywood","grey28"))+
  scale_color_manual(values = c("black","black"))+
  geom_signif(comparisons = list(c("black","barred")), annotations="F(1,105)=0.597, P=0.442",
              y_position = 3.8, tip_length = 0, vjust=0.1)+
  scale_x_discrete(labels = c("barred","black"),name="Species")+
  ylab(label = "velocity (body length/s)")+
  theme_classic()


plot_grid(A,B,C,D,nrow = 2,rel_widths = c(1,1.2),labels = c("A","B","C","D"))

