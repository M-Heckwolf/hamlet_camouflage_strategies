# Melanie J Heckwolf, 2024
# Analysis of escape responses in Hypoplectrus spp reef fishes.

# Manuscript:
# Differences in colour pattern, behaviour and gene expression in the brain suggest divergent camouflage strategies in sympatric reef fish species
# MJ Heckwolf, J Gismann, M González-Santoro, F Coulmance, J Fuß, WO McMillan, O Puebla




# load packages:
library(ggplot2)
library(tidyverse)
library(lme4)
library(cowplot)
library(ggsignif)
library(car)
library(readxl)

# set to your working directory:
setwd("ADD/PATH/TO/YOUR/DATA/FILES/")

# load data:
data <- read.csv2("Supplementary_Table_S1.csv",header=T)

# remove scores where fish reacted to curtain but not hand:
data <- subset(data, reaction != "to_curtain_exclude") # this does not remove complete individuals, but some scores per individuals

# reduce dataframe to varaibles of interest, otherwise na.omit will kick out too much data:
data <- subset(data, select = c(id,escape_video_scored,species,location,duration,quadrants,speed,SDL_cm,escape_response_trial,starting_position,swimming_behaviour))
colnames(data) <- c("id","observer","species","location","duration","quadrants","speed","SDL_cm","escape_response_trial","starting_position","swimming_behaviour")

# summarize to mean response value per individual:
data2 <- data %>% na.omit() %>% group_by(id,observer,species,location) %>% summarize(mean.duration=mean(duration),mean.quadrants=mean(quadrants), mean.speed=mean(speed),SDL_cm=mean(SDL_cm))

# organize factor order:
data2$species <- factor(data2$species, levels=c("black","barred"))











### ------------------ escape response duration ---------------------------------------------------


# Linear mixed model with species as a fixed effect and observer as a random effect
mod.er <- lmer(data = data2, log(mean.duration+1) ~ species + (1 | observer))
Anova(mod.er) # LMM: Chi-squared (1) = 5.73, P = 0.017

# Figure 2A;
X=ggplot(data2,aes(x=species,y=mean.duration,fill=species,col=species))+
  geom_violin(alpha=0.25, width=.9, color = NA)+
  geom_boxplot(width=0.3)+
  #geom_jitter(width=.1)+
  scale_fill_manual(values = c("grey28","burlywood"))+
  scale_color_manual(values = c("black","black"))+
  geom_signif(comparisons = list(c("black","barred")), annotations="*",
              y_position = 2.1, tip_length = 0, vjust=0.1)+
  scale_x_discrete(labels = c("black","barred"),name="Species")+
  ylab(label = "escape duration (s)")+
  theme_classic()+
  theme(legend.position = "none")


# mean response values:
mean(data2$mean.duration[data2$species=="black"])
mean(data2$mean.duration[data2$species=="barred"])


# model validation:
par(mfrow = c(2, 2))
qqnorm(residuals(mod.er))
qqline(residuals(mod.er)) # looks good. Model is robust agains slight deviation from normality
plot(residuals(mod.er)~fitted(mod.er))
boxplot(residuals(mod.er)~data2$species) # homoscedasticity of residuals
boxplot(residuals(mod.er)~data2$observer)

# model validation tests (both should be >0.05)
leveneTest(residuals(mod.er)~data2$species) #0.09
shapiro.test(residuals(mod.er)) #0.521


# Result:
# Laboratory experiments on 50 barred and 47 black hamlets showed that the escape response duration was significantly longer in barred (~0.96s) compared to black hamlets (~0.80s) (Figure 2A; LMM: Chi-squared (1) = 5.73, P = 0.017). 



### --------------------------------- quadrants ------------------------------------------


# GLMM: Chi-squared (1) = 4.89, P = 0.027
mod.qds <- glmer(round(mean.quadrants)~species+(1|observer),data=data2,family = poisson, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
Anova(mod.qds)
summary(mod.qds)


#Figure 2B; 
ylim.prim <- c(0, 0.55)   # in this example, density
ylim.sec <- c(0, 30)    # in this example, count

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1] # there was a bug here

Y=ggplot(data2, aes(x=round(mean.quadrants),fill=species)) +
  geom_histogram(alpha=1, position="dodge", binwidth = 1)+
  scale_y_continuous("count") +
  scale_x_continuous("number of quadrants crossed", breaks = 0:5) +
  scale_fill_manual(values = c("grey28","burlywood"))+
  theme_classic()



# mean response values:
mean(round(data2$mean.quadrants[data2$species=="black"]))
mean(round(data2$mean.quadrants[data2$species=="barred"]))


# model validation:
par(mfrow = c(2, 1))
plot(residuals(mod.qds)~fitted(mod.qds))
boxplot(residuals(mod.qds)~data2$species)



# Result:
# Correspondingly, barred hamlets swam on average across three quadrants during their escape, which is significantly farther than the two quadrants on average in black hamlets (Figure 2B; GLMM: Chi-squared (1) = 4.89, P = 0.027). 




### --------------------------------- plot figure ------------------------------------------

# Plot Figure 2 A and B:
plot_grid(X,Y,nrow = 1,rel_widths = c(.5,1),labels = c("A","B"))








### ---------------------------- test for confounding factors ----------------------------------------------------


#### visually explore other factors that could have introduced a bias:

#------------------------observer effect:

ggplot(data,aes(x=observer,y=duration))+
  geom_boxplot()

ggplot(data,aes(x=observer,y=speed))+
  geom_boxplot()

ggplot(data,aes(fill=as.factor(quadrants),x=observer))+
  geom_bar(position = "stack")

#------------------------location effect:

ggplot(data,aes(x=location,y=duration))+
  geom_boxplot()

ggplot(data,aes(x=location,y=speed))+
  geom_boxplot()

ggplot(data,aes(fill=as.factor(quadrants),x=location))+
  geom_bar(position = "stack")

#-------------------------trial effect:

ggplot(data,aes(x=as.factor(escape_response_trial),y=duration))+
  geom_boxplot()

ggplot(data,aes(x=as.factor(escape_response_trial),y=speed))+
  geom_boxplot()

ggplot(data,aes(fill=as.factor(quadrants),x=as.factor(escape_response_trial)))+
  geom_bar(position = "stack")

#-------------------------starting_pos:

ggplot(data,aes(x=starting_position,y=duration))+
  geom_boxplot()

ggplot(data,aes(x=starting_position,y=speed))+
  geom_boxplot()

ggplot(data,aes(fill=as.factor(quadrants),x=starting_position))+
  geom_bar(position = "stack")


#-------------------------swimming_behaviour:

ggplot(data,aes(x=swimming_behaviour,y=duration))+
  geom_boxplot()

ggplot(data,aes(x=swimming_behaviour,y=speed))+
  geom_boxplot()

ggplot(data,aes(fill=as.factor(quadrants),x=swimming_behaviour))+
  geom_bar(position = "stack")


# factors dont need further consideration.







### ---------------------------- could size differences explain the differences in escape responses? ----------------------------------------------------


#-------------------------- size effect:
ggplot(data2,aes(x=species,y=as.numeric(SDL_cm)))+
  geom_boxplot()

mod.size <- aov(data2$SDL_cm~data2$species)

summary(mod.size)

mean(na.omit(data2$SDL_cm[data2$species=="black"]))
mean(na.omit(data2$SDL_cm[data2$species=="barred"]))


# model follows assumptions for homoscedasticity and normality of residuals:
par(mfrow = c(2, 2))
plot(mod.size)

# model validation tests (both should be >0.05)
leveneTest(residuals(mod.size)~data2$species) #0.564
shapiro.test(residuals(mod.size)) #0.698


# See Discussion: However, the black and barred hamlets from Bocas we used within this study did not differ in size and are representative of the size range of these populations (ANOVA; F 1,95 = 0.118; P=0.732; mean SDL (cm): blacks: 7.34, barred: 7.30). 



