# Melanie J Heckwolf, 2024
# Analysis of consistancy between observers scoring the videos.

# Manuscript:
# Differences in colour pattern, behaviour and gene expression in the brain suggest divergent camouflage strategies in sympatric reef fish species
# MJ Heckwolf, J Gismann, M González-Santoro, F Coulmance, J Fuß, WO McMillan, O Puebla

# Result:
# there is no significant difference between observers in scoring the videos, we can thus distribute this task among us three:
# mean response duration: ANOVA: F (2, 33) = 0.18, P = 0.834; 
# mean distance swum: GLM (Poisson): Chi-square = 0.46, P = 0.795


# load packages:
library(ggplot2)
library(tidyverse)
library(car)

# set to your working directory:
setwd("C:/Users/melan/OneDrive/Dokumente/Results/Hamlets/05a_Behavior/01_analysis/escape_response_tank/observer_comparison/")


# load data
data <- read.csv2("escape_response_data_sheet_pilots_combined_dur_quad.csv",header=T)

# calculated speed:
data$speed <- data$duration/data$quadrants

# summarize trials into mean per individual:
data2 <- data %>% na.omit() %>% group_by(id,observer,species) %>% summarize(mean.duration=mean(duration),mean.quadrants=mean(quadrants), mean.speed=mean(speed))





####-------------------------------------------- visual data inspection:


# visual inspection of all data between observers:
ggplot(data,aes(x=observer,y=duration,col=species))+
  geom_boxplot()

ggplot(data,aes(x=observer,y=speed))+
  geom_boxplot()

ggplot(data,aes(fill=as.factor(quadrants),x=observer))+
  geom_bar(position = "stack")



# visual inspection of mean values between observers:
ggplot(data2,aes(x=observer,y=mean.duration,color=species))+
  geom_boxplot()

ggplot(data2,aes(x=observer,y=mean.speed))+
  geom_boxplot()

ggplot(data2,aes(fill=as.factor(quadrants),x=observer))+
  geom_bar(position = "stack")



####-------------------------------------------- data analysis:


###--------- duration:


# statistical analyses:
obs.dur <- aov(data2$mean.duration~data2$observer)
summary(obs.dur)

# model validation
par(mfrow = c(2, 2))
plot(obs.dur) # datapoint 7 could be an outlier.



# validated that removing datapoint 7 does not change the results:
obs.dur.no7 <- aov(data2$mean.duration[-7]~data2$observer[-7])
summary(obs.dur.no7)

# model validation
par(mfrow = c(2, 2))
plot(obs.dur.no7) #better fit. But 7 does not change the results.



# duration was log+1 transformed for analysis, see if that changes the outcome here:
obs.dur.log <- aov(log(data2$mean.duration+1)~data2$observer)
summary(obs.dur.log)

# model validation
par(mfrow = c(2, 2))
plot(obs.dur.log)




###--------- quadrants:


# statistical analyses:
obs.quad <- glm(round(mean.quadrants)~observer,data=data2,family = poisson)
Anova(obs.quad)

# model validation
par(mfrow = c(2, 1))
plot(residuals(obs.quad)~fitted(obs.quad))
boxplot(residuals(obs.quad)~data2$observer)


