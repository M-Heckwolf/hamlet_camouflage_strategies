# Melanie J Heckwolf, 2024
# Transcriptomics - general read and sequencing statistics

# Manuscript:
# Differences in colour pattern, behaviour and gene expression in the brain suggest divergent camouflage strategies in sympatric reef fish species
# MJ Heckwolf, J Gismann, M González-Santoro, F Coulmance, J Fuß, WO McMillan, O Puebla


# load libraries
library(car)


# set your working directory:
setwd("PATH/TO/YOUR/FILES/")

# load data
htseq <- read.csv2("Transcriptomics_summary_stats.csv",header = T)


# plot data
X=ggplot(htseq, aes(x=species,y=tot.counts,fill=species,col=species))+
  geom_boxplot(width=0.45)+
  geom_jitter(width=.1)+
  scale_fill_manual(values = c("grey28","burlywood"))+
  scale_color_manual(values = c("grey","black"))+
  ylab(label = "total number of reads sequenced")+
  scale_x_discrete(labels = c("black","barred"),name="Species")+
  theme_classic()+
  theme(legend.position = "none")

Y=ggplot(htseq, aes(x=species,y=aligned_not_ambiguous,fill=species,col=species))+
  geom_boxplot(width=0.45)+
  geom_jitter(width=.1)+
  scale_fill_manual(values = c("grey28","burlywood"),
                    labels = c("black", "barred"))+
  scale_color_manual(values = c("grey","black"),
                     labels = c("black", "barred"))+
  ylab(label = "number of filtered and mapped reads")+
  scale_x_discrete(labels = c("black","barred"),name="Species")+
  theme_classic()

plot_grid(X,Y,nrow = 1,rel_widths = c(1,1.3),labels = c("A","B"))


#------------------------------------------- test statistics

#### total read counts
mod.total <- aov(htseq$tot.counts~htseq$species)
summary(mod.total)

# model validation:
par(mfrow = c(2, 2))
plot(mod.total)

# model validation tests (both should be >0.05)
leveneTest(residuals(mod.total)~htseq$species) #0.308
shapiro.test(residuals(mod.total)) #0.395



#### uniquely mapped reads:
mod.um <- aov(htseq$aligned_not_ambiguous~htseq$species)
summary(mod.um)

# model validation:
par(mfrow = c(2, 2))
plot(mod.um)

# model validation tests (both should be >0.05)
leveneTest(residuals(mod.um)~htseq$species) #0.307
shapiro.test(residuals(mod.um)) #0.455



#### comparing the proportions of uniquely mapped reads between species:
ggplot(htseq, aes(x=species,y=aligned_not_ambiguous/tot.counts,fill=species,col=species))+
  geom_boxplot(width=0.45)+
  geom_jitter(width=.1)+
  scale_fill_manual(values = c("grey28","burlywood"),
                    labels = c("black", "barred"))+
  scale_color_manual(values = c("grey","black"),
                     labels = c("black", "barred"))+
  ylab(label = "proportion of filtered and mapped reads")+
  scale_x_discrete(labels = c("black","barred"),name="Species")+
  theme_classic()


# extract values of interest:
barred <- htseq[htseq$species=="Hpue",]
black <- htseq[htseq$species=="Hnig",]

summary(barred$aligned_not_ambiguous); sd(barred$aligned_not_ambiguous)
summary(barred$tot.counts); sd(barred$tot.counts)
summary(barred$aligned_not_ambiguous/barred$tot.counts); sd(barred$aligned_not_ambiguous/barred$tot.counts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.4525  0.4610  0.4635  0.4659  0.4717  0.4793 

summary(black$aligned_not_ambiguous); sd(black$aligned_not_ambiguous)
summary(black$tot.counts); sd(black$tot.counts)
summary(black$aligned_not_ambiguous/black$tot.counts);sd(black$aligned_not_ambiguous/black$tot.counts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.4550  0.4611  0.4727  0.4720  0.4745  0.5272 
