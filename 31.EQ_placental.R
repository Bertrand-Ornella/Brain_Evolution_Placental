
##########################################################################################
##########################   EQ, OLS, permutation tests, plots   #########################
##########################################################################################

library(ape) # read.nexus
library(geiger) # name.check
library(phytools) #open tree
library(ggplot2) #plots
library(ggConvexHull) #convex hull
library(nlme) # GLS analysis
library(coin) #one_way test
library(rcompanion) #pairwisePermutationTest
library(ggpubr) #ggboxplot
library(RRPP) #pairwise comparaisons
library(tibble) #make dataframe
library(ggpubr)#density plot
library(car) #Levene's test
library(onewaytests) #Welsh test

setwd("~/Desktop/Placental_April_2021")

##################################   Analyses EQ   ######################################

#Open data and tree
placental.data<-read.csv("Placental_data_final_Brain.csv", header=T, row.names = 1)

#### PGLS to determine EQ
#Transform data to log10
placental.data$Brain_volume_cm3<-log10(placental.data$Brain_volume_cm3)
names(placental.data)[names(placental.data) == "Brain_volume_cm3"] <- "Brain"
placental.data$Body_mass_g<-log10(placental.data$Body_mass_g)
names(placental.data)[names(placental.data) == "Body_mass_g"] <- "Body"

#Model using OLS (no phylogeny)
OLS<-gls(Brain ~ Body, data=placental.data, method="ML")
summary(OLS)

##################### Phylogenetic encephalozation quotient - PEQ ##################

#values obtained from OLS model above
a <- 0.6508165 #slope
b <--1.3256367 #intercept
Exp_b<-exp(b)
Exp_b

#Open dataset without log10
placental.data<-read.csv("Placental_data_final_Brain.csv", header=T, row.names = 1)

BM <-placental.data$Body_mass_g
Ec <- Exp_b*(BM)^a #Expected brain size
Ei <-placental.data$Brain_volume_cm3
EQ <-Ei/Ec

#Save EQ as column in dataset
EQ_df<-tibble(EQ)
colnames(EQ_df)<-c("EQ")
placental.data$EQ=EQ

#### Now do the predicted and residuals
placental.data$predicted <- predict(OLS)
placental.data$residuals <- residuals(OLS)

#export dataframe
write.csv(placental.data,'Placental_data_final_EQ.csv')

################################ boxplots and permutation test ######################################

#Open data with EQ & residuals
placental.data<-read.csv("Placental_data_final_EQ.csv", header=T, row.names = 1)

##ggplot - boxplot - EQ
ggboxplot(placental.data,x="Group", y="EQ", fill="Group", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='Encephalization quotient')

sink("ANOVA_Welsh_EQ.txt")
shapiro.test(placental.data$EQ)
leveneTest(EQ ~ Group, data = placental.data)
Welsh.test<-welch.test(EQ~Group,data=placental.data)
paircomp(Welsh.test)
sink()

################################### PEQ - density plot #################################

#Open data
placental.data<-read.csv("Placental_data_final_EQ.csv", header=T, row.names = 1)

#Create density plot
ggdensity(placental.data, x = "EQ",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("blueviolet","goldenrod1","orangered1","turquoise"))

#### END








