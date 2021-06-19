
##########################################################################################
##########################   PEQ, PGLS, permutation tests, plots   #######################
##########################################################################################

#Before compiling, check tree and data to import, the slope (a) and intercept (b)
#model to use for PGLS

library(ape) # read.nexus
library(phytools) #open tree
library(ggplot2) #plots
library(nlme) # GLS analysis
library(coin) #one_way test
library(rcompanion) #pairwisePermutationTest
library(ggpubr) #ggboxplot
library(RRPP) #pairwise comparaisons
library(tibble) #make dataframe
library(ggpubr)#density plot
library(geiger)#name.check
library(car) #Levene's test
library(onewaytests) #Welsh test

setwd("~/Desktop/For github")

##############################  Dataset and tree for PEQ  ###############################

#Import placental data to remove NAs from dataset and tree
data<-read.csv("Placental_data_final.csv", header=T)
Endo<-data[c(1,9)] #select column with taxa for the new tree #Brain size (same for PEQ)
Endo.na<-na.omit(Endo)
species<-Endo.na$Species

#Import calibrated tree to take out NAs
tree<-read.nexus("Calibrated_placental1.trees")
par(mfrow=c(1,1))
plot(tree,cex=0.3)
tree_placental<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])

#Open dataset to match with the tree and save tree for analyses
data1<-read.csv("Placental_data_final.csv", header=T, row.names = 1)
Endo1<-data1[c(1:5,7,8)]
placental.data<-na.omit(Endo1)

#Check tree with data
name.check(tree_placental, placental.data, data.names=NULL)

#Brain vs. Body mass dataset
write.csv(placental.data,'Placental_data_final_Brain.csv')

#Tree Brain and PEQ
write.nexus(tree_placental, file = "PEQ_tree.trees")

##################################   Analyses PEQ   ######################################

#Open data and tree
placental.data<-read.csv("Placental_data_final_Brain.csv", header=T, row.names = 1)
tree_placental<-read.nexus("PEQ_tree.trees")
plot(tree_placental,cex=0.3)

#### PGLS to determine PEQ
#Transform data to log10
placental.data$Brain_volume_cm3<-log10(placental.data$Brain_volume_cm3)
names(placental.data)[names(placental.data) == "Brain_volume_cm3"] <- "Brain"
placental.data$Body_mass_g<-log10(placental.data$Body_mass_g)
names(placental.data)[names(placental.data) == "Body_mass_g"] <- "Body"

Lambda<-gls(Brain ~ Body, data=placental.data, 
            correlation=corPagel(value=1,phy=tree_placental), method="ML")
Brownian<-gls(Brain ~ Body, data=placental.data, 
              correlation=corBrownian(1,phy=tree_placental), method="ML")
OU<-gls(Brain ~ Body, data=placental.data, 
        correlation=corMartins(1,phy=tree_placental,fixed = TRUE), method="ML")
Blomberg<-gls(Brain ~ Body, data=placental.data, 
              correlation=corBlomberg(1.5,phy=tree_placental,fixed = TRUE), method="ML") #if 1, this is a Brownian model

anova(Lambda,Brownian,OU, Blomberg) # Lambda is the best
summary(Lambda) #Final model

##################### Phylogenetic encephalozation quotient - PEQ ##################

#values obtained from Lambda model above
a <- 0.6424423 #slope
b <--1.5807302 #intercept
Exp_b<-exp(b)
Exp_b

#Open dataset without log10
placental.data<-read.csv("Placental_data_final_Brain.csv", header=T, row.names = 1)

BM <-placental.data$Body_mass_g
Ec <- Exp_b*(BM)^a #Expected brain size
Ei <-placental.data$Brain_volume_cm3
PEQ <-Ei/Ec

#Save PEQ as column in dataset
PEQ_df<-tibble(PEQ)
colnames(PEQ_df)<-c("PEQ")
placental.data$PEQ=PEQ

#### Now do the predicted and residuals
placental.data$predicted <- predict(Lambda)
placental.data$residuals <- residuals(Lambda)

#export dataframe
write.csv(placental.data,'Placental_data_final_PEQ.csv')

################################ boxplots and permutation test ######################################

#Open data with PEQ & residuals
placental.data<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)

##ggplot - boxplot - PEQ
ggboxplot(placental.data,x="Group", y="PEQ", fill="Group", 
  palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
  order=c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups"),
  xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='Phylogenetic encephalization quotient')

#Test normality and homogeneity of the variance
#Not normal and variances not equal - Welch test ANOVA
sink("ANOVA_Welsh_PEQ.txt")
shapiro.test(placental.data$PEQ)
leveneTest(PEQ ~ Group, data = placental.data)
Welsh.test<-welch.test(PEQ~Group,data=placental.data)
paircomp(Welsh.test)
sink()

#Transform PEQ to log PEQ for test
placental.data$PEQ<-log10(placental.data$PEQ)
names(placental.data)[names(placental.data) == "PEQ"] <- "log_PEQ"

#Fisher-Pitman permutation test #almost not normal so prefer to use leveneTest for log_PEQ
sink("ANOVA_Permutation_log_PEQ.txt")
shapiro.test(placental.data$log_PEQ)
leveneTest(log_PEQ ~ Group, data = placental.data)
oneway_test(log_PEQ~Group,data=placental.data)
placental.data$Group = factor(placental.data$Group, levels = c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups")) 
PT_log_PEQ<-pairwisePermutationTest(log_PEQ~Group,data=placental.data,method="fdr")
PT_log_PEQ
sink()

##ggplot - boxplot - log 10 PEQ
ggboxplot(placental.data,x="Group", y="log_PEQ", fill="Group", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='Log10 (Phylogenetic encephalization quotient)')

################################### PEQ - density plot #################################

#Open data
placental.data<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)

#Create density plot
ggdensity(placental.data, x = "PEQ",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("blueviolet","goldenrod1","orangered1","turquoise"))

### Brain and body size

#Transform data to log10
placental.data$Brain_volume_cm3<-log10(placental.data$Brain_volume_cm3)
names(placental.data)[names(placental.data) == "Brain_volume_cm3"] <- "Brain"
placental.data$Body_mass_g<-log10(placental.data$Body_mass_g)
names(placental.data)[names(placental.data) == "Body_mass_g"] <- "Body"

#Create density plot
ggdensity(placental.data, x = "Brain",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("blueviolet","goldenrod1","orangered1","turquoise"))

#Create density plot
ggdensity(placental.data, x = "Body",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("blueviolet","goldenrod1","orangered1","turquoise"))

#Open data
placental.data<-read.csv("Placental_data_final.csv", header=T, row.names = 1)

placental.data$Body_mass_g<-log10(placental.data$Body_mass_g)
names(placental.data)[names(placental.data) == "Body_mass_g"] <- "Body"

#Create density plot
ggdensity(placental.data, x = "Body",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("blueviolet","goldenrod1","orangered1","turquoise"))

