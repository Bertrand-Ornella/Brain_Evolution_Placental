
##########################################################################################
###############  Olfactory bulbs - Body - PGLS, permutation tests, plots   ###############
##########################################################################################

#Before compiling, check tree and data to import, model to use for PGLS

library(ape) # read.nexus
library(geiger) # name.check
library(phytools) #open tree
library(ggplot2) #plots
library(nlme) # GLS analysis
library(coin) #one_way test
library(rcompanion) #pairwisePermutationTest
library(ggpubr) #ggboxplot
library(RRPP) #pairwise comparaisons
library(car) #Levene's test

setwd("~/Desktop/Placental_April_2021")

########################  Dataset and tree for Olfactory bulbs  ##########################

#Import placental data to remove NAs from dataset and tree
dataOB<-read.csv("Placental_data_final.csv", header=T)
OB<-dataOB[c(1,11)] #select column with taxa for the new tree
OB.na<-na.omit(OB)
species<-OB.na$Species

#Import calibrated tree to take out NAs
tree<-read.nexus("Calibrated_placental1.trees")
tree_placental<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])

#Open dataset to match with the tree and save tree for analyses
data1<-read.csv("Placental_data_final.csv", header=T, row.names = 1)
OB1<-data1[c(1:5,6,9,10,14)]
placental.data<-na.omit(OB1)

#Check tree with data
name.check(tree_placental, placental.data, data.names=NULL)

#Save dataset
write.csv(placental.data,'Placental_data_final_OB.csv')

#Save tree
write.nexus(tree_placental, file = "OB_tree.trees")

##################################   Analyses OB   ######################################

#Open data and tree
placental.data<-read.csv("Placental_data_final_OB.csv", header=T, row.names = 1)
tree_placental<-read.nexus("OB_tree.trees")
par(mfrow=c(1,1))
plot(tree_placental,cex=.3)

####################################   PGLS for OB   ####################################

#Transform data to log10
placental.data$Body.mg<-log10(placental.data$Body_mass_mg)
placental.data$Brain.mg<-log10(placental.data$Brain_volume_mm3)
placental.data$OB<-log10(placental.data$Olfactory_bulbs_mm3)

#OB vs Brain
Lambda_Brain<-gls(OB ~ Brain.mg, data=placental.data, 
            correlation=corPagel(value=0.9,phy=tree_placental), method="ML")
Brownian_Brain<-gls(OB ~ Brain.mg, data=placental.data, 
              correlation=corBrownian(1,phy=tree_placental), method="ML")
OU_Brain<-gls(OB ~ Brain.mg, data=placental.data, 
        correlation=corMartins(1,phy=tree_placental,fixed = TRUE), method="ML")
Blomberg_Brain<-gls(OB ~ Brain.mg, data=placental.data, 
              correlation=corBlomberg(1.5,phy=tree_placental,fixed = TRUE), method="ML") #if 1, this is a Brownian model

anova(Lambda_Brain,Brownian_Brain,OU_Brain, Blomberg_Brain) # Lambda is the best
summary(Lambda_Brain) #Final model

#### Now do the predicted and residuals
placental.data$Predicted_OB_Brain <- predict(Lambda_Brain)
placental.data$Residuals_OB_Brain <- residuals(Lambda_Brain)

#OB vs Body
Lambda_Body<-gls(OB ~ Body.mg, data=placental.data, 
            correlation=corPagel(value=0.9,phy=tree_placental), method="ML")
Brownian_Body<-gls(OB ~ Body.mg, data=placental.data, 
              correlation=corBrownian(1,phy=tree_placental), method="ML")
OU_Body<-gls(OB ~ Body.mg, data=placental.data, 
        correlation=corMartins(1,phy=tree_placental,fixed = TRUE), method="ML")
Blomberg_Body<-gls(OB ~ Body.mg, data=placental.data, 
              correlation=corBlomberg(1.5,phy=tree_placental,fixed = TRUE), method="ML") #if 1, this is a Brownian model

anova(Lambda_Body,Brownian_Body,OU_Body, Blomberg_Body) # Lambda is the best
summary(Lambda_Body) #Final model

#### Now do the predicted and residuals
placental.data$Predicted_OB_Body <- predict(Lambda_Body)
placental.data$Residuals_OB_Body <- residuals(Lambda_Body)

#export dataframe
write.csv(placental.data,'Placental_data_final_OB_res.csv')

################################ boxplots and permutation test ######################################

####Open data with OB residuals
placental.data<-read.csv("Placental_data_final_OB_res.csv", header=T, row.names = 1)

##ggplot - boxplot - Residuals OB - Brain
ggboxplot(placental.data,x="Group", y="Residuals_OB_Brain", fill="Group", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  font("ylab", size = 7.8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='Residuals log10(Olfactory bulb vs Endocranial volume)')

#Test normality and homogeneity of the variance
#Fisher-Pitman permutation test
sink("ANOVA_Permutation_log_OB_Brain.txt")
shapiro.test(placental.data$Residuals_OB_Brain) # NO normal distribution
leveneTest(Residuals_OB_Brain ~ Group, data = placental.data) #YES there is homogeneity of variances
oneway_test(Residuals_OB_Brain~Group,data=placental.data)
placental.data$Group = factor(placental.data$Group, levels = c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups")) 
PT_Residuals_OB_Brain<-pairwisePermutationTest(Residuals_OB_Brain~Group,data=placental.data,method="fdr")
PT_Residuals_OB_Brain
sink()

##ggplot - boxplot - Residuals OB - Body
ggboxplot(placental.data,x="Group", y="Residuals_OB_Body", fill="Group", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  font("ylab", size = 7.8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='Residuals log10(Olfactory bulbs vs Body mass)')

#Test normality and homogeneity of the variance
#Fisher-Pitman permutation test
sink("ANOVA_Permutation_log_OB_Body.txt")
shapiro.test(placental.data$Residuals_OB_Body) # YES normal distribution
bartlett.test(Residuals_OB_Body ~ Group, data = placental.data) #YES there is homogeneity of variances
oneway_test(Residuals_OB_Body~Group,data=placental.data)
placental.data$Group = factor(placental.data$Group, levels = c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups")) 
PT_Residuals_OB_Body<-pairwisePermutationTest(Residuals_OB_Body~Group,data=placental.data,method="fdr")
PT_Residuals_OB_Body
sink()

####### Density plots

#Create density plot
ggdensity(placental.data, x = "Residuals_OB_Brain",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("blueviolet","goldenrod1","orangered1","turquoise"))

#Create density plot
ggdensity(placental.data, x = "Residuals_OB_Body",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("blueviolet","goldenrod1","orangered1","turquoise"))

#################### ggplot - linear regression OB - Brain - Lambda ######################

#open data for analyses # change rowname of X to Species
placental.data<-read.csv("Placental_data_final_OB_res.csv", header=T)
names(placental.data)[names(placental.data) == "X"] <- "Species"

#Open tree
tree_placental<-read.nexus("OB_tree.trees")

#subset dataset per Group
MesG<- placental.data[ which(placental.data$Group=='Mesozoic'), ]
PalG<- placental.data[ which(placental.data$Group=='Paleocene'), ]
EocS<-placental.data[ which(placental.data$Group=='Eocene stem groups'), ]
EocC<-placental.data[ which(placental.data$Group=='Eocene crown groups'), ]

#Find taxa with specific Group
MesG_sp<-subset(placental.data, Group=='Mesozoic', select=Species)
PalG_sp<-subset(placental.data, Group=='Paleocene', select=Species)
EocS_sp<-subset(placental.data, Group=='Eocene stem groups', select=Species)
EocC_sp<-subset(placental.data, Group=='Eocene crown groups', select=Species)

#Select taxa per Group
MesG_species<- unlist (MesG_sp)
PalG_species<- unlist (PalG_sp)
EocS_species<- unlist (EocS_sp)
EocC_species<- unlist (EocC_sp)

#subset tree per Group
MesG_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(MesG_species, tree_placental$tip.label)])
PalG_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(PalG_species, tree_placental$tip.label)])
EocS_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(EocS_species, tree_placental$tip.label)])
EocC_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(EocC_species, tree_placental$tip.label)])

###### For OB vs. Brain size
#Create model PGLS regression line for each Group
MeslineG_Br_B <-gls(OB ~ Brain.mg, correlation=corPagel (0.9,phy=MesG_Tree), data=MesG)
PallineG_Br_B <-gls(OB ~ Brain.mg, correlation=corPagel (0.9,phy=PalG_Tree), data=PalG)
EoclineS_Br_B <-gls(OB ~ Brain.mg, correlation=corPagel (0.9,phy=EocS_Tree), data=EocS)
EoclineC_Br_B <-gls(OB ~ Brain.mg, correlation=corPagel (0.9,phy=EocC_Tree), data=EocC)

#Prepare PGLS for each group
pgls.fit.MesG <- predict(MeslineG_Br_B) 
predframe.MesG <- with(MesG, data.frame(Species, Group, Brain.mg, OB = pgls.fit.MesG))

pgls.fit.PalG <- predict(PallineG_Br_B) 
predframe.PalG <- with(PalG, data.frame(Species, Group, Brain.mg, OB = pgls.fit.PalG))

pgls.fit.EocS <- predict(EoclineS_Br_B) 
predframe.EocS <- with(EocS, data.frame(Species, Group, Brain.mg, OB = pgls.fit.EocS))

pgls.fit.EocC <- predict(EoclineC_Br_B)
predframe.EocC <- with(EocC, data.frame(Species, Group, Brain.mg, OB = pgls.fit.EocC))

#Make graph with PGLS corrected regressions -- OB vs Brain
ggplot(placental.data, aes(Brain.mg, OB, color = Group)) +
  geom_point(data = dplyr::filter(placental.data, Group == "Mesozoic"),
             size = 2, aes(color = "orangered1")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Paleocene"),
             size = 2, aes(color = "turquoise")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Eocene stem groups"),
             size = 2, aes(color = "goldenrod1")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Eocene crown groups"),
             size = 2, aes(color = "blueviolet")) +
  geom_line(data = dplyr::filter(predframe.MesG, Group == "Mesozoic"), color = "orangered1",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.PalG, Group == "Paleocene"), color = "turquoise",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.EocS, Group == "Eocene stem groups"), color = "goldenrod1",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.EocC, Group == "Eocene crown groups"), color = "blueviolet",
            linetype = 1.5) +
  theme_minimal() + 
  #theme(legend.position = "top") +
  scale_color_manual(name = "", values = c("blueviolet","goldenrod1","orangered1","turquoise"),
            labels = c("Eocene crown groups","Eocene stem groups","Mesozoic","Paleocene")) +
  labs(x = "log10(Endocranial volume)", y = "log10(Olfactory bulb volume)") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,
            face = "bold"))

#geom_text(data = dplyr::filter(placental.data, Group == "Mesozoic"), color = "orangered1",
          #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1)+ 
  #geom_text(data = dplyr::filter(placental.data, Group == "Paleocene"), color = "turquoise",
           #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1) +           
  #geom_text(data = dplyr::filter(placental.data, Group == "Eocene stem groups"), color = "goldenrod1",
            #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1) +
  #geom_text(data = dplyr::filter(placental.data, Group == "Eocene crown groups"), color = "blueviolet",
            #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1)

#################### ggplot - linear regression OB - Body - Lambda ######################

#open data for analyses # change rowname of X to Species
placental.data<-read.csv("Placental_data_final_OB_res.csv", header=T)
names(placental.data)[names(placental.data) == "X"] <- "Species"

#Open tree
tree_placental<-read.nexus("OB_tree.trees")

#subset dataset per Group
MesG<- placental.data[ which(placental.data$Group=='Mesozoic'), ]
PalG<- placental.data[ which(placental.data$Group=='Paleocene'), ]
EocS<-placental.data[ which(placental.data$Group=='Eocene stem groups'), ]
EocC<-placental.data[ which(placental.data$Group=='Eocene crown groups'), ]

#Find taxa with specific Group
MesG_sp<-subset(placental.data, Group=='Mesozoic', select=Species)
PalG_sp<-subset(placental.data, Group=='Paleocene', select=Species)
EocS_sp<-subset(placental.data, Group=='Eocene stem groups', select=Species)
EocC_sp<-subset(placental.data, Group=='Eocene crown groups', select=Species)

#Select taxa per Group
MesG_species<- unlist (MesG_sp)
PalG_species<- unlist (PalG_sp)
EocS_species<- unlist (EocS_sp)
EocC_species<- unlist (EocC_sp)

#subset tree per Group
MesG_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(MesG_species, tree_placental$tip.label)])
PalG_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(PalG_species, tree_placental$tip.label)])
EocS_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(EocS_species, tree_placental$tip.label)])
EocC_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(EocC_species, tree_placental$tip.label)])

###### For OB vs. Body mass
#Create model PGLS regression line for each Group
MeslineG_Br_B <-gls(OB ~ Body.mg, correlation=corPagel (0.9,phy=MesG_Tree), data=MesG)
PallineG_Br_B <-gls(OB ~ Body.mg, correlation=corPagel (0.9,phy=PalG_Tree), data=PalG)
EoclineS_Br_B <-gls(OB ~ Body.mg, correlation=corPagel (0.9,phy=EocS_Tree), data=EocS)
EoclineC_Br_B <-gls(OB ~ Body.mg, correlation=corPagel (0.9,phy=EocC_Tree), data=EocC)

#Prepare PGLS for each group
pgls.fit.MesG <- predict(MeslineG_Br_B) 
predframe.MesG <- with(MesG, data.frame(Species, Group, Body.mg, OB = pgls.fit.MesG))

pgls.fit.PalG <- predict(PallineG_Br_B) 
predframe.PalG <- with(PalG, data.frame(Species, Group, Body.mg, OB = pgls.fit.PalG))

pgls.fit.EocS <- predict(EoclineS_Br_B) 
predframe.EocS <- with(EocS, data.frame(Species, Group, Body.mg, OB = pgls.fit.EocS))

pgls.fit.EocC <- predict(EoclineC_Br_B)
predframe.EocC <- with(EocC, data.frame(Species, Group, Body.mg, OB = pgls.fit.EocC))

#Make graph with PGLS corrected regressions -- OB/Body
ggplot(placental.data, aes(Body.mg, OB, color = Group)) +
  geom_point(data = dplyr::filter(placental.data, Group == "Mesozoic"),
               size = 2, aes(color = "orangered1")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Paleocene"),
               size = 2, aes(color = "turquoise")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Eocene stem groups"),
               size = 2, aes(color = "goldenrod1")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Eocene crown groups"),
               size = 2, aes(color = "blueviolet")) +
  geom_line(data = dplyr::filter(predframe.MesG, Group == "Mesozoic"), color = "orangered1",
              linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.PalG, Group == "Paleocene"), color = "turquoise",
              linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.EocS, Group == "Eocene stem groups"), color = "goldenrod1",
              linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.EocC, Group == "Eocene crown groups"), color = "blueviolet",
              linetype = 1.5) +
  theme_minimal() + 
  #theme(legend.position = "top") +
  scale_color_manual(name = "", values = c("blueviolet","goldenrod1","orangered1","turquoise"),
              labels = c("Eocene crown groups","Eocene stem groups","Mesozoic","Paleocene")) +
  labs(x = "log10(Body mass)", y = "log10(Olfactory bulb volume)") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,
              face = "bold"))
  
#geom_text(data = dplyr::filter(placental.data, Group == "Mesozoic"), color = "orangered1",
            #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1)+ 
  #geom_text(data = dplyr::filter(placental.data, Group == "Paleocene"), color = "turquoise",
              #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1) +           
  #geom_text(data = dplyr::filter(placental.data, Group == "Eocene stem groups"), color = "goldenrod1",
              #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1) +
  #geom_text(data = dplyr::filter(placental.data, Group == "Eocene crown groups"), color = "blueviolet",
              #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1)
  
#### END
