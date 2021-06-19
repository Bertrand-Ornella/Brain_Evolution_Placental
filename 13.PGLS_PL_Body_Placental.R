
##########################################################################################
##############  Petrosal lobules - Body mass - PGLS, permutation tests, plots   ##########
##########################################################################################

#Before compiling, check tree and data to import, model to use for PGLS

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
library(car) #Levene's test
library(onewaytests) #Welsh test

setwd("~/Desktop/Placental_April_2021")

########################  Dataset and tree for Petrosal lobules  ##########################

#Import placental data to remove NAs from dataset and tree
dataPL<-read.csv("Placental_data_final.csv", header=T)
PL<-dataPL[c(1,12)] #select column for PL

PL.na<-na.omit(PL)
species<-PL.na$Species

#Import calibrated tree to take out NAs
tree<-read.nexus("Calibrated_placental1.trees")
tree_placental<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])

#Open dataset to match with the tree and save tree for analyses
data1<-read.csv("Placental_data_final.csv", header=T, row.names = 1)

PL.BM<-data1[c(1:5,6,11)] # body mass dataset
placental.data.BM<-na.omit(PL.BM)

#Check tree with data
name.check(tree_placental, placental.data.BM, data.names=NULL)

#Save dataset
write.csv(placental.data.BM,'Placental_data_final_PL_BM.csv')

#Save tree
write.nexus(tree_placental, file = "PL_BM_tree.trees")

##################################   Analyses PL   ######################################

#Open data and tree -body mass
placental.data<-read.csv("Placental_data_final_PL_BM.csv", header=T, row.names = 1)
tree_placental<-read.nexus("PL_BM_tree.trees")
plot(tree_placental,cex=.3)

####################################   PGLS for PL   ####################################

#Transform data to log10
placental.data$Body.mg<-log10(placental.data$Body_mass_mg)
placental.data$PL<-log10(placental.data$Petrosal_lobules_mm3)

#PL vs Body
Lambda_Body<-gls(PL ~ Body.mg, data=placental.data, 
                 correlation=corPagel(value=1,phy=tree_placental), method="ML")
Brownian_Body<-gls(PL ~ Body.mg, data=placental.data, 
                   correlation=corBrownian(1,phy=tree_placental), method="ML")
OU_Body<-gls(PL ~ Body.mg, data=placental.data, 
             correlation=corMartins(1,phy=tree_placental,fixed = TRUE), method="ML")
Blomberg_Body<-gls(PL ~ Body.mg, data=placental.data, 
                   correlation=corBlomberg(1.5,phy=tree_placental,fixed = TRUE), method="ML") #if 1, this is a Brownian model

anova(Lambda_Body,Brownian_Body,OU_Body, Blomberg_Body) # Lambda is the best
summary(Lambda_Body) #Final model

#### Now do the predicted and residuals
placental.data$Predicted_PL_Body <- predict(Lambda_Body)
placental.data$Residuals_PL_Body <- residuals(Lambda_Body)

#export dataframe
write.csv(placental.data,'Placental_data_final_PL_BM_res.csv')

################################ boxplots and permutation test ######################################

####Open data with PL residuals
placental.data<-read.csv("Placental_data_final_PL_BM_res.csv", header=T, row.names = 1)

##ggplot - boxplot - Residuals PL - Body
ggboxplot(placental.data,x="Group", y="Residuals_PL_Body", fill="Group", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  font("ylab", size = 7.8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='Residuals log10(Petrosal lobules vs Body mass)')

#Test normality and homogeneity of the variance
#Fisher-Pitman permutation test
sink("ANOVA_Welsh_log_PL_Body.txt")
shapiro.test(placental.data$Residuals_PL_Body) # NO normal distribution
leveneTest(Residuals_PL_Body ~ Group, data = placental.data) #NO there is homogeneity of variances
Welsh.test<-welch.test(Residuals_PL_Body~Group,data=placental.data)
paircomp(Welsh.test)
sink()

####### Density plots

#Create density plot
ggdensity(placental.data, x = "Residuals_PL_Body",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("blueviolet","goldenrod1","orangered1","turquoise"))

######################## ggplot - linear regression - Lambda ##########################

#open data for analyses # change rowname of X to Species
placental.data<-read.csv("Placental_data_final_PL_BM_res.csv", header=T)
names(placental.data)[names(placental.data) == "X"] <- "Species"

#Open tree
tree_placental<-read.nexus("PL_BM_tree.trees")

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

##### For PL vs. Body mass
#Create model PGLS regression line for each Group -- Body vs. PL
MeslineG_Br_B <-gls(PL ~ Body.mg, correlation=corPagel (1,phy=MesG_Tree), data=MesG)
PallineG_Br_B <-gls(PL ~ Body.mg, correlation=corPagel (1,phy=PalG_Tree), data=PalG)
EoclineS_Br_B <-gls(PL ~ Body.mg, correlation=corPagel (1,phy=EocS_Tree), data=EocS)
EoclineC_Br_B <-gls(PL ~ Body.mg, correlation=corPagel (1,phy=EocC_Tree), data=EocC)

#Prepare PGLS for each group
pgls.fit.MesG <- predict(MeslineG_Br_B) #predict values 
predframe.MesG <- with(MesG, data.frame(Species, Group, Body.mg, PL = pgls.fit.MesG))

pgls.fit.PalG <- predict(PallineG_Br_B) #predict values 
predframe.PalG <- with(PalG, data.frame(Species, Group, Body.mg, PL = pgls.fit.PalG))

pgls.fit.EocS <- predict(EoclineS_Br_B) #predict values 
predframe.EocS <- with(EocS, data.frame(Species, Group, Body.mg, PL = pgls.fit.EocS))

pgls.fit.EocC <- predict(EoclineC_Br_B) #predict values 
predframe.EocC <- with(EocC, data.frame(Species, Group, Body.mg, PL = pgls.fit.EocC))

#Make graph with PGLS corrected regressions - convex hulls -- Body
ggplot(placental.data, aes(Body.mg, PL, color = Group)) +
  geom_point(data = dplyr::filter(placental.data, Group == "Mesozoic"),
             size = 2, aes(color = "orangered1")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Paleocene"),
             size = 2, aes(color = "turquoise")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Eocene stem groups"),
             size = 2, aes(color = "goldenrod1")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Eocene crown groups"),
             size = 2, aes(color = "blueviolet")) +
  theme_minimal() + 
  scale_color_manual(name = "", values = c("blueviolet","goldenrod1","orangered1","turquoise"),labels = c("Eocene crown groups","Eocene stem groups","Mesozoic","Paleocene")) +
  #theme(legend.position = "top") +
  geom_convexhull(data = dplyr::filter(placental.data, Group == "Eocene crown groups"), color = "blueviolet", fill="blueviolet",alpha = 0.2) +
  geom_convexhull(data = dplyr::filter(placental.data, Group == "Eocene stem groups"), color = "goldenrod1", fill="goldenrod1",alpha = 0.2) +
  geom_convexhull(data = dplyr::filter(placental.data, Group == "Paleocene"), color = "turquoise", fill="turquoise",alpha = 0.2) +
  geom_convexhull(data = dplyr::filter(placental.data, Group == "Mesozoic"), color = "orangered1", fill="orangered1",alpha = 0.2) +
  labs(x = "log10(Body mass)", y = "log10(Petrosal lobule volume)") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,face = "bold")) 


#### END

