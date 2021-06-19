
##########################################################################################
#####################   Neocortex - PGLS, permutation tests, plots   #####################
##########################################################################################

#Before compiling, check tree and data to import, and model for PGLS

library(ape) # read.nexus
library(geiger) # name.check
library(phytools) #open tree
library(ggplot2) #plots
library(nlme) # GLS analysis
library(coin) #one_way test
library(rcompanion) #pairwisePermutationTest
library(ggpubr) #ggboxplot
library(RRPP) #pairwise comparaisons

setwd("~/Desktop/Placental_April_2021")

#Import placental data to remove NAs from dataset and tree
data<-read.csv("Placental_data_final.csv", header=T)
Neo<-data[c(1,13)] #select column with taxa for the new tree #Neocortex
Neo.na<-na.omit(Neo)
species<-Neo.na$Species

#Import calibrated tree to take out NAs
tree<-read.nexus("Calibrated_placental1.trees")
tree_placental<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])

#Open dataset to match with the tree and save tree for analyses
data1<-read.csv("Placental_data_final.csv", header=T, row.names = 1)
Neo1<-data1[c(1:5,12,13,16)]
placental.data<-na.omit(Neo1)

#Check tree with data
name.check(tree_placental, placental.data, data.names=NULL)

#Save dataset
write.csv(placental.data,'Placental_data_final_Neocortex.csv')

#Save tree
write.nexus(tree_placental, file = "Neocortex_tree.trees")

##################################### Analyses Neocortex #################################

#Open data and tree
placental.data<-read.csv("Placental_data_final_Neocortex.csv", header=T, row.names = 1)
tree_placental<-read.nexus("Neocortex_tree.trees")
plot(tree_placental,cex=.3)

#Transform data to log10
placental.data$Brain.surf<-log10(placental.data$Brain_surface_mm2)
placental.data$Neocortex<-log10(placental.data$Neocortex_surface_mm2)

#Neocortex vs Brain
Lambda_Neo<-gls(Neocortex ~ Brain.surf, data=placental.data, 
                  correlation=corPagel(value=0.9,phy=tree_placental), method="ML")
Brownian_Neo<-gls(Neocortex ~ Brain.surf, data=placental.data, 
                    correlation=corBrownian(1,phy=tree_placental), method="ML")
OU_Neo<-gls(Neocortex ~ Brain.surf, data=placental.data, 
              correlation=corMartins(1,phy=tree_placental,fixed = TRUE), method="ML")
Blomberg_Neo<-gls(Neocortex ~ Brain.surf, data=placental.data, 
                    correlation=corBlomberg(1.5,phy=tree_placental,fixed = TRUE), method="ML") #if 1, this is a Brownian model

anova(Lambda_Neo,Brownian_Neo,OU_Neo, Blomberg_Neo)
summary(Blomberg_Neo) #Final model

#### Now do the predicted and residuals
placental.data$Predicted_Neocortex <- predict(Blomberg_Neo)
placental.data$Residuals_Neocortex <- residuals(Blomberg_Neo)

#export dataframe
write.csv(placental.data,'Placental_data_final_Neocortex_res.csv')

################################ boxplots and permutation test ######################################

####Open data with Neocortex residuals
placental.data<-read.csv("Placental_data_final_Neocortex_res.csv", header=T, row.names = 1)

##ggplot - boxplot - Neocortex residuals
ggboxplot(placental.data,x="Group", y="Residuals_Neocortex", fill="Group", 
          palette=c("turquoise","goldenrod1","blueviolet"),
          order=c("Paleocene","Eocene stem groups","Eocene crown groups"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  font("ylab", size = 7.5)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='Residuals log10(Neocortex vs Endocranial surface area)')

#Test normality and homogeneity of the variance
#Fisher-Pitman permutation test
sink("ANOVA_Permutation_log_Neocortex.txt")
shapiro.test(placental.data$Residuals_Neocortex) # YES normal distribution
bartlett.test(Residuals_Neocortex ~ Group, data = placental.data) #YES there is homogeneity of variances
oneway_test(Residuals_Neocortex~Group,data=placental.data)
placental.data$Group = factor(placental.data$Group, levels = c("Paleocene","Eocene stem groups","Eocene crown groups")) 
PT_Residuals_Neocortex<-pairwisePermutationTest(Residuals_Neocortex~Group,data=placental.data,method="fdr")
PT_Residuals_Neocortex
sink()

######################## ggplot - linear regression - Lambda ##########################

#open data for analyses # change rowname of X to Species
placental.data<-read.csv("Placental_data_final_Neocortex_res.csv", header=T)
names(placental.data)[names(placental.data) == "X"] <- "Species"

#Open tree
tree_placental<-read.nexus("Neocortex_tree.trees")

#subset dataset per Group
PalG<- placental.data[ which(placental.data$Group=='Paleocene'), ]
EocS<-placental.data[ which(placental.data$Group=='Eocene stem groups'), ]
EocC<-placental.data[ which(placental.data$Group=='Eocene crown groups'), ]

#Find taxa with specific Group
PalG_sp<-subset(placental.data, Group=='Paleocene', select=Species)
EocS_sp<-subset(placental.data, Group=='Eocene stem groups', select=Species)
EocC_sp<-subset(placental.data, Group=='Eocene crown groups', select=Species)

#Select taxa per Group
PalG_species<- unlist (PalG_sp)
EocS_species<- unlist (EocS_sp)
EocC_species<- unlist (EocC_sp)

#subset tree per Group
PalG_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(PalG_species, tree_placental$tip.label)])
EocS_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(EocS_species, tree_placental$tip.label)])
EocC_Tree<-drop.tip(tree_placental,tree_placental$tip.label[-match(EocC_species, tree_placental$tip.label)])

#Create model PGLS regression line for each Group
PallineG_Br_B <-gls(Neocortex ~ Brain.surf, correlation=corBlomberg(value=1.5,phy=PalG_Tree,fixed = TRUE), data=PalG)
EoclineS_Br_B <-gls(Neocortex ~ Brain.surf, correlation=corBlomberg(value=1.5,phy=EocS_Tree,fixed = TRUE), data=EocS)
EoclineC_Br_B <-gls(Neocortex ~ Brain.surf, correlation=corBlomberg(value=1.5,phy=EocC_Tree,fixed = TRUE), data=EocC)

#Prepare PGLS for each group -- Brain
pgls.fit.PalG <- predict(PallineG_Br_B) #predict values for brain size
predframe.PalG <- with(PalG, data.frame(Species, Group, Brain.surf, Neocortex = pgls.fit.PalG))

pgls.fit.EocS <- predict(EoclineS_Br_B) #predict values for brain size
predframe.EocS <- with(EocS, data.frame(Species, Group, Brain.surf, Neocortex = pgls.fit.EocS))

pgls.fit.EocC <- predict(EoclineC_Br_B) #predict values for brain size
predframe.EocC <- with(EocC, data.frame(Species, Group, Brain.surf, Neocortex = pgls.fit.EocC))

#Make graph with PGLS corrected regressions -- Neocortex vs Brain
ggplot(placental.data, aes(Brain.surf, Neocortex, color = Group)) +
  geom_point(data = dplyr::filter(placental.data, Group == "Paleocene"),
             size = 2, aes(color = "turquoise")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Eocene stem groups"),
             size = 2, aes(color = "goldenrod1")) +
  geom_point(data = dplyr::filter(placental.data, Group == "Eocene crown groups"),
             size = 2, aes(color = "blueviolet")) +
  geom_line(data = dplyr::filter(predframe.PalG, Group == "Paleocene"), color = "turquoise",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.EocS, Group == "Eocene stem groups"), color = "goldenrod1",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.EocC, Group == "Eocene crown groups"), color = "blueviolet",
            linetype = 1.5) +
  theme_minimal() + 
  #theme(legend.position = "top") +
  scale_color_manual(name = "", values = c("blueviolet","goldenrod1","turquoise"),labels = c("Eocene crown groups","Eocene stem groups","Paleocene")) +
  labs(x = "log10(Endocranial surface area)", y = "log10(Neocortical surface area)") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,
            face = "bold"))
  #geom_text(data = dplyr::filter(placental.data, Group == "Paleocene"), color = "turquoise",
            #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1) +           
  #geom_text(data = dplyr::filter(placental.data, Group == "Eocene stem groups"), color = "goldenrod1",
            #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1) +
  #geom_text(data = dplyr::filter(placental.data, Group == "Eocene crown groups"), color = "blueviolet",
            #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1)


#### END


