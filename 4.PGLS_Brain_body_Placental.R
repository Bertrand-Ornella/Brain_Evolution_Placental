
##########################################################################################
##################  Brain and Body size, PGLS, permutation tests, plots  #################
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

setwd("~/Desktop/Placental_April_2021")

#Open data
placental.data<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)

#Transform data to log10
placental.data$Brain_volume_cm3<-log10(placental.data$Brain_volume_cm3)
names(placental.data)[names(placental.data) == "Brain_volume_cm3"] <- "Brain"
placental.data$Body_mass_g<-log10(placental.data$Body_mass_g)
names(placental.data)[names(placental.data) == "Body_mass_g"] <- "Body"

########################### boxplots and permutation test - Brain ######################################

##ggplot - boxplot - Brain
ggboxplot(placental.data,x="Group", y="Brain", fill="Group", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='log10(Endocranial volume)')

#Test normality and homogeneity of the variance
#Fisher-Pitman permutation test
sink("ANOVA_Permutation_log_Brain.txt")
shapiro.test(placental.data$Brain) # YES normal distribution
bartlett.test(Brain ~ Group, data = placental.data) #YES there is homogeneity of variances
oneway_test(Brain~Group,data=placental.data)
placental.data$Group = factor(placental.data$Group, levels = c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups")) 
PT_Brain<-pairwisePermutationTest(Brain~Group,data=placental.data,method="fdr")
PT_Brain
sink()

########################### boxplots and permutation test - Body ######################################

##ggplot - boxplot - Body
ggboxplot(placental.data,x="Group", y="Body", fill="Group", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='log10(Body mass)')

#Test normality and homogeneity of the variance
#Fisher-Pitman permutation test
sink("ANOVA_Permutation_log_Body.txt")
shapiro.test(placental.data$Body) # YES normal distribution
bartlett.test(Body ~ Group, data = placental.data) #YES there is homogeneity of variances
oneway_test(Body~Group,data=placental.data)
placental.data$Group = factor(placental.data$Group, levels = c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups")) 
PT_Body<-pairwisePermutationTest(Body~Group,data=placental.data,method="fdr")
PT_Body
sink()

########################## Boxplot residuals- Brain vs Body ##############################

##ggplot - boxplot - Residuals brain vs. body
ggboxplot(placental.data,x="Group", y="residuals", fill="Group", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='Residuals (Endocranial volume vs. Body mass)')

#Fisher-Pitman permutation test
oneway_test(residuals~Group,data=placental.data)
placental.data$Group = factor(placental.data$Group, levels = c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups")) 
PT_Residuals<-pairwisePermutationTest(residuals~Group,data=placental.data,method="fdr")
PT_Residuals

######################## ggplot - linear regression - Lambda ##########################

#open data for analyses # change rowname of X to Species
placental.data<-read.csv("Placental_data_final_PEQ.csv", header=T)
names(placental.data)[names(placental.data) == "X"] <- "Species"

#Open tree
tree_placental<-read.nexus("PEQ_tree.trees")
plot(tree_placental,cex=.3)

#Transform data to log10
placental.data$Brain_volume_cm3<-log10(placental.data$Brain_volume_cm3)
names(placental.data)[names(placental.data) == "Brain_volume_cm3"] <- "Brain"
placental.data$Body_mass_g<-log10(placental.data$Body_mass_g)
names(placental.data)[names(placental.data) == "Body_mass_g"] <- "Body"

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

#Create model PGLS regression line for each Group
MeslineG_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=MesG_Tree), data=MesG)
PallineG_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=PalG_Tree), data=PalG)
EoclineS_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=EocS_Tree), data=EocS)
EoclineC_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=EocC_Tree), data=EocC)

#Prepare PGLS for each group -- Brain
pgls.fit.MesG <- predict(MeslineG_Br_B) #predict values for brain size
predframe.MesG <- with(MesG, data.frame(Species, Group, Body, Brain = pgls.fit.MesG))

pgls.fit.PalG <- predict(PallineG_Br_B) #predict values for brain size
predframe.PalG <- with(PalG, data.frame(Species, Group, Body, Brain = pgls.fit.PalG))

pgls.fit.EocS <- predict(EoclineS_Br_B) #predict values for brain size
predframe.EocS <- with(EocS, data.frame(Species, Group, Body, Brain = pgls.fit.EocS))

pgls.fit.EocC <- predict(EoclineC_Br_B) #predict values for brain size
predframe.EocC <- with(EocC, data.frame(Species, Group, Body, Brain = pgls.fit.EocC))

#Make graph with PGLS corrected regressions -- Brain/Body crown vs. stem
ggplot(placental.data, aes(Body, Brain, color = Group)) +
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
  labs(x = "log10(Body mass)", y = "log10(Endocranial volume)") +
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

MeslineG_Br_B
PallineG_Br_B
EoclineS_Br_B
EoclineC_Br_B

### END

