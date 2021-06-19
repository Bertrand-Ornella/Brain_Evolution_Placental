##########################################################################################
###############################   Boxplot PEQ, PGLS, diet   ##############################
##########################################################################################

library(ggpubr) #ggboxplot
library(coin) #one_way test
library(rcompanion) #pairwisePermutationTest
library(car) #Levene's test
library(onewaytests) #Welsh test

setwd("~/Desktop/Placental_April_2021")

#Open data with PEQ & residuals
data<-read.csv("Categories_boxplots_final.csv", header=T, row.names = 1)

placental.data<-subset(data, data$Group1 != "Euarchontoglires")

##ggplot - boxplot - PEQ
ggboxplot(placental.data,x="Group1", y="PEQ", fill="Group1", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Perissodactyla","Artiodactyla","Carnivoramorpha","Hyaenodonta"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='Phylogenetic encephalization quotient')

#Test normality and homogeneity of the variance
sink("ANOVA_permutation_PEQ_Diet_Order.txt")
shapiro.test(placental.data$PEQ) # not normal
leveneTest(PEQ ~ Group, data = placental.data) # yes homogeneity variances
oneway_test(PEQ~Group1,data=placental.data)
placental.data$Group1 = factor(placental.data$Group1, levels = c("Perissodactyla","Artiodactyla","Carnivoramorpha","Hyaenodonta")) 
PT_PEQ<-pairwisePermutationTest(PEQ~Group1,data=placental.data,method="fdr")
PT_PEQ
sink()

#########################

#Open data with PEQ & residuals
data<-read.csv("Categories_boxplots_final.csv", header=T, row.names = 1)

placental.data<-subset(data, data$Group2 != "Euarchontoglires")

##ggplot - boxplot - PEQ
ggboxplot(placental.data,x="Group2", y="PEQ", fill="Group2", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Paleocene Omnivorous-Carnivorous","Paleocene Herbivorous","Eocene Omnivorous-Carnivorous","Eocene Herbivorous"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 5)+
  rotate_x_text(90)+
  theme(legend.position = "right",legend.text = element_text(size = 5))+
  labs(x='Age', y='Phylogenetic encephalization quotient')

#Fisher-Pitman permutation test
sink("ANOVA_permutation_PEQ_Diet_Age.txt")
shapiro.test(placental.data$PEQ) # not normal
leveneTest(PEQ ~ Group, data = placental.data) # yes homogeneity variances
oneway_test(PEQ~Group2,data=placental.data)
placental.data$Group2 = factor(placental.data$Group2, levels = c("Paleocene Omnivorous-Carnivorous","Paleocene Herbivorous","Eocene Omnivorous-Carnivorous","Eocene Herbivorous")) 
PT_PEQ<-pairwisePermutationTest(PEQ~Group2,data=placental.data,method="fdr")
PT_PEQ
sink()

### END


