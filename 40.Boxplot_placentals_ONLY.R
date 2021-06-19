
##########################################################################################
##########################   PEQ, PGLS, permutation tests, plots   #######################
##########################################################################################

library(coin) #one_way test
library(rcompanion) #pairwisePermutationTest
library(ggpubr) #ggboxplot
library(car) #Levene's test
library(onewaytests) #Welsh test

setwd("~/Desktop/Placental_April_2021")

################################ boxplots and permutation test ######################################

#Open data with PEQ & residuals
data<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)
placental.data<-data[data$Abbreviations != "Pu_an" & data$Abbreviations != "Pt_mo", ]

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
sink("ANOVA_Welsh_PEQ_placental_ONLY.txt")
shapiro.test(placental.data$PEQ)
leveneTest(PEQ ~ Group, data = placental.data)
Welsh.test<-welch.test(PEQ~Group,data=placental.data)
paircomp(Welsh.test)
sink()

#Transform PEQ to log PEQ for test
placental.data$PEQ<-log10(placental.data$PEQ)
names(placental.data)[names(placental.data) == "PEQ"] <- "log_PEQ"

#Fisher-Pitman permutation test #almost not normal so prefer to use leveneTest for log_PEQ
sink("ANOVA_Permutation_log_PEQ_Placental_only_topo1.txt")
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

########

####Open data with OB residuals
data<-read.csv("Placental_data_final_OB_res.csv", header=T, row.names = 1)
placental.data<-data[data$Abbreviations != "Pu_an" & data$Abbreviations != "Pt_mo", ]

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
sink("ANOVA_Permutation_log_OB_Brain_placental_only_topo1.txt")
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
sink("ANOVA_Permutation_log_OB_Body_placental_only_topo1.txt")
shapiro.test(placental.data$Residuals_OB_Body) # YES normal distribution
bartlett.test(Residuals_OB_Body ~ Group, data = placental.data) #YES there is homogeneity of variances
oneway_test(Residuals_OB_Body~Group,data=placental.data)
placental.data$Group = factor(placental.data$Group, levels = c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups")) 
PT_Residuals_OB_Body<-pairwisePermutationTest(Residuals_OB_Body~Group,data=placental.data,method="fdr")
PT_Residuals_OB_Body
sink()

########

####Open data with PL brain residuals
data<-read.csv("Placental_data_final_PL_Br_res.csv", header=T, row.names = 1)
placental.data<-data[data$Abbreviations != "Pu_an" & data$Abbreviations != "Pt_mo", ]

##ggplot - boxplot - Residuals PL - Brain
ggboxplot(placental.data,x="Group", y="Residuals_PL_Brain", fill="Group", 
          palette=c("orangered1","turquoise","goldenrod1","blueviolet"),
          order=c("Mesozoic","Paleocene","Eocene stem groups","Eocene crown groups"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  font("ylab", size = 7.8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Age', y='Residuals log10(Petrosal lobule vs Endocranial volume)')

#Fisher-Pitman permutation test
sink("ANOVA_Welsh_log_PL_Brain_placental_only_topo1.txt")
shapiro.test(placental.data$Residuals_PL_Brain) # NO normal distribution
leveneTest(Residuals_PL_Brain ~ Group, data = placental.data) #NO there is homogeneity of variances
Welsh.test<-welch.test(Residuals_PL_Brain~Group,data=placental.data)
paircomp(Welsh.test)
sink()

########

####Open data with PL body residuals
data<-read.csv("Placental_data_final_PL_BM_res.csv", header=T, row.names = 1)
placental.data<-data[data$Abbreviations != "Pu_an" & data$Abbreviations != "Pt_mo", ]

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
sink("ANOVA_Welsh_log_PL_Body_placental_only_topo1.txt")
shapiro.test(placental.data$Residuals_PL_Body) # NO normal distribution
leveneTest(Residuals_PL_Body ~ Group, data = placental.data) #NO there is homogeneity of variances
Welsh.test<-welch.test(Residuals_PL_Body~Group,data=placental.data)
paircomp(Welsh.test)
sink()

### END
