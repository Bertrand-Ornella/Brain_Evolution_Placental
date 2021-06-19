
##########################################################################################
###############  Prepare tree and age for calibrating phylogeny in MrBayes  ##############
##########################################################################################

#Before compiling, check directory, data and tree to import, position where extant 
#species will be and its tip nb position

library(paleotree) # createMrBayesTipDatingNexus
library(phytools) # bind.tip
library(ape) # read.nexus
library(geiger) # treedata

setwd("~/Desktop/Placental_April_2021")

#Import placental ages
ages_tab<-read.csv("Placental_data_final.csv", header=T, row.names = 1)
n_ages_tab<-ages_tab[,4:5]

#Import tree
tree1<-read.nexus("Placental_tree_1.nex")

#Order the data and make sure the names are the same
data1<-treedata(tree1,n_ages_tab,sort=T,warnings=T)$data
tree1<-treedata(tree1,n_ages_tab,sort=T,warnings=T)$phy

#Adding a recent taxon to help clockless tip-dating analysis to find the root
treef1<-bind.tip(tree1, "Recent_species", edge.length=NULL, where=124, position=0)

#add recent species to dataset
data<-rbind(data1, c(0,0))
row.names(data)[125]<-"Recent_species"

#check tree
plot(treef1, cex=0.3)

#Make MrBayes file
createMrBayesTipDatingNexus(
  tipTimes = data,
  treeConstraints = treef1,
  anchorTaxon="Tetraclaenodon_pliciferus", 
  ageCalibrationType = "uniformRange",
  whichAppearance = "first",
  treeAgeOffset = 10,
  newFile = "Placental1_mrbayes.nex",
  createEmptyMorphMat = TRUE,
  doNotRun = FALSE
)

#Now manually add to nexus file: "sumt;" in the next line of the "sump;" command
#Then you can run this nexus file in MrBayes and get a consensus tree

#If you need Placentalia node at 66 Ma add manually this line to the MrBayes block: 
# "calibrate node13 = fixed (66);" [calibrate the Placentalia node at 66 Ma]
#ADD before [These taxa had fixed tip ages: Tetraclaenodon_pliciferus Recent_species ]
