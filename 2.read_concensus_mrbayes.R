
###########################################################################################
################## Open MrBayes calibrated phylogeny to allow reading in R  ###############
###########################################################################################

############# This code collapses the nodes that MrBayes artificially created #############

#Before compiling, check tree and data to import, the nodes to collapse and trees, 
#location for removing extant species

library(ape) # read.nexus
library(geiger) # treedata
library(strap) # geoscalePhylo
library(phytools) # nodeHeights
library(paleotree)
library(castor)#merge_short_edges
library(strap) # geoscalePhylo

setwd("~/Desktop/Placental_April_2021")

#Import MrBayes tree
tree_after<-read.nexus("Placental1_mrbayes.nex.con.tre")

#Import tree
tree_original<-read.nexus("Placental_tree_1.nex")

#Adding a recent taxon to show same nodes numbers in the original tree
tree_before<-bind.tip(tree_original, "Recent_species", edge.length=NULL, where=124, position=0)

#check tree before
par(mfrow=c(1,1))
plotTree(tree_before,offset=1,fsize=0.3,ftype="i",lwd=1)
tiplabels(cex=0.3,frame="none",adj=c(-0.2,0.3))
nodelabels(cex=0.3,frame="none",adj=c(1.1,-0.4))

#check tree after
par(mfrow=c(1,1))
plotTree(tree_after,offset=1,fsize=0.3,ftype="i",lwd=1)
tiplabels(cex=0.3,frame="none",adj=c(-0.2,0.3))
nodelabels(cex=0.3,frame="none",adj=c(1.1,-0.4))

#Collapse nodes artificially indroduced by MrBayes, does not affect branch lengths 
tree_after_col<-merge_nodes_to_multifurcations(tree_after,nodes_to_merge=c(140,141,147,163,179,202,207,220)-length(tree_after$tip.label),merge_with_parents=TRUE,keep_ancestral_ages=TRUE)$tree

#check new tree after removing articial nodes
plotTree(tree_after_col,offset=1,fsize=0.3,ftype="i",lwd=1)
tiplabels(cex=0.3,frame="none",adj=c(-0.2,0.3))
nodelabels(cex=0.3,frame="none",adj=c(1.1,-0.4))

#compare trees to make sure nodes are the same
comparePhylo(tree_after_col, tree_before, plot=TRUE)

#Check that the branch length are the same
tree_after$edge.length
tree_after_col$edge.length

#Add root time argument with age of the root
tree_after_col$root.time<-max(nodeHeights(tree_after_col))
tree_after_col$root.time # works for PEQ, Brain and Body mass

#Remove the artificial recent tip
tree_drop<-drop.tip(tree_after_col, "Recent_species")

#Import placental ages to check if the tree and data match
placental.data<-read.csv("Placental_data_final.csv", header=T, row.names = 1)

#Order the data and tree, making sure the names are the same
data<-treedata(tree_drop,placental.data,sort=T,warnings=T)$data
tree_final<-treedata(tree_drop,placental.data,sort=T,warnings=T)$phy

name.check(tree_final, data, data.names=NULL)

#check tree
par(mfrow=c(1,1))
plot(tree_final, cex=0.3)

#Plot trees with geological scale
geoscalePhylo(tree=ladderize(tree_final,right=FALSE), units=c("Period", "Epoch", "Age"), boxes="Epoch",
              cex.tip=0.4, cex.age=0.3, cex.ts=0.5, label.offset=0, x.lim=c(-15,235), lwd=2, width=0.5)

#plot with ages
ages1<-as.data.frame(placental.data[,c(4,5)])
ages1$FAD_Ma<-as.numeric(ages1$FAD_Ma)
ages1$LAD_Ma<-as.numeric(ages1$LAD_Ma)
colnames(ages1)<-c("FAD", "LAD")

#Adding fossil ranges (ages=data) 
geoscalePhylo(tree=ladderize(tree_final,right=FALSE), ages=ages1, units=c("Period", "Epoch", "Age"), boxes="Epoch",
              cex.tip=0.35, cex.age=0.3, cex.ts=0.5, label.offset=0, x.lim=c(-15,235), lwd=2, width=0.5)

#save tree into nexus file to open the tree in Bayes Traits
write.nexus(tree_final, file = "Calibrated_placental1.trees")

#### END

