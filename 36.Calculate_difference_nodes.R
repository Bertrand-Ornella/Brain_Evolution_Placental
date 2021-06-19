############################################################################################
############### Calculate difference between nodes to obtain negative values ###############
############################################################################################

setwd("~/Desktop/Placental_April_2021")

library(ape)
library(phytools)
library(strap)

####
### Open data - body mass
Results1<-read.csv("Results_Body_Br_placentals.csv", header=T, row.names = 1)
Results<-Results1[-c(1),]

#open calibrated tree
tree<-read.nexus("PEQ_tree.trees")
tree_scaled<-tree #need to be done because we will use tree to plot the ASR

##Rescale the branches of the timescaled phylogeny by mean scalar
scale<-Results$Median.Scalar/10000
tree_scaled$edge.length<-tree_scaled$edge.length*scale

#Estimate ancestral state reconstructions for this tree using ML 
Body.data1<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)
Body.data <-Body.data1[6]
Body.data <-log10(Body.data)

Body.data<-as.matrix(Body.data)
Body<-as.numeric(Body.data)
Body<-setNames(Body, row.names(Body.data))
anc<-fastAnc(tree_scaled,Body,vars=TRUE,CI=TRUE)
anc2<-anc$ace

#call edge matrix with brannch info
node_values <- tree_scaled$edge

#combine node and tip data
node_values[] <- c(Body[tree_scaled$tip.label], anc2)[tree_scaled$edge]

#calculate the difference
node_values = cbind(node_values, c(node_values[,2] - node_values[,1]))

#####

### Open data - brain size
Results1<-read.csv("Results_Brain_placentals.csv", header=T, row.names = 1)
Results<-Results1[-c(1),]

#open calibrated tree
tree<-read.nexus("PEQ_tree.trees")
tree_scaled<-tree #need to be done because we will use tree to plot the ASR

##Rescale the branches of the timescaled phylogeny by mean scalar
scale<-Results$Median.Scalar/10000
tree_scaled$edge.length<-tree_scaled$edge.length*scale

#Estimate ancestral state reconstructions for this tree using ML 
Brain.data1<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)
Brain.data <-Brain.data1[7]
Brain.data <-log10(Brain.data)

Brain.data<-as.matrix(Brain.data)
Brain<-as.numeric(Brain.data)
Brain<-setNames(Brain, row.names(Brain.data))
anc<-fastAnc(tree_scaled,Brain,vars=TRUE,CI=TRUE)
anc2<-anc$ace

#call edge matrix with brannch info
node_values <- tree_scaled$edge

#combine node and tip data
node_values[] <- c(Brain[tree_scaled$tip.label], anc2)[tree_scaled$edge]

#calculate the difference
node_values = cbind(node_values, c(node_values[,2] - node_values[,1]))

### END
