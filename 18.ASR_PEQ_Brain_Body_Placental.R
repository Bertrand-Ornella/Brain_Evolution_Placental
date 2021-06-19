
##########################################################################################
################### ASR using Bayestraits for PEQ, Brain abd Body mass ###################
##########################################################################################

library(phytools)
library(geoscale)
library(strap)

setwd("~/Desktop/Placental_April_2021")

#tree root - see "2.read_consensus_mrbayes"
Tree.root.PEQ <-201.7734

######## Node numbering ############

tree<-read.nexus("PEQ_tree.trees")
par(mfrow=c(1,1))

plot(tree,no.margin=TRUE,edge.width=1,cex=0.3)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree),
           frame="none",adj=c(1.1,-0.4),cex=0.3)

tree$root.time <- Tree.root.PEQ
time.tree.PEQ<- geoscalePhylo(tree=ladderize(tree,right=FALSE), 
                              units=c("Period", "Epoch", "Age"), boxes="Epoch",
                              cex.tip=0.3, cex.age=0.5, cex.ts=0.5, label.offset=0, 
                              x.lim=c(-15,202), lwd=3, width=1)

tree<-read.nexus("Body_tree.trees")
par(mfrow=c(1,1))

plot(tree,no.margin=TRUE,edge.width=1,cex=0.3)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree),
           frame="none",adj=c(1.1,-0.4),cex=0.3)

tree$root.time <- Tree.root.PEQ
time.tree.PEQ<- geoscalePhylo(tree=ladderize(tree,right=FALSE), 
                              units=c("Period", "Epoch", "Age"), boxes="Epoch",
                              cex.tip=0.3, cex.age=0.5, cex.ts=0.5, label.offset=0, 
                              x.lim=c(-15,202), lwd=3, width=1)

########################  ASR under the variable rates model - PEQ  ########################
#open rates
Results1<-read.csv("Results_PEQ_placentals.csv", header=T, row.names = 1)
Results<-Results1[-c(1),]

#open calibrated tree
tree<-read.nexus("PEQ_tree.trees")
tree_scaled<-tree #need to be done because we will use tree to plot the ASR

##Rescale the branches of the timescaled phylogeny by mean scalar
scale<-Results$Median.Scalar
tree_scaled$edge.length<-tree_scaled$edge.length*scale

#Estimate ancestral state reconstructions for this tree using ML 
PEQ.data1<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)
PEQ.data <-PEQ.data1[8]

PEQ.data<-as.matrix(PEQ.data)
PEQ<-as.numeric(PEQ.data)
PEQ<-setNames(PEQ, row.names(PEQ.data))
anc<-fastAnc(tree_scaled,PEQ,vars=TRUE,CI=TRUE)
anc2<-anc$ace

#save ancestral state reconstruction in excel file
num.anc2<-round(anc2,digits =2)
write.table(num.anc2, file="PEQ_ancestral.csv")

#Plot these ancestral reconstructions in the regular time-scaled tree
plt<-fancyTree(tree, type="scattergram",X=as.matrix(PEQ),
               A=as.matrix(anc2),control=list(spin=FALSE),
               label="horizontal", plot=FALSE)

plt<-setMap(plt, colors=c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1"))
plt<-plt$contMaps[[1]]

#plot tree with node labels and time scale - !!! SAVE as A5 !!!
#Change "0.02" in ylim to change how tall the tree is
par(mfrow=c(1,1))
plot(plt, fsize=0.35,lwd=2.5,legend=FALSE,outline=FALSE,offset=1,
              ylim=c(-2-0.02*(Ntip(plt$tree)-1),Ntip(plt$tree)))

nodelabels(round(anc2,digits =2), cex=0.3, frame="none", adj = c(1.2, -0.75) )
tiplabels(round(PEQ[tree$tip.label],digits =2), cex=0.3, frame="none", adj= c(-0.4,0.3))

add.color.bar(leg=0.3*max(nodeHeights(tree)),plt$cols,title="PEQ",
              lims=plt$lims,digits=2,prompt=FALSE,lwd=5,fsize=0.8,subtitle="", x=5, y=25)

#Change "0.3" in mgp to change position of numbers on the scale bar (lower = higher)
#Change tck to change the length of the bars

axis(side=1, at=Tree.root.PEQ-15:Tree.root.PEQ, labels=15:Tree.root.PEQ,cex.axis=0.5, pos=-2,
     mgp=c(2, 0.2, 0), tck=-0.02)

#Now compare this ancestral reconstructions with reconstructions under a Brownian mode of evolution
anc_b<-fastAnc(tree,PEQ,vars=TRUE,CI=TRUE)
plb<-contMap(tree,PEQ, plot=FALSE)
plb<-setMap(plb, colors=c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1"))

plot(plb, fsize=0.4,lwd=2.5,outline=FALSE)
#nodelabels(round(anc_b$ace,digits =2), cex=0.5, frame="none", adj = c(1.2, -0.75) )

#To see differences in ancestral state reconstructions due the differences in the model assumptions 
anc2-anc_b$ace

########################  ASR under the variable rates model - Brain mass  ########################
#open rates
Results1<-read.csv("Results_Brain_placentals.csv", header=T, row.names = 1)
Results<-Results1[-c(1),]

#open calibrated tree
tree<-read.nexus("PEQ_tree.trees")
tree_scaled<-tree #need to be done because we will use tree to plot the ASR

##Rescale the branches of the timescaled phylogeny by mean scalar
scale<-Results$Median.Scalar
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

#save ancestral state reconstruction in excel file
num.anc2<-round(anc2,digits =2)
write.table(num.anc2, file="Brain_ancestral.csv")

#Plot these ancestral reconstructions in the regular time-scaled tree
plt<-fancyTree(tree,type="scattergram",X=as.matrix(Brain),
               A=as.matrix(anc2),control=list(spin=FALSE),
               label="horizontal", plot=FALSE)

plt<-setMap(plt, colors=c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1"))
plt<-plt$contMaps[[1]]

#plot tree with node labels and time scale - !!! SAVE as A5 !!!
#Change "0.02" in ylim to change how tall the tree is
par(mfrow=c(1,1))
plot(plt, fsize=0.35,lwd=2.5,legend=FALSE,outline=FALSE,offset=1,
     ylim=c(-2-0.02*(Ntip(plt$tree)-1),Ntip(plt$tree)))

nodelabels(round(anc2,digits =2), cex=0.3, frame="none", adj = c(1.2, -0.75) )
tiplabels(round(Brain[tree$tip.label],digits =2), cex=0.3, frame="none", adj= c(-0.4,0.3))

add.color.bar(leg=0.3*max(nodeHeights(tree)),plt$cols,title="Log10 (Endocranial volume)",
              lims=plt$lims,digits=2,prompt=FALSE,lwd=5,fsize=0.4,subtitle="", x=5, y=25)

#Change "0.3" in mgp to change position of numbers on the scale bar (lower = higher)
#Change tck to change the length of the bars

axis(side=1, at=Tree.root.PEQ-15:Tree.root.PEQ, labels=15:Tree.root.PEQ,cex.axis=0.5, pos=-2,
     mgp=c(2, 0.2, 0), tck=-0.02)

#Now compare this ancestral reconstructions with reconstructions under a Brownian mode of evolution
anc_b<-fastAnc(tree,Brain,vars=TRUE,CI=TRUE)
plb<-contMap(tree,Brain, plot=FALSE)
plb<-setMap(plb, colors=c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1"))

plot(plb, fsize=0.4,lwd=2.5,outline=FALSE)
#nodelabels(round(anc_b$ace,digits =2), cex=0.5, frame="none", adj = c(1.2, -0.75) )

#To see differences in ancestral state reconstructions due the differences in the model assumptions 
anc2-anc_b$ace

########################  ASR under the variable rates model - Body mass  ########################
#open rates
Results1<-read.csv("Results_Body_placentals.csv", header=T, row.names = 1)
Results<-Results1[-c(1),]

#open calibrated tree
tree<-read.nexus("Body_tree.trees")
tree_scaled<-tree #need to be done because we will use tree to plot the ASR

##Rescale the branches of the timescaled phylogeny by mean scalar
scale<-Results$Median.Scalar/10000
tree_scaled$edge.length<-tree_scaled$edge.length*scale

#Estimate ancestral state reconstructions for this tree using ML 
Body.data1<-read.csv("Placental_data_final.csv", header=T, row.names = 1)
Body.data <-Body.data1[7]
Body.data <-log10(Body.data)

Body.data<-as.matrix(Body.data)
Body<-as.numeric(Body.data)
Body<-setNames(Body, row.names(Body.data))
anc<-fastAnc(tree_scaled,Body,vars=TRUE,CI=TRUE)
anc2<-anc$ace

#save ancestral state reconstruction in excel file
num.anc2<-round(anc2,digits =2)
write.table(num.anc2, file="Body_ancestral.csv")

#Plot these ancestral reconstructions in the regular time-scaled tree
plt<-fancyTree(tree,type="scattergram",X=as.matrix(Body),
               A=as.matrix(anc2),control=list(spin=FALSE),
               label="horizontal", plot=FALSE)

plt<-setMap(plt, colors=c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1"))
plt<-plt$contMaps[[1]]

#plot tree with node labels and time scale - !!! SAVE as A5 !!!
#Change "0.02" in ylim to change how tall the tree is
par(mfrow=c(1,1))
plot(plt, fsize=0.35,lwd=2.5,legend=FALSE,outline=FALSE,offset=1,
     ylim=c(-2-0.02*(Ntip(plt$tree)-1),Ntip(plt$tree)))

nodelabels(round(anc2,digits =2), cex=0.3, frame="none", adj = c(1.2, -0.75) )
tiplabels(round(Body[tree$tip.label],digits =2), cex=0.3, frame="none", adj= c(-0.4,0.3))

add.color.bar(leg=0.3*max(nodeHeights(tree)),plt$cols,title="Log10 (Body mass)",
              lims=plt$lims,digits=2,prompt=FALSE,lwd=5,fsize=0.6,subtitle="", x=5, y=25)

#Change "0.3" in mgp to change position of numbers on the scale bar (lower = higher)
#Change tck to change the length of the bars

axis(side=1, at=Tree.root.PEQ-15:Tree.root.PEQ, labels=15:Tree.root.PEQ,cex.axis=0.5, pos=-2,
     mgp=c(2, 0.2, 0), tck=-0.02)

#Now compare this ancestral reconstructions with reconstructions under a Brownian mode of evolution
anc_b<-fastAnc(tree,Body,vars=TRUE,CI=TRUE)
plb<-contMap(tree,Body, plot=FALSE)
plb<-setMap(plb, colors=c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1"))

plot(plb, fsize=0.4,lwd=2.5,outline=FALSE)
#nodelabels(round(anc_b$ace,digits =2), cex=0.5, frame="none", adj = c(1.2, -0.75) )

#To see differences in ancestral state reconstructions due the differences in the model assumptions 
anc2-anc_b$ace

### END
