##########################################################################################
############### ASR using Bayestraits for the PL vs.BM and PL vs. Brain  #################
##########################################################################################

library(phytools)
library(geoscale)
library(strap)

setwd("~/Desktop/Placental_April_2021")

#Get root for tree
tree<-read.nexus("PEQ_tree.trees")
tree$root.time<-max(nodeHeights(tree)) # without extant taxon
PEQroot <-tree$root.time

tree<-read.nexus("PL_BM_tree.trees") #same for PL vs. Brain
tree$root.time<-max(nodeHeights(tree)) # without extant taxon
PLroot<-tree$root.time

#Tree root - PEQ, Brain and Body mass - see "2.read_consensus_mrbayes"
Tree.root.PEQ <-201.7734

#Tree root for PL
Tree.root.PL<-Tree.root.PEQ-(PEQroot-PLroot)
Tree.root.PL

######## Node numbering PL ############

tree<-read.nexus("PL_BM_tree.trees")

plot(tree,no.margin=TRUE,edge.width=1,cex=0.3)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree),
           frame="none",adj=c(1.1,-0.4),cex=0.3)

tree$root.time <- Tree.root.PL
time.tree.PL<- geoscalePhylo(tree=ladderize(tree,right=FALSE), 
                             units=c("Period", "Epoch", "Age"), boxes="Epoch",
                             cex.tip=0.3, cex.age=0.5, cex.ts=0.5, label.offset=0, 
                             x.lim=c(-15,202), lwd=3, width=1)

tree<-read.nexus("PL_Br_tree.trees")

plot(tree,no.margin=TRUE,edge.width=1,cex=0.3)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree),
           frame="none",adj=c(1.1,-0.4),cex=0.3)

tree$root.time <- Tree.root.PL
time.tree.PL<- geoscalePhylo(tree=ladderize(tree,right=FALSE), 
                             units=c("Period", "Epoch", "Age"), boxes="Epoch",
                             cex.tip=0.3, cex.age=0.5, cex.ts=0.5, label.offset=0, 
                             x.lim=c(-15,202), lwd=3, width=1)

####################  ASR under the variable rates model - PL-Body  ######################
#open rates
Results1<-read.csv("Results_PL_Body_placentals.csv", header=T, row.names = 1)
Results<-Results1[-c(1),]

#open calibrated tree
tree<-read.nexus("PL_BM_tree.trees")
tree_scaled<-tree #need to be done because we will use tree to plot the ASR

##Rescale the branches of the timescaled phylogeny by mean scalar
scale<-Results$Median.Scalar
tree_scaled$edge.length<-tree_scaled$edge.length*scale

#Estimate ancestral state reconstructions for this tree using ML 
PL.data1<-read.csv("Placental_data_final_PL_BM_res.csv", header=T, row.names = 1)
PL.data <-PL.data1[11] # body

PL.data<-as.matrix(PL.data)
PL<-as.numeric(PL.data)
PL<-setNames(PL, row.names(PL.data))
anc<-fastAnc(tree_scaled,PL,vars=TRUE,CI=TRUE)
anc2<-anc$ace

#save ancestral state reconstruction in excel file
num.anc2<-round(anc2,digits =2)
write.table(num.anc2, file="PL_Body_ancestral.csv")

#Plot these ancestral reconstructions in the regular time-scaled tree
plt<-fancyTree(tree,type="scattergram",X=as.matrix(PL),
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
tiplabels(round(PL[tree$tip.label],digits =2), cex=0.3, frame="none", adj= c(-0.4,0.3))

add.color.bar(leg=0.3*max(nodeHeights(tree)),plt$cols,title="Residuals[Log10 (Petrosal lobules vs. Body mass)]",
              lims=plt$lims,digits=2,prompt=FALSE,lwd=5,fsize=0.4,subtitle="", x=5, y=10)

#Change "0.3" in mgp to change position of numbers on the scale bar (lower = higher)
#Change tck to change the length of the bars

axis(side=1, at=Tree.root.PL-15:Tree.root.PL, labels=15:Tree.root.PL,cex.axis=0.5, pos=-2,
     mgp=c(2, 0.2, 0), tck=-0.02)

#Now compare this ancestral reconstructions with reconstructions under a Brownian mode of evolution
anc_b<-fastAnc(tree,PL,vars=TRUE,CI=TRUE)
plb<-contMap(tree,PL, plot=FALSE)
plb<-setMap(plb, colors=c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1"))

plot(plb, fsize=0.4,lwd=2.5,outline=FALSE)
#nodelabels(round(anc_b$ace,digits =2), cex=0.5, frame="none", adj = c(1.2, -0.75) )

#To see differences in ancestral state reconstructions due the differences in the model assumptions 
anc2-anc_b$ace

####################  ASR under the variable rates model - PL-Brain  ######################
#open rates
Results1<-read.csv("Results_PL_Brain_placentals.csv", header=T, row.names = 1)
Results<-Results1[-c(1),]

#open calibrated tree
tree<-read.nexus("PL_Br_tree.trees")
tree_scaled<-tree #need to be done because we will use tree to plot the ASR

##Rescale the branches of the timescaled phylogeny by mean scalar
scale<-Results$Median.Scalar
tree_scaled$edge.length<-tree_scaled$edge.length*scale

#Estimate ancestral state reconstructions for this tree using ML 
PL.data1<-read.csv("Placental_data_final_PL_Br_res.csv", header=T, row.names = 1)
PL.data <-PL.data1[12]

PL.data<-as.matrix(PL.data)
PL<-as.numeric(PL.data)
PL<-setNames(PL, row.names(PL.data))
anc<-fastAnc(tree_scaled,PL,vars=TRUE,CI=TRUE)
anc2<-anc$ace

#save ancestral state reconstruction in excel file
num.anc2<-round(anc2,digits =2)
write.table(num.anc2, file="PL_Brain_ancestral.csv")

#Plot these ancestral reconstructions in the regular time-scaled tree
plt<-fancyTree(tree,type="scattergram",X=as.matrix(PL),
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
tiplabels(round(PL[tree$tip.label],digits =2), cex=0.3, frame="none", adj= c(-0.4,0.3))

add.color.bar(leg=0.3*max(nodeHeights(tree)),plt$cols,title="Residuals[Log10 (Petrosal lobules vs. Endocranial volume)]",
              lims=plt$lims,digits=2,prompt=FALSE,lwd=5,fsize=0.35,subtitle="", x=5, y=10)

#Change "0.3" in mgp to change position of numbers on the scale bar (lower = higher)
#Change tck to change the length of the bars

axis(side=1, at=Tree.root.PL-15:Tree.root.PL, labels=15:Tree.root.PL,cex.axis=0.5, pos=-2,
     mgp=c(2, 0.2, 0), tck=-0.02)

#Now compare this ancestral reconstructions with reconstructions under a Brownian mode of evolution
anc_b<-fastAnc(tree,PL,vars=TRUE,CI=TRUE)
plb<-contMap(tree,PL, plot=FALSE)
plb<-setMap(plb, colors=c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1"))

plot(plb, fsize=0.4,lwd=2.5,outline=FALSE)
#nodelabels(round(anc_b$ace,digits =2), cex=0.5, frame="none", adj = c(1.2, -0.75) )

#To see differences in ancestral state reconstructions due the differences in the model assumptions 
anc2-anc_b$ace

### END
