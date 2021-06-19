
##########################################################################################
####################   Evolutionary rate through time - Body mass ########################
##########################################################################################

library(ape)
library(geoscale)
library(strap)
library(tibble)
library(ggplot2)

setwd("~/Desktop/Placental_April_2021")

#Open data
tree<-read.nexus("Body_tree.trees") # body mass

#plot tree
plot(tree, cex=0.3)

#Add root
tree$root.time <- 201.7734
time.tree<- geoscalePhylo(tree=ladderize(tree,right=FALSE), units=c("Period", "Epoch", "Age"), boxes="Epoch",
                          cex.tip=0.3, cex.age=0.5, cex.ts=0.5, label.offset=0, x.lim=c(-15,202), lwd=3, width=1)

#edgelabels(cex=0.3, frame="none", adj = c(1.2, -0.7) )
#nodelabels(cex=0.3, frame="none", adj = c(1.2, -0.7) )

### Open data - body mass
Results<-read.csv("Results_Body_placentals.csv", header=T, row.names = 1)

#Select data and log10 Median Scalar
Results$Median.Scalar<-log10(Results$Median.Scalar)
names(Results)[names(Results) == "Median.Scalar"] <- "Log.Median.Scalar"
Results1<-Results[-c(1),]
Trait_scalar<-Results1[7]

#Define bins, get the length between root and 66Mya
bindiv<-tree$root.time-66
#Choose number of bins and divide by this number to get the bin length for this section of time
bin_size<-bindiv/13
#sequence lower bin first number
seqlower<-tree$root.time-bin_size
#Make bins #18 total
UpperBin<-seq(tree$root.time, 24, by=-bin_size) 
LowerBin <-seq(seqlower, 13, by=-bin_size)

#######
treeedges<- tree$edge
lengthEdges<-nrow(treeedges)
treelengths<-tree$edge.length

#Create dataframe that will hold the results 
#branch length -- 1) branch number, 2) start position, 3) end position
edgelengtharray<-array(1:lengthEdges,dim=c(lengthEdges,3))

#Create loop
for (val in 1:dim(edgelengtharray)[1]){
  #get the start position of the edge which is the sum of its previous parent nodes
  #on the 1st column we always have the parent node
  #on the 2nd column we have the branch extending from that node
  #it is a new node branch if its value is >102 or the num of tips.
  #some node branches can create multiple nodes
  #how can we tell going backwards if something was its parent node?
  #we look at its parent node, then fill an array with the parent node of that node
  #and so on, until we reach the origin node value (103)
  
  #print(val)
  branchToInspect<-val
  
  #length of this particular branch
  branchlength<-treelengths[val]
  
  parentNodeOfBranch<-treeedges[branchToInspect,1]
  parentNodeOfBranch
  #totalLengthsOfPriorBranches
  TLOPB<-0
  #find row where parentNodeOfBranch is in column 2.
  #keep doing this until there are none found
  #add those lengths to TLOPB, that is start position
  lengthOfNode<-0
  indexOfParentNode<-0
  wExiter<-0
  while(wExiter==0){
    indexOfParentNode<-which(treeedges[,2]==parentNodeOfBranch)
    if(length(indexOfParentNode) !=0)
    {
      lengthOfNode<-treelengths[indexOfParentNode]
      TLOPB<-TLOPB+lengthOfNode
      parentNodeOfBranch<- treeedges[indexOfParentNode]
    }
    if(length(indexOfParentNode) ==0)
    {
      #print(TLOPB)
      edgelengtharray[val,2]<-TLOPB
      edgelengtharray[val,3]<-TLOPB+branchlength
      wExiter<-1
    }
  }
  
}

# Result showing 1) branch number, 2) start position, 3) end position
#edgelengtharray

#Defining the number of bins
bin_num <- c(1:18) #number of bins
desired_length <- 18 #number of bins
edge_num <- vector(mode = "list", length = desired_length)

#Create dataframe that will hold the results (branch row number in each bins)
myDF <- data.frame(bin_num=bin_num, edge_num=cbind(edge_num))
Binage<- nrow(edgelengtharray)

#Create the loop
for(loop_val in 1:Binage){
  #check the left
  #check the right
  #add to all the bins in between
  leftval <- edgelengtharray[loop_val,2]
  rightval <- edgelengtharray[loop_val,3]
  
  
  for(innerLoopVal in 1:(length(bin_num)+1))
  {
    bottomBinRange<- (innerLoopVal-1)*bin_size
    topBinRange <- (innerLoopVal *bin_size)
    
    #three conditions
    #left is greater than min, but right is within max
    #left is lower than min, but right is bigger than bottom
    
    if(rightval <=topBinRange & leftval >= bottomBinRange)
    {
      #add loop_val (the branch number) to the entry for innerLoopVal (bin 1, 2, 3, etc,)
      
      getListOfBin<-myDF$edge_num[innerLoopVal]
      unlistedVector<-unlist(getListOfBin)
      #adding branch num from outer loop
      unlistedVector<-append(unlistedVector,loop_val)
      #replacing list in myDF with the new list with this added branch value
      myDF$edge_num[innerLoopVal] <- list(unlistedVector)
      
    }
    else {
      if(leftval<bottomBinRange & rightval>=bottomBinRange)
      {
        #add loop_val (the branch number) to the entry for innerLoopVal (bin 1, 2, 3, etc,)
        
        getListOfBin<-myDF$edge_num[innerLoopVal]
        unlistedVector<-unlist(getListOfBin)
        #adding branch num from outer loop
        unlistedVector<-append(unlistedVector,loop_val)
        #replacing list in myDF with the new list with this added branch value
        myDF$edge_num[innerLoopVal] <- list(unlistedVector) 
      }
    }
  }
  
}

#Result contains the branch numbers inside each bin of 5 Millions of years
#myDF

###### make the graph
desired_length <- 18 #number of bins 
Table1 <- data.frame(UpperBin=UpperBin, LowerBin=LowerBin)
Upper<-Table1$UpperBin
Lower<-Table1$LowerBin
Average_Bin<-(Upper+Lower)/2
Average_Bin_df<-data.frame(Average_Bin)

Table1[,"Average_Bin"] <- NA
colnames(Average_Bin_df)<-c("Average_Bin_df")
Table1$Average_Bin=Average_Bin
Table1 <- data.frame(UpperBin=UpperBin, LowerBin=LowerBin, Average_Bin=Average_Bin)

Mean2 <- vector(mode = "list", length = desired_length)
Left2 <- vector(mode = "list", length = desired_length)
Right2 <- vector(mode = "list", length = desired_length)
SD2 <-vector(mode = "list", length = desired_length)
N2 <-vector(mode = "list", length = desired_length)

Mean<-as.numeric(unlist(Mean2))
Left<-as.numeric(unlist(Left2))
Right<-as.numeric(unlist(Right2))
SD<-as.numeric(unlist(SD2))
N<-as.numeric(unlist(N2))

Bins<-myDF$edge_num[1:18] #bin number
LengthBins<-length(Bins)

#Create loop
for(i in 1:LengthBins){
  
  Vals<-(Bins[[i]])
  
  Mean[i]<-mean(Trait_scalar[c(Vals),1])
  SD[i]<- sd(Trait_scalar[c(Vals),1])
  N[i]<-length(Vals)
  
  Left[i]<- Mean[i]-(qnorm(0.975)*SD[i]/sqrt(N[i]))
  Right[i]<- Mean[i]+(qnorm(0.975)*SD[i]/sqrt(N[i]))
}

#save data in table
Mean_df<-tibble(Mean)
colnames(Mean_df)<-c("Mean")
Table1$Mean=Mean

Left_df<-tibble(Left)
colnames(Left_df)<-c("Left")
Table1$Left=Left

Right_df<-tibble(Right)
colnames(Right_df)<-c("Right")
Table1$Right=Right

Table2<-Table1[c(1:17),]
write.csv(Table2,'Rate_Body_timebins.csv')

#Add dashed line for the different bins PEQ, Body mass, Brain size 
vline.data <- data.frame(z = seq(tree$root.time,20, by=-bin_size))

#Make plot 
ggplot(Table2, aes(x=Average_Bin,y=Mean))+
  scale_x_reverse(breaks=seq(200,20, by=-20))+
  geom_vline(aes(xintercept = z), vline.data, linetype="dashed", color = "grey", size=0.1)+
  geom_vline(aes(xintercept = 66), linetype="dashed", color = "black", size=0.2)+
  geom_point()+
  geom_ribbon(aes(ymin=Left, ymax=Right), linetype=2, alpha=0.1)+
  geom_line(color="red", size=0.5)+
  theme_classic(base_size = 10)+
  labs(x='Millions of years ago', y='log10 (Body mass rate of evolution)')

#### END
