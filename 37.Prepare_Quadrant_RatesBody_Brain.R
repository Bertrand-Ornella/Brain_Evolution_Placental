
############################################################################################
####################### Prepare data for Rate Body vs. Brain graph #########################
############################################################################################

setwd("~/Desktop/Placental_April_2021")

library(phytools)

#Open data
tree<-read.nexus("PEQ_tree.trees") #same for Brain size

par(mfrow=c(1,1))
plotTree(tree,offset=1,fsize=0.3,ftype="i",lwd=1)
tiplabels(cex=0.3,frame="none",adj=c(-0.2,0.3))
nodelabels(cex=0.3,frame="none",adj=c(1.1,-0.4))

#plot tree
plot(tree, cex=0.3)
edgelabels(cex=0.3,frame="none",adj=c(1,-0.5))

### Open data - body mass
Results<-read.csv("Results_Body_Br_placentals.csv", header=T, row.names = 1)

#Select data and log10 Median Scalar
Results$Median.Scalar<-log10(Results$Median.Scalar)
names(Results)[names(Results) == "Median.Scalar"] <- "Log.Median.Scalar.Body"
ResultsBody<-Results[-c(1),]
write.csv(ResultsBody,'Log_Body_Mass_rate_negative.csv')

### go add negative sign using data produced from 'Calculate_difference_nodes'

### Open data - brain size
Results1<-read.csv("Results_Brain_placentals.csv", header=T, row.names = 1)

Results1$Median.Scalar<-log10(Results1$Median.Scalar)
names(Results1)[names(Results1) == "Median.Scalar"] <- "Log.Median.Scalar.Brain"
ResultsBrain<-Results1[-c(1),]
write.csv(ResultsBrain,'Log_Brain_size_rate_negative.csv')

### go add negative sign using data produced from 'Calculate_difference_nodes'

##### Make dataset for graph
#Open data with negative values
ResultsBody<-read.csv("Log_Body_Mass_rate_negative.csv", header=T, row.names = 1)
Rate_Body<-ResultsBody[7]

#Open data with negative values
ResultsBrain<-read.csv("Log_Brain_size_rate_negative.csv", header=T, row.names = 1)

Rate_Brain<-ResultsBrain[7] 
List_taxa<-ResultsBrain[21]

#Make table including both median scalar
Data<- data.frame(Rate_Body=Rate_Body,Rate_Brain=Rate_Brain,List_taxa=List_taxa)
write.csv(Data,'Rates_Body_Brain_negative.csv')

### END


