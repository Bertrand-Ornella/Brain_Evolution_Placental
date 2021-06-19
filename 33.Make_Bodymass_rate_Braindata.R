
##########################################################################################
##########################   Make dataset for Bayestraits   ##############################
##########################################################################################

library(phytools)
library(geiger) #name.check

setwd("~/Desktop/Placental_April_2021")

#########################  Bayesian dataset and tree for PEQ  #############################

#### Open data with Body mass
Body.data1<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)
Body.data <-Body.data1[6]

#Save Body data
write.table(x = Body.data, file = "Body_Br_ratedata.txt", sep = "\t")

#Open this file
data_Body2<-read.table("Body_Br_ratedata.txt", sep="")

#Save data to be used in Bayestraits
write.table(data_Body2,"Body_Br_rate.txt",col.names = FALSE, row.names = TRUE, sep=" ",fileEncoding="UTF-8", quote = FALSE)

#### Make tree for Bayestraits for Body mass
#Open tree
PEQ_tree<-read.nexus("PEQ_tree.trees")

#Check if the data and tree have the same names
name.check(PEQ_tree, data_Body2, data.names=NULL)

