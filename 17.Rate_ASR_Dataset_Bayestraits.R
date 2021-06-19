
##########################################################################################
##########################   Make dataset for Bayestraits   ##############################
##########################################################################################

library(phytools)
library(geiger) #name.check

setwd("~/Desktop/Placental_April_2021")

########################  Bayesian dataset PEQ and Brain size ###########################

###!!"Placental_data_final_PEQ.csv" + "PEQ_tree.trees" created in "PGLS_PEQ_Placental"!!

#### Open data PEQ
PEQ.data1<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)
PEQ.data <-PEQ.data1[8]

#Open tree
PEQ_tree<-read.nexus("PEQ_tree.trees")

#Save PEQ data
write.table(x = PEQ.data, file = "PEQ_ratedata.txt", sep = "\t")

#Open this file
data_PEQ2<-read.table("PEQ_ratedata.txt", sep="")

#Save data to be used in Bayestraits
write.table(data_PEQ2,"PEQ_rate.txt",col.names = FALSE, row.names = TRUE, sep=" ",fileEncoding="UTF-8", quote = FALSE)

#Check if the data and tree have the same names
name.check(PEQ_tree, data_PEQ2, data.names=NULL)

#### Open data with Brain size
Brain.data1<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)
Brain.data <-Brain.data1[7]

#Save Brain data
write.table(x = Brain.data, file = "Brain_ratedata.txt", sep = "\t")

#Open this file
data_Brain2<-read.table("Brain_ratedata.txt", sep="")

#Save data to be used in Bayestraits
write.table(data_Brain2,"Brain_rate.txt",col.names = FALSE, row.names = TRUE, sep=" ",fileEncoding="UTF-8", quote = FALSE)

#Check if the data and tree have the same names
name.check(PEQ_tree, data_Brain2, data.names=NULL)

#######################  Bayesian dataset and tree for Body mass  #########################

#### Open data with Body mass
Body.data1<-read.csv("Placental_data_final.csv", header=T, row.names = 1)
Body.data <-Body.data1[7]

#Save Body data
write.table(x = Body.data, file = "Body_ratedata.txt", sep = "\t")

#Open this file
data_Body2<-read.table("Body_ratedata.txt", sep="")

#Save data to be used in Bayestraits
write.table(data_Body2,"Body_rate.txt",col.names = FALSE, row.names = TRUE, sep=" ",fileEncoding="UTF-8", quote = FALSE)

#### Make tree for Bayestraits for Body mass
Body_tree<-read.nexus("Calibrated_placental1.trees")

#save tree into nexus file to open the tree in Bayes Traits
write.nexus(Body_tree, file = "Body_tree.trees")

#Check if the data and tree have the same names
name.check(Body_tree, data_Body2, data.names=NULL)

##########################  Bayesian dataset for Olfactory bulbs ##########################

##!!"Placental_data_final_OB_res.csv" + "OB_tree.trees" created in "PGLS_OB_Placental"!!

#Open data with OB size
OB.data1<-read.csv("Placental_data_final_OB_res.csv", header=T, row.names = 1)
OB_Brain.data <-OB.data1[14] # Residuals - Olfactory bulb - Brain size
OB_Body.data <-OB.data1[16] # Residuals - Olfactory bulb - Body mass

#Save OB data
write.table(x = OB_Body.data, file = "OB_Body_ratedata.txt", sep = "\t")
write.table(x = OB_Brain.data, file = "OB_Brain_ratedata.txt", sep = "\t")

#Open this file
data_OB_Body2<-read.table("OB_Body_ratedata.txt", sep="")
data_OB_Brain2<-read.table("OB_Brain_ratedata.txt", sep="")

#Save data to be used in Bayestraits
write.table(data_OB_Body2,"OB_Body_rate.txt",col.names = FALSE, row.names = TRUE, sep=" ",fileEncoding="UTF-8", quote = FALSE)
write.table(data_OB_Brain2,"OB_Brain_rate.txt",col.names = FALSE, row.names = TRUE, sep=" ",fileEncoding="UTF-8", quote = FALSE)

#Open tree
OB_tree<-read.nexus("OB_tree.trees")

#Check if the data and tree have the same names
name.check(OB_tree, data_OB_Body2, data.names=NULL)

#Check if the data and tree have the same names
name.check(OB_tree, data_OB_Brain2, data.names=NULL)

##########################  Bayesian dataset for Petrosal lobules ##########################

##!!"Placental_data_final_PL_BM_res.csv" and "Placental_data_final_PL_Br_res.csv" 
## + "PL_BM_tree.trees" and "PL_Br_tree.trees" created in "PGLS_PL_Placental"!!

#Open data with PL size
PL.data_Body<-read.csv("Placental_data_final_PL_BM_res.csv", header=T, row.names = 1)
PL.data_Brain<-read.csv("Placental_data_final_PL_Br_res.csv", header=T, row.names = 1)

PL_Body.data <-PL.data_Body[11] # Residuals - Petrosal lobule - Body mass
PL_Brain.data <-PL.data_Brain[12] # Residuals - Petrosal lobule - Brain

#Save PL data
write.table(x = PL_Body.data, file = "PL_Body_ratedata.txt", sep = "\t")
write.table(x = PL_Brain.data, file = "PL_Brain_ratedata.txt", sep = "\t")

#Open this file
data_PL_Body2<-read.table("PL_Body_ratedata.txt", sep="")
data_PL_Brain2<-read.table("PL_Brain_ratedata.txt", sep="")

#Save data to be used in Bayestraits
write.table(data_PL_Body2,"PL_Body_rate.txt",col.names = FALSE, row.names = TRUE, sep=" ",fileEncoding="UTF-8", quote = FALSE)
write.table(data_PL_Brain2,"PL_Brain_rate.txt",col.names = FALSE, row.names = TRUE, sep=" ",fileEncoding="UTF-8", quote = FALSE)

#Open tree
PL_tree_Body<-read.nexus("PL_BM_tree.trees")
PL_tree_Brain<-read.nexus("PL_Br_tree.trees")

#Check if the data and tree have the same names
name.check(PL_tree_Body, data_PL_Body2, data.names=NULL)

#Check if the data and tree have the same names
name.check(PL_tree_Brain, data_PL_Brain2, data.names=NULL)

### END





