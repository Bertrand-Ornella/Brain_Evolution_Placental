
##########################################################################################
############################# PEQ per bins through time   ################################
##########################################################################################

#Before compiling, check data to import

library(ape)
library(geoscale)
library(strap)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/Placental_April_2021")

#Open data
placental.data<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)
#data<-read.csv("Placental_data_final_PEQ_Jackknife.csv", header=T, row.names = 1)

#Take out certain categories - jackknife
#placental.data<-data[data$Clade4 != "Perrissodacyla" & data$Clade4 != "Ferae", ]
#placental.data<-data[data$Clade4 != "Perrissodacyla", ]
#placental.data<-data[data$Clade4 != "Ferae", ]

#Obtain average for each FAD and LAD
FAD<-placental.data$FAD_Ma
LAD<-placental.data$LAD_Ma
Mean_age<-(FAD+LAD)/2
Age_df<-data.frame(Mean_age) 

max(Mean_age) #oldest bin average

#Create new columns for the mean age
placental.data[,"Age"] <- NA
colnames(Age_df)<-c("Age")
placental.data$Age=Age_df

#Age used to obtain the species that are in specific bins
myvars <- c("Age")
Ages <- placental.data[myvars]

#Number of bins
bin_num <- c(1:19)
desired_length <- 19
edge_num <- vector(mode = "list", length = desired_length)

#Create dataframe that will hold the results (species row number in each bins)
myDF <- data.frame(bin_num=bin_num, edge_num=cbind(edge_num))
AgeBin<- nrow(Ages)

#Create loop for every 5 millions of years, starting at the tree root (= 195.05 Mya)
for(loop_val in 1:AgeBin){
  
  ageval <- Ages[loop_val,1]
  
  for(innerLoopVal in 1:(length(bin_num)+1)){
    
    oldBinRange <- 195.05 - (innerLoopVal-1)*9.926923
    recentBinRange <- 195.05 - (innerLoopVal)*9.926923
  
    if(ageval >=recentBinRange & ageval <= oldBinRange)
    {
      
      getListOfBin<-myDF$edge_num[innerLoopVal]
      unlistedVector<-unlist(getListOfBin)
      
      unlistedVector<-append(unlistedVector,loop_val)
      myDF$edge_num[innerLoopVal] <- list(unlistedVector)
      
    }
  }
 }

### Reconstruct trait through time -- data should already be loaded
Trait<-placental.data[8]

#Define bins, get the length between root and 66Mya
195.05-66

#Choose number of bins and divide by this number to get the bin length for this section of time
129.05/13
#Make bins #19 total
UpperBin<-seq(195.05, 16.365386, by=-9.926923) 
LowerBin <-seq(185.123077, 6.438463, by=-9.926923)

#Make table to hold results
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

#Define variables for the loop (19 bins total)
Bins<-myDF$edge_num[1:19] #bin number
LengthBins<-length(Bins)

#Create loop
for(i in 1:LengthBins){
  
  Vals<-(Bins[[i]])
  
  Mean[i]<-mean(Trait[c(Vals),1])
  SD[i]<- sd(Trait[c(Vals),1])
  N[i]<-length(Vals)
  
  Left[i]<- Mean[i]-(qnorm(0.975)*SD[i]/sqrt(N[i]))
  Right[i]<- Mean[i]+(qnorm(0.975)*SD[i]/sqrt(N[i]))
 }

#save data in table to make graph and check numbers
Mean_df<-tibble(Mean)
colnames(Mean_df)<-c("Mean")
Table1$Mean=Mean

Left_df<-tibble(Left)
colnames(Left_df)<-c("Left")
Table1$Left=Left

Right_df<-tibble(Right)
colnames(Right_df)<-c("Right")
Table1$Right=Right

N_df<-tibble(N)
colnames(N_df)<-c("N")
Table1$N=N

#make subset of data to elimnate NaN
Table2<-Table1[c(1,4,7,12:17),] #19 bins ~10Mya

#Make a line for each bins on graph  
vline.data <- data.frame(z = seq(195.05,30, by=-9.926923)) # change to the lowest bin value

#Make plot
ggplot(Table2, aes(x=Average_Bin,y=Mean))+
  scale_x_reverse(breaks=seq(200,30, by=-20))+
  geom_vline(aes(xintercept = z), vline.data, linetype="dashed", color = "grey", size=0.1)+
  geom_vline(aes(xintercept = 66), linetype="dashed", color = "black", size=0.2)+
  geom_point()+
  geom_ribbon(aes(ymin=Left, ymax=Right), linetype=2, alpha=0.1)+
  geom_line(color="red", size=0.5)+
  labs(x='Millions of years ago', y='Phylogenetic Encephalization Quotient')+
  theme_classic()

#check taxa in bins with example for bin 8
#placental.data[c(8),]

################################## Jackknife dataset ####################################
#Extract data from Bins
myDF

#Bin 14
dat<-placental.data[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 16, 20, 21, 79, 80, 104, 111),]
nrow(dat)
#Bin 15
data<-placental.data[c(10, 11, 15, 17, 18, 19, 22, 23, 24, 26, 27, 28, 29, 30, 31, 32, 33, 49, 50, 51, 53, 54, 57, 58, 62, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 76, 77, 81, 82, 85, 86, 90, 91),]
nrow(data)
#Bin 16
data1<-placental.data[c(25, 34, 35, 36, 38, 40, 43, 44, 45, 47, 55, 56, 59, 63, 75, 78, 83, 87, 93, 94, 97, 98, 99, 100, 103),]
nrow(data1)
#Bin 17
data2<-placental.data[c(37, 39, 41, 42, 46, 48, 52, 60, 61, 64, 84, 88, 89, 92, 95, 96, 101, 102),]
nrow(data2)

###################################################
#Extract bin 14 data
Bin14<-dat[c(3,8,12:21)]
write.csv(Bin14,'Bin14.csv')

#Extract bin 15 data (groups and PEQ)
Bin15<-data[c(3,8,11:21)]
write.csv(Bin15,'Bin15.csv')

#Extract bin 16 data (groups and PEQ)
Bin16<-data1[c(3,8,11:21)]
write.csv(Bin16,'Bin16.csv')

#Extract bin 17 data (groups and PEQ)
Bin17<-data2[c(3,8,11:21)]
write.csv(Bin17,'Bin17.csv')

#### END
