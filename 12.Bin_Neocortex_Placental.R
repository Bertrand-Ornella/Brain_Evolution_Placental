
##########################################################################################
##########################   Neocortex per bins through time   ############################
##########################################################################################

#Before compiling, check data to import

#Check out tree
library(ape)
library(geoscale)
library(strap)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/Placental_April_2021")

#Open data
placental.data<-read.csv("Placental_data_final_Neocortex_res.csv", header=T, row.names = 1)

#Obtain average for each FAD and LAD
FAD<-placental.data$FAD_Ma
LAD<-placental.data$LAD_Ma
Mean_age<-(FAD+LAD)/2
Age_df<-data.frame(Mean_age)

#Create new columns
placental.data[,"Age"] <- NA
colnames(Age_df)<-c("Age")
placental.data$Age=Age_df

#Age used to obtain the species that are in specific bins
myvars <- c("Age")
Ages <- placental.data[myvars]

#Number of bins
bin_num <- c(1:9)
desired_length <- 9
edge_num <- vector(mode = "list", length = desired_length)

#Create dataframe that will hold the results (species row number in each bins)
myDF <- data.frame(bin_num=bin_num, edge_num=cbind(edge_num))
AgeBin<- nrow(Ages)

for(loop_val in 1:AgeBin){
  
  ageval <- Ages[loop_val,1]
  
  for(innerLoopVal in 1:(length(bin_num)+1)){
    
    oldBinRange <- 66 - (innerLoopVal-1)*5
    recentBinRange <- 66 - (innerLoopVal)*5
    
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
Trait<-placental.data[12] # Neocortex residuals

#Make bins #9 total
UpperBin<-seq(66, 26, by=-5) 
LowerBin <-seq(61, 21, by=-5)

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

#Define variables for the loop
Bins<-myDF$edge_num[1:9] #bin number
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
Table2<-Table1[c(1:7),]

#Make a line for each bins on graph   
vline.data <- data.frame(z = seq(66,21, by=-5)) # change to the lowest bin value

#Make plot
ggplot(Table2, aes(x=Average_Bin,y=Mean))+
  scale_x_reverse(breaks=seq(70,30, by=-5))+
  geom_vline(aes(xintercept = z), vline.data, linetype="dashed", color = "grey", size=0.1)+
  geom_vline(aes(xintercept = 66), linetype="dashed", color = "black", size=0.2)+
  geom_point()+
  geom_ribbon(aes(ymin=Left, ymax=Right), linetype=2, alpha=0.1)+
  geom_line(color="red", size=0.5)+
  labs(x='Millions of years ago', y='Residuals log10(Neocortex vs. Endocranial surface area)')+
  theme_classic(base_size = 9.5)

#check taxa in bins - example
#myDF
#placental.data[c(3),]

#### END
