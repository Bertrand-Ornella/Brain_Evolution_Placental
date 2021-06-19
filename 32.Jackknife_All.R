##########################################################################################
####################################### PEQ Jackknife ####################################
##########################################################################################

library(bootstrap)

setwd("~/Desktop/Placental_April_2021")

################################# Crown vs. Stem taxa ##################################
#import bin data
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data into crown vs stem
stem.data<- Bin15.data[ which(Bin15.data$Group=='Eocene stem groups'), ]
crown.data<- Bin15.data[ which(Bin15.data$Group=='Eocene crown groups'), ]

stem.data1<- Bin16.data[ which(Bin16.data$Group=='Eocene stem groups'), ]
crown.data1<- Bin16.data[ which(Bin16.data$Group=='Eocene crown groups'), ]

stem.data2<- Bin17.data[ which(Bin17.data$Group=='Eocene stem groups'), ]
crown.data2<- Bin17.data[ which(Bin17.data$Group=='Eocene crown groups'), ]

#ANOVA of original data Bin 15
t.all <- t.test(crown.data$PEQ, stem.data$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all

#ANOVA of original data Bin 16
t.all1 <- t.test(crown.data1$PEQ, stem.data1$PEQ, alternative="greater")
t.all.p1 <- t.all1$p.value
t.all1

#ANOVA of original data Bin 17
t.all2 <- t.test(crown.data2$PEQ, stem.data2$PEQ, alternative="greater")
t.all.p2 <- t.all2$p.value
t.all2

## Results: Only in Bin15, crown and stem are significantly different

#check sample size in each category
nrow(stem.data)
nrow(crown.data)

nrow(stem.data1)
nrow(crown.data1)

nrow(stem.data2)
nrow(crown.data2)

#preset empty numeric vector to store output statistic in
#number of replicates is number of leave-one-out jackknife samples times 1000 bootstraps
#using the smallest subset
replicates_Bin15 <- numeric(length = 1000*nrow(crown.data))
replicates_Bin16 <- numeric(length = 1000*nrow(stem.data1))
replicates_Bin17 <- numeric(length = 1000*nrow(stem.data2))

#-take crown data (the smaller subset) and loop through leave-one-out combinations
#-for each jackknife set, bootstrap sample the larger dataset without replacement 
#with same sample size.
#-then calculate t test and record p-value

#Bin 15
for (i in 1:nrow(crown.data)) {
  jackknife_sample <- crown.data[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(stem.data$PEQ, size=nrow(crown.data)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_Bin15[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_Bin15, breaks=20)
abline(v=t.all.p, col="red")

#Bin 16
for (i in 1:nrow(stem.data1)) {
  jackknife_sample <- stem.data1[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(crown.data1$PEQ, size=nrow(stem.data1)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_Bin16[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_Bin16, breaks=20)
abline(v=t.all.p1, col="red")

#Bin 17
for (i in 1:nrow(stem.data2)) {
  jackknife_sample <- stem.data2[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(crown.data2$PEQ, size=nrow(stem.data2)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_Bin17[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_Bin17, breaks=20)
abline(v=t.all.p2, col="red")

########################################################################################
######### Crown 15 vs. 16 vs. 17

#import bin data
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data into crown vs stem
crown.data<- Bin15.data[ which(Bin15.data$Group=='Eocene crown groups'), ]
crown.data1<- Bin16.data[ which(Bin16.data$Group=='Eocene crown groups'), ]
crown.data2<- Bin17.data[ which(Bin17.data$Group=='Eocene crown groups'), ]

#ANOVA of original data 15-16
t.all <- t.test(crown.data$PEQ, crown.data1$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all

#ANOVA of original data 16-17
t.all1 <- t.test(crown.data1$PEQ, crown.data2$PEQ, alternative="greater")
t.all.p1 <- t.all1$p.value
t.all1

#check sample size in each category
nrow(crown.data)
nrow(crown.data1)
nrow(crown.data2)

#using the smallest subset
replicates_15_16_Crown <- numeric(length = 1000*nrow(crown.data))
replicates_16_17_Crown <- numeric(length = 1000*nrow(crown.data2))

## Replicates Bin 15 vs. Bin 16
for (i in 1:nrow(crown.data)) {
  jackknife_sample <- crown.data[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(crown.data1$PEQ, size=nrow(crown.data)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_16_Crown[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_16_Crown, breaks=20)
abline(v=t.all.p, col="red")

## Replicates Bin 16 vs. Bin 17
for (i in 1:nrow(crown.data2)) {
  jackknife_sample <- crown.data2[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(crown.data1$PEQ, size=nrow(crown.data2)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_16_17_Crown[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_16_17_Crown, breaks=20)
abline(v=t.all.p1, col="red")

#########################################################
######### Stem 15 vs. 16 vs. 17

#import bin data
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data into crown vs stem
stem.data<- Bin15.data[ which(Bin15.data$Group=='Eocene stem groups'), ]
stem.data1<- Bin16.data[ which(Bin16.data$Group=='Eocene stem groups'), ]
stem.data2<- Bin17.data[ which(Bin17.data$Group=='Eocene stem groups'), ]

#ANOVA of original data 15-16
t.all <- t.test(stem.data$PEQ, stem.data1$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all

#ANOVA of original data 16-17
t.all1 <- t.test(stem.data1$PEQ, stem.data2$PEQ, alternative="greater")
t.all.p1 <- t.all1$p.value
t.all1

#check sample size in each category
nrow(stem.data)
nrow(stem.data1)
nrow(stem.data2)

#using the smallest subset
replicates_15_16_Stem <- numeric(length = 1000*nrow(stem.data1))
replicates_16_17_Stem <- numeric(length = 1000*nrow(stem.data2))

## Replicates Bin 15 vs. Bin 16
for (i in 1:nrow(stem.data1)) {
  jackknife_sample <- stem.data1[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(stem.data$PEQ, size=nrow(stem.data1)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_16_Stem[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_16_Stem, breaks=20)
abline(v=t.all.p, col="red")

## Replicates Bin 16 vs. Bin 17
for (i in 1:nrow(stem.data2)) {
  jackknife_sample <- stem.data2[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(stem.data1$PEQ, size=nrow(stem.data2)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_16_17_Stem[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_16_17_Stem, breaks=20)
abline(v=t.all.p1, col="red")

############################### Clade1A Crown Cetartiodactyla ############################

#import bin data
Bin14.data<-read.csv("Bin14.csv", header=T, row.names = 1)
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data
Cet.15 <-Bin15.data[ which(Bin15.data$Clade1A=='Cetartiodactyla'), ]
Cet.16 <-Bin16.data[ which(Bin16.data$Clade1A=='Cetartiodactyla'), ]
Cet.17 <-Bin17.data[ which(Bin17.data$Clade1A=='Cetartiodactyla'), ]

#ANOVA of original data Bin 15 - 16
t.all <- t.test(Cet.15$PEQ, Cet.16$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all

#ANOVA of original data Bin 16 - 17
t.all1 <- t.test(Cet.16$PEQ, Cet.17$PEQ, alternative="greater")
t.all.p1 <- t.all1$p.value
t.all1

### No significant difference between consecutive bins for Crown Cetartiodactyla

#check sample size in each category
nrow(Cet.15) #2 Not enough taxa to it
nrow(Cet.16)
nrow(Cet.17)

#using the smallest subset
replicates_16_17_Cet <- numeric(length = 1000*nrow(Cet.17))

for (i in 1:nrow(Cet.17)) {
  jackknife_sample <- Cet.17[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Cet.16$PEQ, size=nrow(Cet.17)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_16_17_Cet[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_16_17_Cet, breaks=20)
abline(v=t.all.p, col="red")
### No significant difference from bin 16 to 17 for Crown Cetartiodactyla

########################## Clade1A Crown Perrissodactyla #####################

#import bin data
Bin14.data<-read.csv("Bin14.csv", header=T, row.names = 1)
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data
Per.15 <-Bin15.data[ which(Bin15.data$Clade1B=='Perrissodacyla'), ]
Per.16 <-Bin16.data[ which(Bin16.data$Clade1B=='Perrissodacyla'), ]
Per.17 <-Bin17.data[ which(Bin17.data$Clade1B=='Perrissodacyla'), ]

#check sample size in each category
nrow(Per.15) 
nrow(Per.16) #1 not enough 
nrow(Per.17) #0 
### Not enough observation to run test

################# Clade 2A - Stem ungulates Cet vs. Crown Cetartiodactyla ################

#import bin data
Bin14.data<-read.csv("Bin14.csv", header=T, row.names = 1)
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data
SUC.14 <-Bin14.data[ which(Bin14.data$Clade2A=='Stem ungulate'), ]
SUC.15 <-Bin15.data[ which(Bin15.data$Clade2A=='Stem ungulate'), ]

CC.15 <-Bin15.data[ which(Bin15.data$Clade2A=='Cetartiodactyla'), ]
CC.16 <-Bin16.data[ which(Bin16.data$Clade2A=='Cetartiodactyla'), ]
CC.17 <-Bin17.data[ which(Bin17.data$Clade2A=='Cetartiodactyla'), ]

#check sample size in each category
nrow(CC.15) #2
nrow(CC.16) # no SUC for this bin
nrow(CC.17) # no SUC for this bin
nrow(SUC.14) # no CC for this bin
nrow(SUC.15) #1
### Not enough observation to run test

################# Clade 2B - Stem ungulates Per vs. Crown Perrissodactyla ################

#import bin data
Bin14.data<-read.csv("Bin14.csv", header=T, row.names = 1)
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data
CP.15 <-Bin15.data[ which(Bin15.data$Clade2B=='Perrissodacyla'), ]
CP.16 <-Bin16.data[ which(Bin16.data$Clade2B=='Perrissodacyla'), ]

SUP.14 <-Bin14.data[ which(Bin14.data$Clade2B=='Stem ungulate'), ]
SUP.15 <-Bin15.data[ which(Bin15.data$Clade2B=='Stem ungulate'), ]
SUP.16 <-Bin16.data[ which(Bin16.data$Clade2B=='Stem ungulate'), ]
SUP.17 <-Bin17.data[ which(Bin17.data$Clade2B=='Stem ungulate'), ]

#check sample size in each category
nrow(CP.15) #6
nrow(CP.16) #1
nrow(SUP.14) #2
nrow(SUP.15) #6
nrow(SUP.16) #2
nrow(SUP.17) #1

#ANOVA of original data
t.all <- t.test(CP.15$PEQ, SUP.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all## YES - Stem and Crown perrissodactyles are significantly different in Bin 15

#using the smallest subset
replicates_15_CP <- numeric(length = 1000*nrow(CP.15))
replicates_15_SUP <- numeric(length = 1000*nrow(SUP.15))

#### CP.15 considered the smaller sample
for (i in 1:nrow(CP.15)) {
  jackknife_sample <- CP.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(SUP.15$PEQ, size=nrow(CP.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_CP[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_CP, breaks=20)
abline(v=t.all.p, col="red")

#### SUP.15 considered the smaller sample
for (i in 1:nrow(SUP.15)) {
  jackknife_sample <- SUP.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(CP.15$PEQ, size=nrow(SUP.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_SUP[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_SUP, breaks=20)
abline(v=t.all.p, col="red")

### Does not appears to hold with the Jackknife but maybe yes because all values are
### below 0.05
#When we choose Stem as the smaller values, it is different.

################# Clade 2C - All Cetartiodactyla vs. All Perrissodactyla ################

#import bin data
Bin14.data<-read.csv("Bin14.csv", header=T, row.names = 1)
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data
Per.14 <-Bin14.data[ which(Bin14.data$Clade2C=='Perrissodacyla'), ]
Per.15 <-Bin15.data[ which(Bin15.data$Clade2C=='Perrissodacyla'), ]
Per.16 <-Bin16.data[ which(Bin16.data$Clade2C=='Perrissodacyla'), ]
Per.17 <-Bin17.data[ which(Bin17.data$Clade2C=='Perrissodacyla'), ]

Cet.14 <-Bin14.data[ which(Bin14.data$Clade2C=='Cetartiodactyla'), ]
Cet.15 <-Bin15.data[ which(Bin15.data$Clade2C=='Cetartiodactyla'), ]
Cet.16 <-Bin16.data[ which(Bin16.data$Clade2C=='Cetartiodactyla'), ]
Cet.17 <-Bin17.data[ which(Bin17.data$Clade2C=='Cetartiodactyla'), ]

#check sample size in each category
nrow(Per.14) #2 - not enough for jackknife
nrow(Per.15) #12
nrow(Per.16) #3
nrow(Per.17) #1 - not enough

nrow(Cet.14) #7
nrow(Cet.15) #3
nrow(Cet.16) #13
nrow(Cet.17) #9

#ANOVA of original data Bin 15
t.all <- t.test(Per.15$PEQ, Cet.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_15_Per_Cet <- numeric(length = 1000*nrow(Cet.15))

for (i in 1:nrow(Cet.15)) {
  jackknife_sample <- Cet.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Per.15$PEQ, size=nrow(Cet.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_Per_Cet[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_Per_Cet, breaks=20)
abline(v=t.all.p, col="red")

#ANOVA of original data Bin 16
t.all <- t.test(Per.16$PEQ, Cet.16$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_16_Per_Cet <- numeric(length = 1000*nrow(Per.16))

for (i in 1:nrow(Per.16)) {
  jackknife_sample <- Per.16[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Cet.16$PEQ, size=nrow(Per.16)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_16_Per_Cet[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_16_Per_Cet, breaks=20)
abline(v=t.all.p, col="red")

################################## All ungulates vs. other mammals ###################### 

#import bin data
Bin14.data<-read.csv("Bin14.csv", header=T, row.names = 1)
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data
Ung.14 <-Bin14.data[ which(Bin14.data$Clade2D=='Ungulate'), ]
Ung.15 <-Bin15.data[ which(Bin15.data$Clade2D=='Ungulate'), ]
Ung.16 <-Bin16.data[ which(Bin16.data$Clade2D=='Ungulate'), ]
Ung.17 <-Bin17.data[ which(Bin17.data$Clade2D=='Ungulate'), ]

Mam.14 <-Bin14.data[ which(Bin14.data$Clade2D=='Mammalia'), ]
Mam.15 <-Bin15.data[ which(Bin15.data$Clade2D=='Mammalia'), ]
Mam.16 <-Bin16.data[ which(Bin16.data$Clade2D=='Mammalia'), ]
Mam.17 <-Bin17.data[ which(Bin17.data$Clade2D=='Mammalia'), ]

#check sample size in each category
nrow(Ung.14) #9
nrow(Ung.15) #15
nrow(Ung.16) #16
nrow(Ung.17) #10

nrow(Mam.14) #8
nrow(Mam.15) #28
nrow(Mam.16) #9
nrow(Mam.17) #8

#ANOVA of original data Bin 14
t.all <- t.test(Ung.14$PEQ, Mam.14$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_14_Ung_Mam <- numeric(length = 1000*nrow(Mam.14))

for (i in 1:nrow(Mam.14)) {
  jackknife_sample <- Mam.14[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Ung.14$PEQ, size=nrow(Mam.14)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_14_Ung_Mam[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_14_Ung_Mam, breaks=20)
abline(v=t.all.p, col="red")

#ANOVA of original data Bin 15
t.all <- t.test(Ung.15$PEQ, Mam.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_15_Ung_Mam <- numeric(length = 1000*nrow(Ung.15))

for (i in 1:nrow(Ung.15)) {
  jackknife_sample <- Ung.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Mam.15$PEQ, size=nrow(Ung.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_Ung_Mam[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_Ung_Mam, breaks=20)
abline(v=t.all.p, col="red")

#ANOVA of original data Bin 16
t.all <- t.test(Ung.16$PEQ, Mam.16$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_16_Ung_Mam <- numeric(length = 1000*nrow(Mam.16))

for (i in 1:nrow(Mam.16)) {
  jackknife_sample <- Mam.16[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Ung.16$PEQ, size=nrow(Mam.16)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_16_Ung_Mam[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_16_Ung_Mam, breaks=20)
abline(v=t.all.p, col="red")

#ANOVA of original data Bin 17
t.all <- t.test(Ung.17$PEQ, Mam.17$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_17_Ung_Mam <- numeric(length = 1000*nrow(Mam.17))

for (i in 1:nrow(Mam.17)) {
  jackknife_sample <- Mam.17[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Ung.17$PEQ, size=nrow(Mam.17)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_17_Ung_Mam[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_17_Ung_Mam, breaks=20)
abline(v=t.all.p, col="red")

################################ Stem ungulates vs. Crown ungulates ###################### 

#import bin data
Bin14.data<-read.csv("Bin14.csv", header=T, row.names = 1)
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data
SUng.14 <-Bin14.data[ which(Bin14.data$Clade3=='Stem ungulate'), ]
SUng.15 <-Bin15.data[ which(Bin15.data$Clade3=='Stem ungulate'), ]
SUng.16 <-Bin16.data[ which(Bin16.data$Clade3=='Stem ungulate'), ]
SUng.17 <-Bin17.data[ which(Bin17.data$Clade3=='Stem ungulate'), ]

CUng.14 <-Bin14.data[ which(Bin14.data$Clade3=='Crown ungulate'), ]
CUng.15 <-Bin15.data[ which(Bin15.data$Clade3=='Crown ungulate'), ]
CUng.16 <-Bin16.data[ which(Bin16.data$Clade3=='Crown ungulate'), ]
CUng.17 <-Bin17.data[ which(Bin17.data$Clade3=='Crown ungulate'), ]

#check sample size in each category
nrow(SUng.14) #9
nrow(SUng.15) #7
nrow(SUng.16) #2 - not enough
nrow(SUng.17) #1 - not enough

nrow(CUng.14) #0 - not enough
nrow(CUng.15) #8
nrow(CUng.16) #14
nrow(CUng.17) #9

#ANOVA of original data Bin 15
t.all <- t.test(SUng.15$PEQ, CUng.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_15_SUng_CUng <- numeric(length = 1000*nrow(SUng.15))

for (i in 1:nrow(SUng.15)) {
  jackknife_sample <- SUng.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(CUng.15$PEQ, size=nrow(SUng.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_SUng_CUng[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_SUng_CUng, breaks=20)
abline(v=t.all.p, col="red")

################################ Other Stem mammals vs. Crown mammals ###################### 

#import bin data
Bin14.data<-read.csv("Bin14.csv", header=T, row.names = 1)
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data
SMam.14 <-Bin14.data[ which(Bin14.data$Clade3=='Stem mammal'), ]
SMam.15 <-Bin15.data[ which(Bin15.data$Clade3=='Stem mammal'), ]
SMam.16 <-Bin16.data[ which(Bin16.data$Clade3=='Stem mammal'), ]
SMam.17 <-Bin17.data[ which(Bin17.data$Clade3=='Stem mammal'), ]

CMam.14 <-Bin14.data[ which(Bin14.data$Clade3=='Crown mammal'), ]
CMam.15 <-Bin15.data[ which(Bin15.data$Clade3=='Crown mammal'), ]
CMam.16 <-Bin16.data[ which(Bin16.data$Clade3=='Crown mammal'), ]
CMam.17 <-Bin17.data[ which(Bin17.data$Clade3=='Crown mammal'), ]

#check sample size in each category
nrow(SMam.14) #8
nrow(SMam.15) #18
nrow(SMam.16) #3
nrow(SMam.17) #3

nrow(CMam.14) #0 - not enough
nrow(CMam.15) #10
nrow(CMam.16) #6
nrow(CMam.17) #5

#ANOVA of original data Bin 15
t.all <- t.test(SMam.15$PEQ, CMam.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_15_SMam_CMam <- numeric(length = 1000*nrow(CMam.15))

for (i in 1:nrow(CMam.15)) {
  jackknife_sample <- CMam.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(SMam.15$PEQ, size=nrow(CMam.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_SMam_CMam[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_SMam_CMam, breaks=20)
abline(v=t.all.p, col="red")

#ANOVA of original data Bin 16
t.all <- t.test(SMam.16$PEQ, CMam.16$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_16_SMam_CMam <- numeric(length = 1000*nrow(SMam.16))

for (i in 1:nrow(SMam.16)) {
  jackknife_sample <- SMam.16[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(CMam.16$PEQ, size=nrow(SMam.16)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_16_SMam_CMam[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_16_SMam_CMam, breaks=20)
abline(v=t.all.p, col="red")

#ANOVA of original data Bin 17
t.all <- t.test(SMam.17$PEQ, CMam.17$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_17_SMam_CMam <- numeric(length = 1000*nrow(SMam.17))

for (i in 1:nrow(SMam.17)) {
  jackknife_sample <- SMam.17[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(CMam.17$PEQ, size=nrow(SMam.17)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_17_SMam_CMam[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_17_SMam_CMam, breaks=20)
abline(v=t.all.p, col="red")

###################### Euarchontoglires, Ferae, all Artio, all all Perisso ###################### 

#import bin data
Bin14.data<-read.csv("Bin14.csv", header=T, row.names = 1)
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data
Cet.14 <-Bin14.data[ which(Bin14.data$Clade4=='Cetartiodactyla'), ]
Cet.15 <-Bin15.data[ which(Bin15.data$Clade4=='Cetartiodactyla'), ]
Cet.16 <-Bin16.data[ which(Bin16.data$Clade4=='Cetartiodactyla'), ]
Cet.17 <-Bin17.data[ which(Bin17.data$Clade4=='Cetartiodactyla'), ]

Eur.14 <-Bin14.data[ which(Bin14.data$Clade4=='Euarchontoglires'), ]
Eur.15 <-Bin15.data[ which(Bin15.data$Clade4=='Euarchontoglires'), ]
Eur.16 <-Bin16.data[ which(Bin16.data$Clade4=='Euarchontoglires'), ]
Eur.17 <-Bin17.data[ which(Bin17.data$Clade4=='Euarchontoglires'), ]

Per.14 <-Bin14.data[ which(Bin14.data$Clade4=='Perrissodacyla'), ]
Per.15 <-Bin15.data[ which(Bin15.data$Clade4=='Perrissodacyla'), ]
Per.16 <-Bin16.data[ which(Bin16.data$Clade4=='Perrissodacyla'), ]
Per.17 <-Bin17.data[ which(Bin17.data$Clade4=='Perrissodacyla'), ]

Fer.14 <-Bin14.data[ which(Bin14.data$Clade4=='Ferae'), ]
Fer.15 <-Bin15.data[ which(Bin15.data$Clade4=='Ferae'), ]
Fer.16 <-Bin16.data[ which(Bin16.data$Clade4=='Ferae'), ]
Fer.17 <-Bin17.data[ which(Bin17.data$Clade4=='Ferae'), ]

#check sample size in each category
nrow(Cet.14) #7
nrow(Cet.15) #3
nrow(Cet.16) #19
nrow(Cet.17) #9

nrow(Eur.14) #2 - not enough
nrow(Eur.15) #11
nrow(Eur.16) #4
nrow(Eur.17) #2 - not enough

nrow(Per.14) #2 - not enough
nrow(Per.15) #12
nrow(Per.16) #3
nrow(Per.17) #1 - not enough

nrow(Fer.14) #0 - not enough
nrow(Fer.15) #9
nrow(Fer.16) #2 - not enough
nrow(Fer.17) #3 

##### All Cetartiodactyles vs. Euarchontoglires

nrow(Cet.15) #3
nrow(Eur.15) #11
nrow(Cet.16) #19
nrow(Eur.16) #4

#ANOVA of original data Bin 15
t.all <- t.test(Cet.15$PEQ, Eur.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_15_Cet_Eur <- numeric(length = 1000*nrow(Cet.15))

for (i in 1:nrow(Cet.15)) {
  jackknife_sample <- Cet.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Eur.15$PEQ, size=nrow(Cet.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_Cet_Eur[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_Cet_Eur, breaks=20)
abline(v=t.all.p, col="red")

#ANOVA of original data Bin 16
t.all <- t.test(Cet.16$PEQ, Eur.16$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_16_Cet_Eur <- numeric(length = 1000*nrow(Eur.16))

for (i in 1:nrow(Eur.16)) {
  jackknife_sample <- Eur.16[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Cet.16$PEQ, size=nrow(Eur.16)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_16_Cet_Eur[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_16_Cet_Eur, breaks=20)
abline(v=t.all.p, col="red")

##### All Cetartiodactyles vs. Ferae

nrow(Cet.15) #3
nrow(Fer.15) #9
nrow(Cet.17) #9
nrow(Fer.17) #3

#ANOVA of original data Bin 15
t.all <- t.test(Cet.15$PEQ, Fer.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_15_Cet_Fer <- numeric(length = 1000*nrow(Cet.15))

for (i in 1:nrow(Cet.15)) {
  jackknife_sample <- Cet.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Fer.15$PEQ, size=nrow(Cet.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_Cet_Fer[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_Cet_Fer, breaks=20)
abline(v=t.all.p, col="red")

#ANOVA of original data Bin 17
t.all <- t.test(Cet.17$PEQ, Fer.17$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_17_Cet_Fer <- numeric(length = 1000*nrow(Fer.17))

for (i in 1:nrow(Fer.17)) {
  jackknife_sample <- Fer.17[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Cet.17$PEQ, size=nrow(Fer.17)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_17_Cet_Fer[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_17_Cet_Fer, breaks=20)
abline(v=t.all.p, col="red")

##### Ferae vs. Euarchontoglires

nrow(Eur.15) #11
nrow(Fer.15) #9

#ANOVA of original data Bin 15
t.all <- t.test(Fer.15$PEQ, Eur.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_15_Fer_Eur <- numeric(length = 1000*nrow(Fer.15))

for (i in 1:nrow(Fer.15)) {
  jackknife_sample <- Fer.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Eur.15$PEQ, size=nrow(Fer.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_Fer_Eur[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_Fer_Eur, breaks=20)
abline(v=t.all.p, col="red")

##### All Perrissodactyles vs. Euarchontoglires

nrow(Eur.15) #11
nrow(Per.15) #12
nrow(Eur.16) #4
nrow(Per.16) #3

#ANOVA of original data Bin 15
t.all <- t.test(Per.15$PEQ, Eur.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_15_Per_Eur <- numeric(length = 1000*nrow(Eur.15))

for (i in 1:nrow(Eur.15)) {
  jackknife_sample <- Eur.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Per.15$PEQ, size=nrow(Eur.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_Per_Eur[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_Per_Eur, breaks=20)
abline(v=t.all.p, col="red")

#ANOVA of original data Bin 16
t.all <- t.test(Per.16$PEQ, Eur.16$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Not significantly different

#using the smallest subset
replicates_16_Per_Eur <- numeric(length = 1000*nrow(Per.16))

for (i in 1:nrow(Per.16)) {
  jackknife_sample <- Per.16[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Eur.16$PEQ, size=nrow(Per.16)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_16_Per_Eur[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_16_Per_Eur, breaks=20)
abline(v=t.all.p, col="red")

###### Ferae vs. All Perissodactyles

nrow(Per.15) #12
nrow(Fer.15) #9

#ANOVA of original data Bin 15
t.all <- t.test(Fer.15$PEQ, Per.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Significantly different

#using the smallest subset
replicates_15_Fer_Per <- numeric(length = 1000*nrow(Fer.15))

for (i in 1:nrow(Fer.15)) {
  jackknife_sample <- Fer.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(Per.15$PEQ, size=nrow(Fer.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_Fer_Per[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_Fer_Per, breaks=20)
abline(v=t.all.p, col="red")

###################### Euarchontoglires, Ferae, all Artio, all all Perisso ###################### 

#import bin data
Bin14.data<-read.csv("Bin14.csv", header=T, row.names = 1)
Bin15.data<-read.csv("Bin15.csv", header=T, row.names = 1)
Bin16.data<-read.csv("Bin16.csv", header=T, row.names = 1)
Bin17.data<-read.csv("Bin17.csv", header=T, row.names = 1)

#parse data
EurS.14 <-Bin14.data[ which(Bin14.data$Clade5B=='Stem Euarchontoglires'), ]
EurS.15 <-Bin15.data[ which(Bin15.data$Clade5B=='Stem Euarchontoglires'), ]
EurS.16 <-Bin16.data[ which(Bin16.data$Clade5B=='Stem Euarchontoglires'), ]
EurS.17 <-Bin17.data[ which(Bin17.data$Clade5B=='Stem Euarchontoglires'), ]

EurC.14 <-Bin14.data[ which(Bin14.data$Clade5B=='Crown Euarchontoglires'), ]
EurC.15 <-Bin15.data[ which(Bin15.data$Clade5B=='Crown Euarchontoglires'), ]
EurC.16 <-Bin16.data[ which(Bin16.data$Clade5B=='Crown Euarchontoglires'), ]
EurC.17 <-Bin17.data[ which(Bin17.data$Clade5B=='Crown Euarchontoglires'), ]

FerS.14 <-Bin14.data[ which(Bin14.data$Clade5B=='Stem Ferae'), ]
FerS.15 <-Bin15.data[ which(Bin15.data$Clade5B=='Stem Ferae'), ]
FerS.16 <-Bin16.data[ which(Bin16.data$Clade5B=='Stem Ferae'), ]
FerS.17 <-Bin17.data[ which(Bin17.data$Clade5B=='Stem Ferae'), ]

FerC.14 <-Bin14.data[ which(Bin14.data$Clade5B=='Crown Ferae'), ]
FerC.15 <-Bin15.data[ which(Bin15.data$Clade5B=='Crown Ferae'), ]
FerC.16 <-Bin16.data[ which(Bin16.data$Clade5B=='Crown Ferae'), ]
FerC.17 <-Bin17.data[ which(Bin17.data$Clade5B=='Crown Ferae'), ]

#check sample size in each category
nrow(EurS.14) #2
nrow(EurS.15) #3
nrow(EurS.16) #0 - not enough
nrow(EurS.17) #0 - not enough

nrow(EurC.14) #0 - not enough
nrow(EurC.15) #8
nrow(EurC.16) #4
nrow(EurC.17) #2 - not enough

nrow(FerS.14) #0 - not enough
nrow(FerS.15) #3
nrow(FerS.16) #1 - not enough
nrow(FerS.17) #2 - not enough

nrow(FerC.14) #0 - not enough
nrow(FerC.15) #5 
nrow(FerC.16) #1 - not enough
nrow(FerC.17) #1 - not enough

##### Crown Euarchontoglires vs. Stem Euarchontoglires

nrow(EurS.15) #3
nrow(EurC.15) #8

#ANOVA of original data Bin 15
t.all <- t.test(EurC.15$PEQ, EurS.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all### Yes significantly different

#using the smallest subset
replicates_15_EurC_EurS <- numeric(length = 1000*nrow(EurS.15))

for (i in 1:nrow(EurS.15)) {
  jackknife_sample <- EurS.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(EurC.15$PEQ, size=nrow(EurS.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_EurC_EurS[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_EurC_EurS, breaks=20)
abline(v=t.all.p, col="red")

###### Crown Ferae vs. Stem Ferae

nrow(FerC.15) #5 
nrow(FerS.15) #3

#ANOVA of original data Bin 15
t.all <- t.test(FerC.15$PEQ, FerS.15$PEQ, alternative="greater")
t.all.p <- t.all$p.value
t.all ### Not significantly different

#using the smallest subset
replicates_15_FerC_FerS <- numeric(length = 1000*nrow(FerS.15))

for (i in 1:nrow(FerS.15)) {
  jackknife_sample <- FerS.15[-i,]$PEQ
  for (j in 1:1000) {
    bootstrap_sample <- sample(FerC.15$PEQ, size=nrow(FerS.15)-1, replace=F)
    t.test.outcome <- t.test(jackknife_sample, bootstrap_sample, alternative="greater")
    replicates_15_FerC_FerS[i*j] <- t.test.outcome$p.value
  }
}

hist(replicates_15_FerC_FerS, breaks=20)
abline(v=t.all.p, col="red")

### END
