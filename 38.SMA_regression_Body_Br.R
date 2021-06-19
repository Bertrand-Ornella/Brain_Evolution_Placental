############################################################################################
####################### Standard major axis/reduced major axis regression #########################
############################################################################################

setwd("~/Desktop/For github")

library(ggplot2)
library(smatr)

#Open data
data1<-read.csv("Rates_Body_Brain_negative_group_internal.csv", header=T, row.names = 1)

#Take internal branch out to calculate slopes
data<-subset(data1, data1$Group != "Internal")

#subset quadrants
Q1 <- data[which(data$Log.Median.Scalar.Body <= 0 & data$Log.Median.Scalar.Brain >= 0),]
Q2 <- data[which(data$Log.Median.Scalar.Body >= 0 & data$Log.Median.Scalar.Brain >= 0),]
Q3 <- data[which(data$Log.Median.Scalar.Body <= 0 & data$Log.Median.Scalar.Brain <= 0),]
Q4 <- data[which(data$Log.Median.Scalar.Body >= 0 & data$Log.Median.Scalar.Brain <= 0),]

#SMA regression on all data
sma_model <- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body-1, data = data, method = "SMA", alpha = 0.05)
sma_model
plot(sma_model)

#SMA regression with groups on all data
sma_model_grouped <- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body*Group-1, data = data, method = "SMA", alpha = 0.05)
sma_model_grouped
summary(sma_model_grouped)

#SMA regression on Quad 1
sma_model_quad1 <- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body-1, data = Q1, method = "SMA", alpha = 0.05)
sma_model_quad1

#Quad 1 with groups
sma_model_quad1_grouped <- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body*Group-1, data = Q1, method = "SMA", alpha = 0.05)
sma_model_quad1_grouped
summary(sma_model_quad1_grouped)

#SMA regression on Quad 2
sma_model_quad2 <- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body-1, data = Q2, method = "SMA", alpha = 0.05)
sma_model_quad2

#Quad 2 with groups
sma_model_quad2_grouped <- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body*Group-1, data = Q2, method = "SMA", alpha = 0.05)
sma_model_quad2_grouped
summary(sma_model_quad2_grouped)

sink("sma_model_quad2_grouped.txt")
print(summary(sma_model_quad2_grouped))
sink()

#SMA regression on Quad 3
sma_model_quad3 <- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body-1, data = Q3, method = "SMA", alpha = 0.05)
sma_model_quad3

#Quad 3 with groups
sma_model_quad3_grouped <- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body*Group-1, data = Q3, method = "SMA", alpha = 0.05)
sma_model_quad3_grouped
summary(sma_model_quad3_grouped)

#SMA regression on Quad 4
sma_model_quad4 <- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body-1, data = Q4, method = "SMA", alpha = 0.05)
sma_model_quad4

#Quad 4 with groups
sma_model_quad4_grouped <- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body*Group-1, data = Q4, method = "SMA", alpha = 0.05)
sma_model_quad4_grouped
summary(sma_model_quad4_grouped)

#Open data for graph - includes internal branches
data<-read.csv("Rates_Body_Brain_negative_group_internal.csv", header=T, row.names = 1)

#Data
Group<-data$Group

#Add new column to identify branches in graph
data$Branch<-seq(1:214)

#Rate
ggplot(data, aes(Log.Median.Scalar.Body, Log.Median.Scalar.Brain, color = Group)) +
  geom_point(size = 2)+
  labs (x = "Body mass rate", y = "Endocranial volume rate")+
  scale_shape_manual(values = c("Eocene crown groups","Eocene stem groups","Mesozoic","Paleocene", "Internal"))+
  scale_color_manual(values = c("blueviolet","goldenrod1","grey","orangered1","turquoise")) +
  geom_segment(aes(x = 0, xend = 10, y = 0, yend = 0.4029607*10),color="turquoise",size = 0.2) +
  geom_segment(aes(x = 0, xend = 10, y = 0, yend = 0.3900990*10),color="orangered1",size = 0.2) +
  geom_segment(aes(x = 0, xend = 10, y = 0, yend = 0.4429696*10),color="goldenrod1",size = 0.2) +
  geom_segment(aes(x = 0, xend = 10, y = 0, yend = 0.5626080*10),color="blueviolet",size = 0.2) +
  scale_x_continuous(expand=c(0,0), limits=c(-11,11)) +
  scale_y_continuous(expand=c(0,0), limits=c(-11,11)) +
  coord_cartesian(xlim=c(-11,11), ylim=c(-11,11)) +
  geom_hline(yintercept = 0,size = 0.2)+
  geom_vline(xintercept = 0, size = 0.2)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", size = 0.2) +
  #geom_text(aes(label=Branch),hjust=0, vjust=0)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))+
  theme_minimal()

####### Test for significant differences for slope ########

sink("sma_Q2_Pal.txt")
sma_Q2_pal<- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body*Group-1, 
              data = Q2, method = "SMA", alpha = 0.05, slope.test = 0.4029607)
summary(sma_Q2_pal)
sink()

sink("sma_Q2_Mes.txt")
sma_Q2_Mes<- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body*Group-1, 
                 data = Q2, method = "SMA", alpha = 0.05, slope.test = 0.3900990)
summary(sma_Q2_Mes)
sink()

sink("sma_Q2_EcS.txt")
sma_Q2_EcS<- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body*Group-1, 
                 data = Q2, method = "SMA", alpha = 0.05, slope.test = 0.4429696)
summary(sma_Q2_EcS)
sink()

sink("sma_Q2_EcC.txt")
sma_Q2_EcC<- sma(Log.Median.Scalar.Brain ~ Log.Median.Scalar.Body*Group-1, 
                 data = Q2, method = "SMA", alpha = 0.05, slope.test = 0.5626080)
summary(sma_Q2_EcC)
sink()

### END
