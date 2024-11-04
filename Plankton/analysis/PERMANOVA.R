#Load the packages

library(vegan)
library(tidyverse)
library(readxl)

#Import Data

set.seed(100)
data<-read_excel(file.choose())

perm_data<-data[,3:9]

#Transform

perm_data_trans<-sqrt(perm_data)

#Run NMDS Model for Visualizing the composition

nmds_result<-metaMDS (perm_data_trans, distance = "bray")


#Extract NMDS Scores 
nmds_scores <-as.data.frame(scores(nmds_result)$sites)

#Find out the centroids

group_centroids <- data.frame(
  Location = c("Lower", "Upper"),
  Centroid_X = c(mean(nmds_scores$NMDS1[data$Location == "Lower"]),
                 mean(nmds_scores$NMDS1[data$Location == "Upper"])),
  Centroid_Y = c(mean(nmds_scores$NMDS2[data$Location == "Lower"]),
                 mean(nmds_scores$NMDS2[data$Location == "Upper"]))
)


#Create data frame for ggplot

plot_data<-data.frame(
    Location = data$Location,
    NMDS1=nmds_scores$NMDS1,
    NMDS2=nmds_scores$NMDS2,
    xend=c(rep( group_centroids[1,2],10),rep(group_centroids[2,2],10)),
    yend=c(rep(group_centroids[1,3],10),rep( group_centroids[2,3],10)))


#Plot the data

ggplot(plot_data, aes(NMDS1,NMDS2)) + 
  geom_point(aes(color = Location),size=2)+ 
  stat_ellipse(geom = "polygon", alpha = 0.04, aes(group = Location), 
                             color = "black",fill="blue")+ 
  geom_point(data = group_centroids, aes(x = Centroid_X, y = Centroid_Y), 
             color = "black", size = 2, shape = 7)+
  geom_segment(data = plot_data, aes(x =NMDS1, y = NMDS2, 
              xend = xend, yend = yend, color = Location), alpha = 0.5)+
  scale_color_manual(name=" Locations",labels= unique(plot_data$Location),
                   
                                    values= c("darkgreen", "darkred"))+ theme_bw()
                     




##############PERMANOVA##############################


#Distance Matrix

perm_dist<-vegdist(perm_data_trans, method='bray')


#Assumptions

dispersion<-betadisper(perm_dist, group=data$Location,type = "centroid")

plot(dispersion)

anova(dispersion)


#Test

Perma_result<-adonis2( perm_dist~as.factor(plot_data$Location), data=perm_dist,
                       permutations=9999)

Perma_result

