# this script was used to analysis Lompolo vegetation cluster pattern in 2018
# the main purpose is mapping vegetation
library(tidyverse)  # data manipulation
library(tidyr)

# Load packages and functions
library(vegan) # load the vegan package to the work space
library(dplyr)
library(MASS)
library(ggplot2)
library(ggrepel)

library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(indicspecies)

library(ggforce) # annotate ellipses around clusters

# set the working directory
setwd("Z:/RS-Data/Veg_mapping/data")

# ordination and clustering
lom.sps <- read.csv("lompolo_2018_plot_species_v1.csv", header=T, sep =",", stringsAsFactors = T) #201&68
# v2:excludes non-veg components
lom.sps <- read.csv("lompolo_2018_plot_species_v2.csv", header=T, sep =",", stringsAsFactors = T) #201&62

names(lom.sps)

#######1 ordination#######
set.seed(123)
# Stress: 0.16
NMDS.lom.sps <- metaMDS(lom.sps[,-1], distance = "bray", k = 4, try = 20)

# plot(NMDS.lom.sps)
# text(NMDS.lom.sps, display = "species", cex=0.7, col="black")

df.lom.plot <- data.frame(lom.sps$Plot, NMDS.lom.sps$points) #
df.lom.spe <- data.frame(NMDS.lom.sps$species) #
write.csv(df.lom.plot,'lom18_plot_nmds_v2.csv')
write.csv(df.lom.spe,'lom18_species_nmds_v2.csv')

# Significant test,0.82;0.7978
lom.sig <- adonis2(lom.sps[,1]~ ., data=lom.sps, permutations = 999, method="bray", by = NULL)
#######1 ordination#######

#######2 clustering#######
fviz_nbclust(df.lom.plot, pam, method = "silhouette") #k=4
set.seed(123)
lom.fanny.plot <- fanny(df.lom.plot, k=4, memb.exp = 1.5, metric = "euclidean")

fviz_cluster(lom.fanny.sps, data = df.lom.plot[,-c(2,3)], geom = "point", 
             # ellipse.type = "norm",
             xlab = 'MDS3',ylab = 'MDS4',
             ggtheme = theme_bw())

da.lom.plot <- data.frame(df.lom.plot,'Cluster' = as.factor(lom.fanny.sps$clustering))
df.lom.spe

ggplot() + 
  geom_point(data = da.lom.plot, aes(x = MDS1, y = MDS2, color = Cluster), size = 2) +
  stat_ellipse(data = da.lom.plot, aes(x = MDS1, y = MDS2, color = Cluster), type = "norm", linetype = 2, lwd = 1.2) +
  # geom_point(df.lom.spe, aes(x = MDS1, y = MDS2), shape = 10, size = 1.5) +
  geom_text_repel(data = df.lom.spe, aes(x = MDS1, y = MDS2, label = row.names(df.lom.spe), color = 'black', fontface="italic")) +
  guides(color = 'none') +
  theme_bw(base_size = 15)


# # cluster species
# fviz_nbclust(df.lom.spe, pam, method = "silhouette") #k=5
# set.seed(123)
# lom.fanny.sps <- fanny(df.lom.spe, k=5, memb.exp = 1.5, metric = "euclidean")
# 
# fviz_cluster(lom.fanny.sps, df.lom.spe, 
#              # ellipse.type = "norm",
#              ggtheme = theme_bw())
# 
# da.lom.sps <- data.frame(df.lom.spe,'Cluster' = as.factor(lom.fanny.sps$clustering))
# 
# ggplot(da.lom.sps, aes(x = MDS3, y = MDS4, color = Cluster, label = rownames(da.lom.sps))) + 
#   geom_point(size = 2) +
#   geom_text(aes(fontface=3)) +
#   # stat_ellipse(type = "norm", linetype = 2, lwd = 1.2) +
#   # guides(color = 'none') +
#   theme_bw(base_size = 15)


table(lom.fanny.sps$clustering)
# 1  2  3  4
# 42 46 59 54

# export the clustering results
lom.cluster <- data.frame(df.lom.plot,lom.fanny.sps$clustering,lom.fanny.sps$membership)
write.csv(lom.cluster,'lom18_fcm_cluster.csv',row.names=FALSE)

lom.cluster <- data.frame(df.lom.plot,lom.fanny.sps$clustering,lom.fanny.sps$membership)
write.csv(lom.cluster,'lom18_fcm_cluster_v2.csv',row.names=FALSE)
#######2 clustering#######

#######3 indicator species#######
lom.clu <- read.csv("lom18_fcm_cluster.csv", header=T, sep =",", stringsAsFactors = T)
names(lom.clu)[c(1,6)] <- c('Plot','Cluster')
lom.clu$Cluster <- as.factor(lom.clu$Cluster)

lom.sps <- read.csv("lompolo_2018_plot_species_v1.csv", header=T, sep =",", stringsAsFactors = T) #201&68
lom.indsps <- multipatt(lom.sps[,-1], lom.clu$Cluster, control = how(nperm=999))
summary(lom.indsps, indvalcomp=TRUE, alpha = 0.05)

# Total number of species: 67
# Selected number of species: 20

# Only analysis the veg-component
lom.clu2 <- read.csv("lom18_fcm_cluster_v2.csv", header=T, sep =",", stringsAsFactors = T)
names(lom.clu2)[c(1,6)] <- c('Plot','Cluster')
lom.clu2$Cluster <- as.factor(lom.clu2$Cluster)

lom.sps2 <- read.csv("lompolo_2018_plot_species_v2.csv", header=T, sep =",", stringsAsFactors = T) #201&62
lom.indsps2 <- multipatt(lom.sps2[,-1], lom.clu2$Cluster, control = how(nperm=999))
summary(lom.indsps2, indvalcomp=TRUE, alpha = 0.05)
summary(lom.indsps2, indvalcomp=TRUE, alpha = 1) #all species list

#######3 indicator species#######

#######4 figures#######
####### all land cover #######
# nmds of all species
lom.spe <- read.csv("lom18_species_nmds.csv",header=T, row.names=1, sep =",") #67*4
# 20 indicator species
lom.indspe <- read.csv("lom18_indispe_nmds.csv", header=T, row.names=1, sep =",") #20*4

# fcm clustering
lom.clu <- read.csv("lom18_fcm_cluster.csv", header=T, sep =",", stringsAsFactors = T)
names(lom.clu)[c(1,6)] <- c('Plot','Cluster')
lom.clu$Cluster <- as.factor(lom.clu$Cluster)

# plot NMDS
# 1 Plotting whole species
library(ggrepel)
p1 <- ggplot(data = lom.clu) + 
  geom_point(aes(x = MDS1, y = MDS2, shape = Cluster), size = 2) +
  scale_shape_manual(values=c(15:18)) +
  theme_bw(base_size = 14)
  # geom_mark_ellipse(aes(x = MDS1, y = MDS2, fill = Cluster)) +

p1 + geom_point(data = lom.spe, aes(x = MDS1, y = MDS2), shape = 10, size = 1.5) +
  geom_text_repel(data = lom.spe, aes(x = MDS1, y = MDS2, label = row.names(lom.spe), 
                                      color = 'red', fontface="italic")) + 
  theme_bw(base_size = 14)
                     
  # geom_point(data = lom.spe, aes(x = MDS1, y = MDS2), shape = 42, size = 3) +
  # geom_text(data = lom.spe, aes(x = MDS1, y = MDS2, label = row.names(lom.spe)),
  #           check_overlap = T, vjust = "inward", hjust = 0, color = 'red', fontface="italic") +
  

# 2 Only Plotting indicator species
# Change the order of clusters
lom.clu$cluster <- factor(lom.clu$cluster, levels=c(1,2,3,4),
                          labels=c("String margin", "Riparian fen", "Flark", "wet flark"))

ggplot(data = lom.clu) + 
  geom_point(aes(x = MDS1, y = MDS2, shape = Cluster), size = 2) +
  scale_shape_manual(values=c(15:18)) +
  # geom_mark_ellipse(aes(x = MDS1, y = MDS2, fill = Cluster)) +
  geom_point(data = lom.indspe, aes(x = MDS1, y = MDS2), color = 'red', shape = 10, size = 1.5) + 
  geom_text_repel(data = lom.indspe, aes(x = MDS1, y = MDS2, label = row.names(lom.indspe)), 
                                      color = 'red', fontface="italic") + 
  ggtitle("Lompolojankka_2018") +
  guides(color = 'none') +
  theme_bw(base_size = 14)
####### all land cover #######

####### veg-component #######
lom.spe2 <- read.csv("lom18_species_nmds_v2.csv",header=T, row.names=1, sep =",") #61*4
# 20 indicator species
lom.indspe2 <- read.csv("lom18_indispe_nmds_v2.csv", header=T, row.names=1, sep =",") #18*4
lom.indspe2$Cluster <- as.factor(lom.indspe2$Cluster)

# fcm clustering
lom.clu2 <- read.csv("lom18_fcm_cluster_v2.csv", header=T, sep =",", stringsAsFactors = T)
names(lom.clu2)[c(1,6)] <- c('Plot','Cluster')
lom.clu2$Cluster <- as.factor(lom.clu2$Cluster)

# Plot
ggplot(data = lom.clu2) + 
  geom_point(aes(x = MDS3, y = MDS4, color = Cluster), size = 2) +
  # scale_shape_manual(values=c(15:18)) +
  stat_ellipse(aes(x = MDS3, y = MDS4, color = Cluster), type = "norm", linetype = 2, lwd = 1.2) +
  geom_point(data = lom.indspe2, aes(x = MDS3, y = MDS4), color = 'black', shape = 10, size = 1.5) + 
  geom_text_repel(data = lom.indspe2, aes(x = MDS3, y = MDS4, label = row.names(lom.indspe2), fontface="italic")) + 
  # facet_wrap(~ Cluster) +
  # guides(color = 'none') +
  theme_bw(base_size = 14)

ggplot(lom.indspe2, aes(x = MDS1, y = MDS2, color = Cluster)) + 
  geom_point(shape = 10, size = 1.5) + 
  geom_text_repel(aes(label = row.names(lom.indspe2), fontface="italic")) + 
  # stat_ellipse(aes(x = MDS1, y = MDS2, color = Cluster), type = "norm", linetype = 2, lwd = 1.2) +
  # guides(color = 'none') +
  theme_bw(base_size = 14)

##### Version 2, more indicator species #####
# 46 indicator species
lom.indspe3 <- read.csv("lom18_indispe_nmds_v3.csv", header=T, row.names=1, sep =",") #25*4
lom.indspe3$Cluster <- as.factor(lom.indspe3$Cluster)

# fcm clustering
lom.clu2 <- read.csv("lom18_fcm_cluster_v2.csv", header=T, sep =",", stringsAsFactors = T)
names(lom.clu2)[c(1,6)] <- c('Plot','Cluster')
lom.clu2$Cluster <- as.factor(lom.clu2$Cluster)

# Plot
ggplot(data = lom.clu2) + 
  # geom_point(aes(x = MDS3, y = MDS4, color = Cluster), size = 2) +
  # scale_shape_manual(values=c(15:18)) +
  # stat_ellipse(aes(x = MDS3, y = MDS4, color = Cluster), type = "norm", linetype = 2, lwd = 1.2) +
  geom_point(data = lom.indspe3, aes(x = MDS1, y = MDS2, color = Cluster), shape = 10, size = 1.5) + 
  geom_text_repel(data = lom.indspe3, aes(x = MDS1, y = MDS2, label = row.names(lom.indspe3),color = Cluster), fontface="italic") + 
  # facet_wrap(~ Cluster) +
  # guides(color = 'none') +
  theme_bw(base_size = 14)

ggplot(lom.indspe3, aes(x = MDS1, y = MDS2, color = Cluster)) + 
  geom_point(shape = 10, size = 1.5) + 
  geom_text_repel(aes(label = row.names(lom.indspe3)),fontface="italic") + 
  # facet_wrap(~ Cluster) +
  # stat_ellipse(aes(x = MDS1, y = MDS2, color = Cluster), type = "norm", linetype = 2, lwd = 1.2) +
  # guides(color = 'none') +
  theme_bw(base_size = 14)

####### veg-component #######