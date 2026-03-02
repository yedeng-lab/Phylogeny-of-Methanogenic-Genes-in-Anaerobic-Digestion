##### Taxonomic and phylogenetic β-diversity of AD microbiome #####
rm(list =ls())
library(ggplot2)
library(ade4) 
library(vegan)
library(RColorBrewer)
library(plyr)

setwd("/AD_Methane_Data")

plot.theme = theme(plot.title=element_text(size=20, color="black", family  = "serif", face= "bold",vjust=0.5,hjust=0),
                   axis.line=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text= element_text(size=20, color="black", family  = "serif", face= "bold", vjust=0.5, hjust=0.5),
                   axis.title = element_text(size=20, color="black", family  = "serif",face= "bold", vjust=0.5, hjust=0.5),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.title = element_text( face = "bold.italic", family = "serif",hjust=0.5,size=20),
                   legend.text = element_text(colour = 'black', size = 18,  family  = "serif",face = 'bold'),
                   panel.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank())

##16S
##BC-d
otu<- read.delim('AD_OTU_Table.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,row.names = 1)
otu <- t(otu)
otu <- as.data.frame(otu)
Grouplist<- read.delim('Metagenome_Sample_List.txt',header = T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
samp.otu <- rownames(otu)
samp.group.list <- Grouplist$Sample
matched <- match(samp.otu,samp.group.list)
Grouplist <- Grouplist[matched,]
KO.dist = vegdist(otu,method='bray')
data_pcoa <- cmdscale(KO.dist,k=(nrow(otu)-1),eig=TRUE) 
library(ape)
bray_pcoa <- pcoa(KO.dist)
taxonomic_pcoa1 <- bray_pcoa$vectors[,1]
data_pcoa_eig <- data_pcoa$eig 
data_pcoa_exp <- data_pcoa_eig/sum(data_pcoa_eig)  
pcoa1 <- paste(round(100*data_pcoa_exp[1],2),'%')
pcoa2 <- paste(round(100*data_pcoa_exp[2],2),'%')
sample_site <- data.frame(data_pcoa$points)[1:2]
sample_site$level<-factor(Grouplist[,"Site"],levels = c('CS','FS','WZ','JZ','BJ','QQ','QH'))
names(sample_site)[1:3] <- c('PCoA1', 'PCoA2','level')
head(sample_site)
summary(sample_site)
level_order <- c('CS','FS','WZ','JZ','BJ','QQ','QH')
level_order <- factor(1:length(level_order),labels = level_order)
sample_site$level <- factor(sample_site$level,levels = levels(level_order))
PCoA1_mean <- tapply(sample_site$PCoA1,sample_site$level, mean)
PCoA2_mean <- tapply(sample_site$PCoA2,sample_site$level, mean)
mean_point <- rbind(PCoA1_mean,PCoA2_mean)
mean_point <- data.frame(mean_point)
t1 <- t(mean_point)
t2 <- as.data.frame(t1)
t2
sample_site_reordered <- sample_site[order(sample_site$level),]
colnames(sample_site_reordered)[3] <- "Site"
pcoa_plot <- ggplot(sample_site_reordered, aes(PCoA1, PCoA2, fill=Site,color=Site,shape=Site)) +
  plot.theme+        
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  stat_ellipse(aes(color = Site,fill = Site),type = "t", geom ="polygon", alpha=0.1, level = 0.95, show.legend = FALSE, linetype = 2) +
  geom_point(size = 4,aes(color = Site)) +
  scale_shape_manual(values = c(20,20,20,20,20,20,20)) +
  scale_color_manual(values = c("#38A3A5" ,"#EFC86E","#97C684", "#43C59E", "#5862A0" , "#1F6E9C" ,"#808FE1")) +
  scale_fill_manual(values = c("#A4DFE0" ,"#F3D591","#B8D8AB", "#A1E2CF", "#B0B5D4" , "#99CDEB" ,"#ACB6EC")) +
  theme(panel.grid = element_line(color = 'black', linetype = 1, size = 0.1),
        panel.background = element_rect(color = 'black', fill = 'transparent')
  )+labs(x = paste("PCoA1: ",pcoa1,sep=""), y = paste("PCoA2: ",pcoa2,sep=""),title = "Bray-Curtis distance")
pcoa_plot


#Unifrac
rm(list = ls()[-which(ls()[]=="plot.theme")])
Unifrac<- read.delim('UniFrac_Matrix/Community_Unweighted_Unifrac_matrix.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,row.names = 1)
u.dist <- as.dist(Unifrac)
Grouplist<- read.delim('Metagenome_Sample_List.txt',header = T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
samp.otu <- rownames(Unifrac)
samp.group.list <- Grouplist$Sample
matched <- match(samp.otu,samp.group.list)
Grouplist <- Grouplist[matched,]
KO.dist = u.dist
data_pcoa <- cmdscale(KO.dist,k=(nrow(Unifrac)-1),eig=TRUE) 
library(ape)
bray_pcoa <- pcoa(KO.dist)
phylogentic_pcoa1 <- bray_pcoa$vectors[,1]
data_pcoa_eig <- data_pcoa$eig 
data_pcoa_exp <- data_pcoa_eig/sum(data_pcoa_eig)  
pcoa1 <- paste(round(100*data_pcoa_exp[1],2),'%')
pcoa2 <- paste(round(100*data_pcoa_exp[2],2),'%')
sample_site <- data.frame(data_pcoa$points)[1:2]
sample_site$level<-factor(Grouplist[,"Site"],levels = c('CS','FS','WZ','JZ','BJ','QQ','QH'))
names(sample_site)[1:3] <- c('PCoA1', 'PCoA2','level')
head(sample_site)
summary(sample_site)
level_order <- c('CS','FS','WZ','JZ','BJ','QQ','QH')
level_order <- factor(1:length(level_order),labels = level_order)
sample_site$level <- factor(sample_site$level,levels = levels(level_order))
PCoA1_mean <- tapply(sample_site$PCoA1,sample_site$level, mean)
PCoA2_mean <- tapply(sample_site$PCoA2,sample_site$level, mean)
mean_point <- rbind(PCoA1_mean,PCoA2_mean)
mean_point <- data.frame(mean_point)
t1 <- t(mean_point)
t2 <- as.data.frame(t1)
t2
sample_site_reordered <- sample_site[order(sample_site$level),]
colnames(sample_site_reordered)[3] <- "Site"
pcoa_plot <- ggplot(sample_site_reordered, aes(PCoA1, PCoA2, fill=Site,color=Site,shape=Site)) +
  plot.theme+        
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  stat_ellipse(aes(color = Site,fill = Site),type = "t", geom ="polygon", alpha=0.1, level = 0.95, show.legend = FALSE, linetype = 2) +
  geom_point(size = 4,aes(color = Site)) +
  scale_shape_manual(values = c(20,20,20,20,20,20,20)) +
  scale_color_manual(values = c("#38A3A5" ,"#EFC86E","#97C684", "#43C59E", "#5862A0" , "#1F6E9C" ,"#808FE1")) +
  scale_fill_manual(values = c("#A4DFE0" ,"#F3D591","#B8D8AB", "#A1E2CF", "#B0B5D4" , "#99CDEB" ,"#ACB6EC")) +
  theme(panel.grid = element_line(color = 'black', linetype = 1, size = 0.1),
        panel.background = element_rect(color = 'black', fill = 'transparent')
  )+
  labs(x = paste("PCoA1: ",pcoa1,sep=""), y = paste("PCoA2: ",pcoa2,sep=""),title = "Unweighted Unifrac distance")

pcoa_plot
