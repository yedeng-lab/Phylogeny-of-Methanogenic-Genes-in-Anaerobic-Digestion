##### Partial Mantel tests between the phylogenetic variation (geographic distance-corrected) of each functional gene and different categories of factors #####

rm(list=ls())

plot.theme = theme(plot.title=element_text(size=15, color="black", family  = "serif", face= "bold",vjust=0.5,hjust=0),
                   axis.line=element_line(linewidth = 0.5, colour="black"),
                   axis.ticks = element_line(linewidth = 1),
                   axis.ticks.length = unit(0.3,"cm"),
                   axis.text= element_text(size=18, color="black", family  = "serif", face= "bold", vjust=0.5, hjust=0.5),
                   axis.title = element_text(size=18, color="black", family  = "serif",face= "bold", vjust=0.5, hjust=0.5),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text = element_text(face = "bold.italic", family = "serif",size = 18),
                   legend.title = element_text(face = "bold.italic", family = "serif",hjust=0.5,size = 18),
                   panel.background = element_rect(colour = "black", linewidth = 0.5, linetype = "solid",fill = NULL),
                   panel.grid=element_blank())

setwd("/AD_Methane_Data")

C_annotation_data<- read.delim('Methanogenic_Genes_Annotation.txt',header = T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
C_annotation_data[C_annotation_data$Pathway=="Hydrogenotrophic methanogenesis",]$Pathway <- "Hydrogenotrophic\nmethanogenesis"
C_annotation_data[C_annotation_data$Pathway=="Aceticlastic methanogenesis",]$Pathway <- "Aceticlastic\nmethanogenesis"
C_annotation_data[C_annotation_data$Pathway=="Methylotrophic methanogenesis",]$Pathway <- "Methylotrophic\nmethanogenesis"
C_gene_order <- C_annotation_data$Gene
C_pathway <- data.frame(Gene=C_annotation_data$Gene,Pathway=C_annotation_data$Pathway)
pathway_order <- unique(C_annotation_data$Pathway) 

data_label2site<- read.delim('Metagenome_Sample_List.txt',header = T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

otu_dir_C <- "Methanogenic_Genes_Abundance_Table/"
file.list <- list()
for (g in C_gene_order){
  data_file_path <- paste0(otu_dir_C,"/",g,"_normalized_abun_int.txt")
  tmp <- read.delim(data_file_path, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,row.names = 1)
  tmp <- as.data.frame(t(tmp))
  rownames(tmp) <- data_label2site[match(rownames(tmp),data_label2site$`Data lable`),]$Sample
  file.list[[g]] <- tmp
}

#####read env+isc data#####
#filtered env
env_filtered <-  read.delim('Environmental_Factors.txt', row.names = 1)
#filter na
env_filtered2 <- env_filtered%>%select_if(~!any(is.na(.)))
isc <-  read.delim('Performance_Indicators.txt', row.names = 1)
selected.isc <- c("MC","MP","AMP")
CLIM <- read.delim('Climatic_Factors.txt', row.names = 1)
selected.clim <- c("T2M", "WS2M", "RH2M", "PRECTOTCORR")
env <- cbind(CLIM[selected.clim],env_filtered2,isc[selected.isc])
rownames(env) <- data_label2site[match(rownames(env),data_label2site$Sample),]$`Data lable` 

category <- list(CLIM=c("T2M","WS2M","RH2M","PRECTOTCORR"),
                 Influent_ENV=c("I_TN","I_NH3_N","I_TS","I_SCOD","I_pH","I_Salinity","I_S_carbohydrate" ,"I_S_protein"),
                 Reactor_ENV=c("R_TN","R_NH3_N","R_pH","R_Salinity","R_T","R_TS","R_SCOD","R_S_carbohydrate","R_S_protein"),
                 Performance=c("MC","MP","AMP"))

# uw-Unifrac
un_dir_C <- "UniFrac_Matrix/"
file.list <- list()
for (g in C_gene_order){
  data_file_path <- paste0(un_dir_C,"/",g,"_Unweighted_Unifrac_matrix.txt")
  tmp <- read.delim(data_file_path, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,row.names = 1)
  file.list[[g]] <- tmp
}

mantel.prepare.data <- function(env_data,taxon,taxon_name,partial,distance_matrix,sorensen=F){
  if(is.null(taxon)){
    taxon.dis = distance_matrix
  }else{
    taxon.dis = vegdist(taxon, method="bray",binary=sorensen)
  }
  geo <- env_data[c('Longitude', 'Latitude')]
  env_data <- decostand(env_data,method = "standardize",MARGIN = 2)
  env_data2 <- env_data[,!(names(env_data) %in% c("Longitude", "Latitude"))]
  sitedis <- geosphere::distm(geo[c('Longitude', 'Latitude')])
  sitedis[sitedis==0] = 1 
  sitedis = as.dist(log10(sitedis))
  envdis <- vegdist(env_data2,method = "euclidean", na.rm=T)
  m <- ecodist::mantel(taxon.dis~sitedis)
  m.p <- ecodist::mantel(taxon.dis~sitedis+envdis)
  mantel <- data.frame(env="Geo_Distance",r=m["mantelr"],p=m["pval1"])
  mantel.p <- data.frame(env="Geo_Distance",r=m.p["mantelr"],p=m.p["pval1"])
  for(i in names(category)){
    ii <- category[[i]]
    envdis <- vegdist(env_data[,ii],method = "euclidean", na.rm=T)
    envdis.ctrl <- vegdist(env_data[,!(colnames(env_data)%in%ii)],method = "euclidean", na.rm=T)
    m <- ecodist::mantel(taxon.dis~envdis)
    m.p <- ecodist::mantel(taxon.dis~envdis+sitedis)
    mantel <- rbind(mantel,c(i,m["mantelr"],m["pval1"]))
    mantel.p <- rbind(mantel.p,c(i,m.p["mantelr"],m.p["pval1"]))
  }
  rownames(mantel) <- c(1:nrow(mantel))
  rownames(mantel.p) <- c(1:nrow(mantel.p))
  mantel <- cbind(spec=rep(taxon_name,nrow(mantel)),mantel)
  mantel.p <- cbind(spec=rep(taxon_name,nrow(mantel.p)),mantel.p)
  if(partial=="TRUE"){return(mantel.p)}
  else{return(mantel)}
}

df.m <- data.frame()
for(g in C_gene_order){
  tmp <- file.list[[g]]
  samp.tmp <- rownames(tmp)
  samp.env <- rownames(env)
  env_tmp <- env[match(samp.tmp,samp.env),]
  u.dist <- as.dist(tmp)
  mp <- mantel.prepare.data(env_tmp,NULL,g,"TRUE",u.dist)
  df.m <- rbind(df.m,mp)
}

library(ggplot2)
library(reshape2)
library(ggtext)

data <- df.m
data$Significance <- rep("",nrow(data))
data[data$p>=0.05,]$Significance <- "NS"
data[data$p<0.05,]$Significance <- "P < 0.05"
data$r <- as.numeric(data$r)
data$r2 <- data$r
data[data$p>=0.05,]$r2 <- 0
data$env <- factor(data$env,levels = c("Geo_Distance","CLIM","Influent_ENV","Reactor_ENV","Performance"))
data$spec <- factor(data$spec,levels = rev(C_gene_order))
data$Pathway <- C_pathway[match(data$spec,C_pathway$Gene),]$Pathway
data$Pathway <- factor(data$Pathway,levels = pathway_order)
data <- data[!data$env=="Geo_Distance",]

#plot
p <- ggplot(data, aes(x = env, y = spec)) +
  geom_tile(aes(fill = Significance),color = "gray")+
  scale_fill_manual(values=c("lightgrey","blue"))+
  geom_point(aes(color=abs(r)),size=4)+
  facet_grid(rows = vars(Pathway), scales = "free",space = "free", switch = "y")+
  scale_color_gradient2(low = "white", high = "red", limit = c(0, max(data$r)),
                        space = "Lab")+
  scale_y_discrete(position = "right")+
  scale_x_discrete(position = "bottom")+
  theme_classic()+
  theme(plot.title = element_text(face = "bold", family = "sans"),
        axis.text = element_text(size = 20,family = "serif",color = 'black'),
        axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y.right = element_text(face = "italic"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "#F5F5F5"),
        strip.text = element_text(face = "bold", family = "serif",size = 20),
        # strip.placement = "outside", 
        legend.text = element_text(face = "bold.italic", family = "serif"),
        legend.title = element_text(face = "bold.italic", family = "serif",hjust=0.5),
        panel.background = element_rect(colour = "black", linewidth = 0.5, linetype = "solid",fill = NULL),
        panel.grid = element_blank())+
  guides(color = guide_colorbar(title = "Mantel's |r|"),fill="none")

p

# stat 
stat <- data%>%dplyr::group_by(env)%>%summarise(`P < 0.05`=sum(p<0.05),NS=sum(p>=0.05))
stat$env <- factor(stat$env,levels = c("CLIM","Influent_ENV","Reactor_ENV","Performance"))
stat.longer <- stat%>%pivot_longer(c(`P < 0.05`,NS),names_to = "Mantel's P",values_to = "Count")
stat.longer$`Mantel's P` <- factor(stat.longer$`Mantel's P`,levels = c("NS","P < 0.05"))
stat.longer$Percent <- stat.longer$Count/114*100

p2 <- ggplot(stat.longer, aes(x = env, y = Percent, fill = `Mantel's P`)) +
  geom_bar(width = 0.8,stat = "identity",position='stack') +
  scale_x_discrete(position = "bottom") +
  scale_color_manual(values=c("lightgrey","blue"))+
  scale_fill_manual(values=c("lightgrey","blue"))+
  theme_classic()+theme( axis.ticks.x = element_blank(),axis.line.x  = element_blank(),
                         strip.text = element_blank(),strip.background = element_blank(), 
                         axis.text.x = element_blank(),axis.title.x = element_blank(),
                         axis.title.y = element_text(size = 20, face = "bold", family = "serif"),
                         axis.text.y = element_text(size = 18,family = "serif",face = "bold"),
                         legend.text = element_text(family = "serif",face = "bold"),
                         legend.title = element_text( face = "bold.italic", family = "serif",hjust=0.5)
  )+ylab("Percentage (%)")
p2


library(patchwork)
(p3 <- p2/p+plot_layout(heights = c(1,14),guides = "collect"))
