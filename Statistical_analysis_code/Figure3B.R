#####  Partial Mantel tests between the functional compositional variation (geographic distance-corrected) of methanogenic pathways and grouped factors. #####

rm(list=ls())

library(ggplot2)

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
C_gene_order <- C_annotation_data$Gene

data_label2site<- read.delim('Metagenome_Sample_List.txt',header = T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

generate_total_ab_table <- function(ri_ab_dir){
  list.ra <- list.files(ri_ab_dir,full.names = TRUE)
  files.ra <- list.ra[!file.info(list.ra)$isdir]
  file_list.ra <- lapply(files.ra, read.delim,header = T,row.names=1)
  file.list.ra <- list()
  for (i in 1:length(file_list.ra)){
    file.name <- basename(files.ra[i])
    gene <- unlist(strsplit(file.name,"_"))[1]
    data_file <- file_list.ra[[i]]
    file_list.ra[[i]] <- data.frame(sample=colnames(data_file),abundance=colSums(data_file),richness=sapply(data_file, function(x) sum(x != 0)))
    colnames(file_list.ra[[i]])[which(colnames(file_list.ra[[i]])=="abundance")] <- gene
    rownames(file_list.ra[[i]]) <- file_list.ra[[i]]$sample
    file.list.ra[[gene]] <- file_list.ra[[i]][,c(1:2)][-1]
  }
  merged_abundance <- data.frame()
  for (i in names(file.list.ra)){
    tmp <- file.list.ra[[i]]
    if(length(merged_abundance)==0){
      merged_abundance <- tmp
    }else{
      s.rownames <- rownames(merged_abundance)
      q.rownames <- rownames(tmp)
      tmp[,i] <- tmp[match(s.rownames,q.rownames),i]
      merged_abundance <- cbind(merged_abundance,tmp)
    }
  }
  
  return(merged_abundance)
}

ri_ab_dir_C <- "Methanogenic_Genes_Abundance_Table/"
merged_abundance_C <- generate_total_ab_table(ri_ab_dir_C)
merged_abundance_C <- merged_abundance_C[,C_gene_order]
merged_abundance_C <- merged_abundance_C[data_label2site$`Data lable`,]
rownames(merged_abundance_C) <- data_label2site$Sample
merged_abundance_C <- data.frame(t(merged_abundance_C))
merged_abundance_C <- merged_abundance_C[,order(colnames(merged_abundance_C))]

C_pathway <- data.frame(Gene=C_annotation_data$Gene,Pathway=C_annotation_data$Pathway)


library(dplyr)
library(linkET)
library(ggcorrplot)
library(vegan)
library(viridis)
library(ecodist)

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

category <- list(CLIM=c("T2M","WS2M","RH2M","PRECTOTCORR"),
                      Influent_ENV=c("I_TN","I_NH3_N","I_TS","I_SCOD","I_pH","I_Salinity","I_S_carbohydrate" ,"I_S_protein"),
                      Reactor_ENV=c("R_TN","R_NH3_N","R_pH","R_Salinity","R_T","R_TS","R_SCOD","R_S_carbohydrate","R_S_protein"),
                      Performance=c("MC","MP","AMP"))

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
for(p in unique(C_annotation_data$Pathway)[c(4,3,2,1)]){
  tmp <- merged_abundance_C[C_annotation_data[C_annotation_data$Pathway==p,]$Gene,]
  tmp <- as.data.frame(t(tmp))
  samp.tmp <- rownames(tmp)
  samp.env <- rownames(env)
  env_tmp <- env[match(samp.tmp,samp.env),]
  mp <- mantel.prepare.data(env_tmp,tmp,p,"TRUE",NULL)
  df.m <- rbind(df.m,mp)
}


# overall functional composition
tmp <- merged_abundance_C
tmp <- as.data.frame(t(tmp)) 
samp.tmp <- rownames(tmp)
samp.env <- rownames(env)
env_tmp <- env[match(samp.tmp,samp.env),]
mp <- mantel.prepare.data(env_tmp,tmp,"Overall methanogenic function","TRUE",NULL)
df.m <- rbind(df.m,mp)

library(ggplot2)
library(reshape2)
library(ggtext)

data <- df.m
data$significance <- rep("",nrow(data))
data[data$p>0.05,]$significance <- ""
data[data$p<=0.001,]$significance <- "***"
data[data$p<=0.01&data$p>0.001,]$significance <- "**"
data[data$p<=0.05&data$p>0.01,]$significance <- "*"
data$r <- as.numeric(data$r)
data$r2 <- data$r
data[data$p>0.05,]$r2 <- 0
data$env <- factor(data$env,levels = c("Geo_Distance","CLIM","Influent_ENV","Reactor_ENV","Performance"))
data$spec <- factor(data$spec,levels = rev(c("Overall methanogenic function","Central methanogenic pathway","Aceticlastic methanogenesis","Hydrogenotrophic methanogenesis","Methylotrophic methanogenesis")))

#plot
p <- ggplot(data, aes(x = env, y = spec)) +
  geom_square(aes(r0=abs(r),fill=abs(r)),color=NA)+
  geom_square(aes(x = env, y = spec,r0=abs(r2)), fill = "NA",colour = "#457033",size = 5, linetype = "dashed")+
  geom_text(aes(label = significance), size = 5)+
  scale_fill_gradient2(low = "white", high = "red", limit = c(0, 1),
                       space = "Lab")+
  theme_classic()+
  theme(plot.title = element_text(face = "bold", family = "sans"),
        axis.text = element_text(family = "serif",color = 'black',size = 10),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "#F5F5F5"),
        strip.text = element_text(face = "bold", family = "serif",size = 10),
        legend.text = element_text(face = "bold.italic", family = "serif"),
        legend.title = element_text(face = "bold.italic", family = "serif",hjust=0.5),
        panel.background = element_rect(colour = "black", linewidth = 0.5, linetype = "solid",fill = NULL),
        panel.grid = element_blank())+
  guides(fill = guide_colorbar(title = "Mantel's |r|"))

p
