#####  Partial Mantel tests between the functional compositional variation (geographic distance-corrected) of methanogenic pathways and individual factors. #####
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
library(ggplot2)
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

mantel.prepare.data <- function(env_data,taxon,taxon_name,partial,distance_matrix,sorensen=F){
  if(is.null(taxon)){
    taxon.dis = distance_matrix
  }else{
    taxon.dis = vegdist(taxon, method="bray",binary=sorensen)
  }
  geo <- env_data[c('Longitude', 'Latitude')]
  env_data <- decostand(env_data,method = "standardize",MARGIN = 2)
  env_data2 <- env_data[,!(names(env_data) %in% c("Longitude"))]
  sitedis <- geosphere::distm(geo[c('Longitude', 'Latitude')])
  sitedis[sitedis==0] = 1 
  sitedis = as.dist(log10(sitedis))
  envdis <- vegdist(env_data2,method = "euclidean", na.rm=T)
  m <- ecodist::mantel(taxon.dis~sitedis)
  m.p <- ecodist::mantel(taxon.dis~sitedis+envdis)
  mantel <- data.frame(env="Geo_Distance",r=m["mantelr"],p=m["pval1"])
  mantel.p <- data.frame(env="Geo_Distance",r=m.p["mantelr"],p=m.p["pval1"])
  for(i in colnames(env_data2)){
    envdis <- vegdist(env_data[,i],method = "euclidean", na.rm=T)
    envdis.ctrl <- vegdist(env_data[,-which(colnames(env_data)==i)],method = "euclidean", na.rm=T)
    m <- ecodist::mantel(taxon.dis~envdis)
    m.p <- ecodist::mantel(taxon.dis~envdis+sitedis+envdis.ctrl)
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

# mantel <- df
mantel_plot <- function(mantel,env){
  for(i in 3:ncol(mantel)){mantel[,i] <- as.numeric(mantel[,i])}
  mantel <- mutate(mantel,
                   rd = cut(r, breaks = c(-Inf, 0.1, 0.2, Inf),
                            labels = c("< 0.1", "0.1 - 0.2", ">= 0.2")),
                   pd = cut(p, breaks = c(-Inf, 0.05, Inf),
                            labels = c("< 0.05",">= 0.05")))
  
  env2 <- na.omit(env)
  correlate(env2,method = "pearson",engine = "Hmisc")
  #筛选p值小于等于0.05
  cor_p <- Hmisc::rcorr(as.matrix(env),type = "pearson")
  cor_p[[1]]#correlation
  cor_p[[3]]#p value
  
  #if p value>0.05, the cor is null
  cor_p2 <- cor_p[[1]] %>%
    reshape2::melt() %>%
    as_tibble() %>%
    left_join(cor_p[[3]] %>%
                reshape2::melt() %>%
                as_tibble(),
              by = c("Var1","Var2")) %>%
    dplyr::mutate(value.x = as.character(value.x),
                  value.x = if_else(value.y > 0.05,"",value.x),
                  value.x = as.numeric(value.x))
  
  cor_p3 <- cor_p2 %>%
    select(Var1,Var2,value.x) %>%
    reshape2::acast(Var1 ~ Var2,value.var = "value.x") %>%
    data.frame()
  
  cor_p3[is.na(cor_p3)]=0
  #绘图
  p <- qcorrplot(cor_p3, type = "lower", diag = FALSE) +
    geom_square() +
    geom_couple(aes(colour = pd, size = rd), 
                data = mantel,
                label.family = "serif",
                label.fontface = "bold.italic",
                node.colour = c("#457033", "grey"),
                node.fill = c("#97C684", "grey"),
                node.size = c(3.5, 2),
                curvature = nice_curvature()) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(min(cor_p3), max(cor_p3)),
                         space = "Lab") +
    scale_size_manual(values = c(1,2,3)) +
    scale_colour_manual(values = c("#256947","#B0E1C3")) +
    guides(size = guide_legend(title = "Mantel's r",
                               override.aes = list(colour = "grey35"), 
                               order = 2),
           colour = guide_legend(title = "Mantel's p", 
                                 override.aes = list(size = 3), 
                                 order = 1),
           fill = guide_colorbar(title = "Pearson's r", order = 3))+ 
    theme(
      axis.text = element_text(family = "serif",face = "bold.italic",size = 10),axis.text.x.top = element_text(angle = 0,hjust = 0.5,vjust=0),
      legend.text = element_text(face = "bold.italic",family = "serif",size=10), legend.title = element_text( face = "bold.italic", family = "serif",hjust=0.5,size=10)
    )
  return(p) 
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

#plot
df.m <- df.m[df.m$env!="Geo_Distance",]
# p <- mantel_plot(df,subset(env,select = -c(Latitude,Longtitude)))
p <- mantel_plot(df.m,env[!(names(env) %in% c("Longitude"))])
p
