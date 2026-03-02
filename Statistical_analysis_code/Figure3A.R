##### Spearman correlations between the abundance of methanogenic pathways and environmental factors #####
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

#####read env+isc data#####
group1_y <-c("T2M", "WS2M", "RH2M", "PRECTOTCORR")
group2_y <- c("I_TN","I_NH3_N","I_pH","I_Salinity","I_TS","I_SCOD","I_S_carbohydrate","I_S_protein")
group3_y <- c("R_TN","R_NH3_N","R_pH","R_Salinity","R_T","R_TS","R_SCOD","R_S_carbohydrate","R_S_protein")
#filtered env
env_filtered <-  read.delim('Environmental_Factors.txt', row.names = 1)
env_filtered <- env_filtered[, !(names(env_filtered) %in% c("Longitude", "Latitude"))]
#filter na
env_filtered2 <- env_filtered%>%select_if(~!any(is.na(.)))
CLIM <- read.delim('Climatic_Factors.txt', row.names = 1)
selected.clim <- c("T2M", "WS2M", "RH2M", "PRECTOTCORR")
env <- cbind(CLIM[selected.clim],env_filtered2)

library(ggplot2)
library(cowplot)
library(corrplot)
library(tidyr)
library(MetBrewer)
library(linkET)
library(dplyr)


calculate_cor_matrix <- function(community_data,env_data){
  M <- list()
  M.p <- data.frame(matrix(ncol = ncol(env_data), nrow = ncol(community_data)))
  env_data <- env_data[match(rownames(community_data),rownames(env_data)),]
  rownames(M.p) <- colnames(community_data)
  colnames(M.p) <- colnames(env_data)
  M.r <- M.p
  for(i in colnames(community_data)){
    tmp <- community_data[[i]]
    if(sum(!is.na(tmp))==0){
      M.p[i,j] <- NA
      M.r[i,j] <- NA
      next
    }
    for (j in colnames(env_data)){
      cor <- cor.test(tmp,env_data[,j],alternative="two.sided",method="spearman",exact = F)
      M.p[i,j] <- cor$p.value
      M.r[i,j] <- cor$estimate
    }
  }
  M[["p"]] <- M.p
  M[["r"]] <- M.r
  return(M)
}

p_adjust <- function(cor.matrix.ri){
  M.p <- cor.matrix.ri[["p"]]
  M.r <- cor.matrix.ri[["r"]]
  p_vector <- as.vector(as.matrix(M.p))
  adjusted_p_vector <- p.adjust(p_vector, method = "fdr")
  adjusted_p_values <- matrix(adjusted_p_vector, nrow = nrow(M.p), byrow = FALSE)
  rownames(adjusted_p_values) <- rownames(M.p)
  colnames(adjusted_p_values) <- colnames(M.p)
  adjusted_p_values <- as.data.frame(adjusted_p_values)
  M.p2 <- adjusted_p_values
  print(sum(M.p<0.05))
  print(sum(M.p2<0.05))
  M.r2 <- M.r
  for (i in colnames(M.p2)){
    for (j in rownames(M.p2)){
      if (M.p2[j,i]>=0.05){
        M.r2[j,i] <- 0
      }
    }
  }
  cor.matrix.ri[["p"]] <- M.p2
  cor.matrix.ri[["r"]] <- M.r2
  return(cor.matrix.ri)
}

process_table <- function(M.r2,gene_order2){
  M.r2$Gene <- rownames(M.r2)
  M.r2 <- pivot_longer(M.r2, cols = -c(Gene),names_to = "ENV",values_to = "r")
  M.r2$ENV_CLASS <- rep("",nrow(M.r2))
  for (i in 1:nrow(M.r2)) {
    if(M.r2[i,]$ENV %in% group1_y){
      M.r2[i,]$ENV_CLASS <- "CLIM"
    }else if(M.r2[i,]$ENV %in% group2_y){
      M.r2[i,]$ENV_CLASS <- "Influent_ENV"
    }else if(M.r2[i,]$ENV %in% group3_y){
      M.r2[i,]$ENV_CLASS <- "Reactor_ENV"
    }
  }
  M.r2$Gene <- factor(M.r2$Gene,levels = rev(gene_order2))
  M.r2$ENV_CLASS <- factor(M.r2$ENV_CLASS,levels = c("CLIM","Influent_ENV","Effluent_ENV"))
  M.r2$ENV <- factor(M.r2$ENV,levels = c(group1_y,group2_y,group3_y))
  return(M.r2)
}

env_dat <- env
env_dat <- decostand(env_dat,MARGIN = 2,method = "standardize")
merged_abundance_C$Pathway <- C_pathway[match(rownames(merged_abundance_C),C_pathway$Gene),]$Pathway
merged_abundance_C.new <- merged_abundance_C%>%dplyr::group_by(Pathway)%>%summarise(across(everything(), sum, na.rm = TRUE))
merged_abundance_C.new <- as.data.frame(merged_abundance_C.new)
rownames(merged_abundance_C.new )<- merged_abundance_C.new$Pathway
merged_abundance_C.new <- merged_abundance_C.new[-1]
merged_abundance_C.new["Overall methanogenic function",] <- colSums(merged_abundance_C.new)
otu <- merged_abundance_C.new
otu <- as.data.frame(t(otu))

cor.matrix.ri <- calculate_cor_matrix(otu,env_dat)
M.r <- cor.matrix.ri[["r"]]
cor.matrix.ri.adj <- p_adjust(cor.matrix.ri)
M.r2 <- cor.matrix.ri.adj[["r"]]

gene_order2 <-  c("Overall methanogenic function","Central methanogenic pathway","Aceticlastic methanogenesis","Hydrogenotrophic methanogenesis","Methylotrophic methanogenesis")
# 用于标注显著性
M.r2 <- process_table(M.r2,gene_order2)
M.p2 <- process_table(cor.matrix.ri.adj[["p"]],gene_order2)
M.p2$significance <- rep("",nrow(M.p2))
M.p2[M.p2$r>0.05,]$significance <- ""
M.p2[M.p2$r<=0.001,]$significance <- "***"
M.p2[M.p2$r<=0.01&M.p2$r>0.001,]$significance <- "**"
M.p2[M.p2$r<=0.05&M.p2$r>0.01,]$significance <- "*"
M.r3 <- M.r
M.r3 <- process_table(M.r3,gene_order2)

p <- ggplot(M.r3, aes(x = ENV, y = Gene)) +
  geom_square(aes(r0=r,fill=r),color=NA)+
  geom_square(data = M.r2,aes(x = ENV, y = Gene,r0=r), fill = "NA",colour = "#457033",size = 5, linetype = "dashed")+
  geom_text(data = M.p2,aes(label = significance), size = 5)+
  facet_grid(cols =  vars(ENV_CLASS), scales = "free",space = "free")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1),
                       space = "Lab")+
  theme_classic()+
  theme(plot.title = element_text(face = "bold", family = "sans"),
        axis.text = element_text(family = "serif",color = 'black'),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "#F5F5F5"),
        strip.text = element_text(face = "bold", family = "serif",size = 12),
        legend.text = element_text(face = "bold.italic", family = "serif"),
        legend.title = element_text(face = "bold.italic", family = "serif",hjust=0.5),
        panel.background = element_rect(colour = "black", linewidth = 0.5, linetype = "solid",fill = NULL),
        panel.grid = element_blank())+
  guides(fill = guide_colorbar(title = "Spearman's r"))

p
