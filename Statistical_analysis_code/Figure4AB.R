##### Phylogenetic features of methanogenesis functional genes in relation to gas production performance #####
rm(list = ls())

setwd("/AD_Methane_Data")

library(ggplot2)
library(cowplot)
library(corrplot)
library(tidyr)
library(dplyr)
library(MetBrewer)
library(linkET)

#####read env+isc data#####
group1_y <-c("T2M", "WS2M", "RH2M", "PRECTOTCORR")
group2_y <- c("I_TN","I_NH3_N","I_pH","I_Salinity","I_TS","I_SCOD","I_S_carbohydrate","I_S_protein")
group3_y <- c("R_TN","R_NH3_N","R_pH","R_Salinity","R_T","R_TS","R_SCOD","R_S_carbohydrate","R_S_protein")
group4_y <- c("MC","GP","MP","AGP","AMP")

#filtered env
env_filtered <-  read.delim('Environmental_Factors.txt', row.names = 1)
#filter na
env_filtered2 <- env_filtered%>%select_if(~!any(is.na(.)))
isc <-  read.delim('Performance_Indicators.txt', row.names = 1)
selected.isc <- c("MC","MP","AMP","GP","AGP")
CLIM <- read.delim('Climatic_Factors.txt', row.names = 1)
selected.clim <- c("T2M", "WS2M", "RH2M", "PRECTOTCORR")
env <- cbind(CLIM[selected.clim],env_filtered2,isc[selected.isc])
env_dat <- env[,!(names(env) %in% c("Longitude","Latitude"))]
data_label2site<- read.delim('Metagenome_Sample_List.txt',header = T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
rownames(env_dat) <- data_label2site[match(rownames(env_dat),data_label2site$Sample),]$`Data lable`

C_annotation_data<- read.delim('Methanogenic_Genes_Annotation.txt',header = T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
C_annotation_data$Gene2 <- rep("",nrow(C_annotation_data))
C_annotation_data$Gene2 <- C_annotation_data$Gene
C_annotation_data[C_annotation_data$Pathway=="Hydrogenotrophic methanogenesis",]$Pathway <- "Hydrogenotrophic\nmethanogenesis"
C_annotation_data[C_annotation_data$Pathway=="Aceticlastic methanogenesis",]$Pathway <- "Aceticlastic\nmethanogenesis"
C_annotation_data[C_annotation_data$Pathway=="Methylotrophic methanogenesis",]$Pathway <- "Methylotrophic\nmethanogenesis"

gene_order <- C_annotation_data$Gene
gene_order2 <- C_annotation_data$Gene2
pathway_order <- unique(C_annotation_data$Pathway) 
 
# 定义 x 轴和 y 轴的分组
group1_x <- C_annotation_data[C_annotation_data$Pathway=="Central methanogenic pathway",]$Gene2
group2_x <- C_annotation_data[C_annotation_data$Pathway=="Hydrogenotrophic methanogenesis",]$Gene2
group3_x <- C_annotation_data[C_annotation_data$Pathway=="Aceticlastic methanogenesis",]$Gene2
group4_x <- C_annotation_data[C_annotation_data$Pathway=="Methylotrophic methanogenesis",]$Gene2

calculate_cor_matrix <- function(community_data_list,env_data,community_feature){
  M <- list()
  M.p <- data.frame(matrix(ncol = ncol(env_dat), nrow = length(community_data_list)))
  rownames(M.p) <- names(community_data_list)
  colnames(M.p) <- colnames(env_dat)
  M.r <- M.p
  for(i in names(community_data_list)){
    tmp <- community_data_list[[i]]
    if(sum(!is.na(tmp[,community_feature]))==0){
      M.p[i,j] <- NA
      M.r[i,j] <- NA
      next
    }
    samp.dat <- tmp$sample
    samp.env <- rownames(env_dat)
    matched <- match(samp.dat,samp.env)
    env_dat <- env_dat[matched,]
    for (j in colnames(env_dat)){
      cor <- cor.test(tmp[,community_feature],env_dat[,j],alternative="two.sided",method="spearman")
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
  M.r2$Pathway <- C_annotation_data[match(rownames(M.r2),C_annotation_data$Gene2),]$Pathway
  M.r2 <- pivot_longer(M.r2, cols = -c(Gene,Pathway),names_to = "ENV",values_to = "r")
  M.r2$ENV_CLASS <- rep("",nrow(M.r2))
  for (i in 1:nrow(M.r2)) {
    if(M.r2[i,]$ENV %in% group1_y){
      M.r2[i,]$ENV_CLASS <- "CLIM"
    }else if(M.r2[i,]$ENV %in% group2_y){
      M.r2[i,]$ENV_CLASS <- "Influent_ENV"
    }else if(M.r2[i,]$ENV %in% group3_y){
      M.r2[i,]$ENV_CLASS <- "Reactor_ENV"
    }else if(M.r2[i,]$ENV %in% group4_y){
      M.r2[i,]$ENV_CLASS <- "Performance"
    }
  }
  M.r2$Gene <- factor(M.r2$Gene,levels = rev(gene_order2))
  M.r2$Pathway <- factor(M.r2$Pathway,levels = pathway_order)
  M.r2$ENV_CLASS <- factor(M.r2$ENV_CLASS,levels = c("CLIM","Influent_ENV","Reactor_ENV","Performance"))
  M.r2$ENV <- factor(M.r2$ENV,levels = c(group1_y,group2_y,group3_y,group4_y))
  return(M.r2)
}

#####PD####
PD_dir <- "Methanogenic_Genes_PD/"
list.pd <- list.files(PD_dir,full.names = TRUE)
files.pd <- list.pd[!file.info(list.pd)$isdir]
file_list.pd <- lapply(files.pd, read.delim,header = T,sep = "\t")
file.list.pd <- list()
for (i in 1:length(file_list.pd)){
  file.name <- basename(files.pd[i])
  gene <- unlist(strsplit(file.name,"_"))[1]
  if(gene %in% gene_order){
    tmp <- file_list.pd[[i]]
    tmp <- data.frame(sample=rownames(tmp),PD=tmp$PD)
    tmp$PD <- as.numeric(tmp$PD)
    gene2 <- C_annotation_data[C_annotation_data$Gene==gene,]$Gene2
    file.list.pd[[gene2]] <- tmp
  }
  
}

cor.matrix.ri <- calculate_cor_matrix(file.list.pd,env_dat,"PD")
M.r <- cor.matrix.ri[["r"]]
if(length(which(rowSums(is.na(M.r))>0))>0) {M.r <- M.r[-which(rowSums(is.na(M.r))>0),]}
cor.matrix.ri[["r"]] <- M.r
M.p <- cor.matrix.ri[["p"]]
if(length(which(rowSums(is.na(M.p))>0))>0) {M.p <- M.p[-which(rowSums(is.na(M.p))>0),]}
cor.matrix.ri[["p"]] <- M.p

cor.matrix.ri.adj <- p_adjust(cor.matrix.ri)
M.r2 <- cor.matrix.ri.adj[["r"]]

gene_order2.new <- rownames(M.r)[match(gene_order2,rownames(M.r))]
gene_order2.new <- gene_order2.new[!is.na(gene_order2.new)]

M.r2 <- process_table(M.r2,gene_order2.new)
M.r3 <- M.r
M.r3 <- process_table(M.r3,gene_order2.new)
M.r3 <- M.r3[!M.r3$Pathway=="Anaerobic oxidation\nof methane (AOM)",]
M.r2 <- M.r2[!M.r2$Pathway=="Anaerobic oxidation\nof methane (AOM)",]

p <- ggplot(M.r3, aes(x = ENV, y = Gene)) +
  geom_square(aes(r0=r,fill=r),color=NA)+
  geom_square(data = M.r2,aes(x = ENV, y = Gene,r0=r), fill = "NA",colour = "#457033",size = 5, linetype = "dashed")+
  facet_grid(rows =  vars(Pathway), cols = vars(ENV_CLASS), scales = "free",space = "free", switch = "y")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1),
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
        legend.text = element_text(face = "bold.italic", family = "serif"),
        legend.title = element_text(face = "bold.italic", family = "serif",hjust=0.5),
        panel.background = element_rect(colour = "black", linewidth = 0.5, linetype = "solid",fill = NULL),
        panel.grid = element_blank())+
  guides(fill = guide_colorbar(title = "Spearman's r"))

p

# stat 
stat <- M.r2%>%dplyr::group_by(ENV_CLASS,ENV)%>%summarise(Positive=sum(r>0),Negative=sum(r<0),NS=sum(r==0))
stat$ENV_CLASS <- factor(stat$ENV_CLASS,levels = c("CLIM","Influent_ENV","Reactor_ENV","Performance"))
stat$ENV <- factor(stat$ENV,levels = c(group1_y,group2_y,group3_y,group4_y))
stat.longer <- stat%>%pivot_longer(c(Positive,Negative,NS),names_to = "Correlation",values_to = "Count")
stat.longer$Correlation <- factor(stat.longer$Correlation,levels = c("NS","Positive","Negative"))
stat.longer$Percent <- stat.longer$Count/114*100

p2 <- ggplot(stat.longer, aes(x = ENV, y = Percent, fill = Correlation)) +
  geom_bar(width = 0.8,stat = "identity",position='stack') +
  facet_grid(cols =  vars(ENV_CLASS), scales = "free",space = "free",switch = "x")+
  scale_x_discrete(position = "bottom") +
  # scale_y_reverse()+
  labs(title = "PD")+
  scale_color_manual(values=c("lightgrey","red","blue"))+
  scale_fill_manual(values=c("lightgrey","red","blue"))+
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


#####MPD####
PD_dir <- "Methanogenic_Genes_PD/"
list.pd <- list.files(PD_dir,full.names = TRUE)
files.pd <- list.pd[!file.info(list.pd)$isdir]
file_list.pd <- lapply(files.pd, read.delim,header = T,sep = "\t")
file.list.pd <- list()
for (i in 1:length(file_list.pd)){
  file.name <- basename(files.pd[i])
  gene <- unlist(strsplit(file.name,"_"))[1]
  if(gene %in% gene_order){
    tmp <- file_list.pd[[i]]
    if(sum(is.infinite(tmp$MNTD))>0){
      tmp[apply(tmp,1, function(x) any(is.infinite(x))),]$MNTD <- NA
    }
    tmp$sample <- rownames(tmp)
    rownames(tmp) <- c(1:nrow(tmp))
    gene2 <- C_annotation_data[C_annotation_data$Gene==gene,]$Gene2
    file.list.pd[[gene2]] <- tmp
  }
}

cor.matrix.ri <- calculate_cor_matrix(file.list.pd,env_dat,"MPD")

M.r <- cor.matrix.ri[["r"]]
M.r <- M.r[-which(rowSums(is.na(M.r))>0),]
cor.matrix.ri[["r"]] <- M.r
M.p <- cor.matrix.ri[["p"]]
M.p <- M.p[-which(rowSums(is.na(M.p))>0),]
cor.matrix.ri[["p"]] <- M.p

cor.matrix.ri.adj <- p_adjust(cor.matrix.ri)
M.r2 <- cor.matrix.ri.adj[["r"]]

gene_order2.new <- rownames(M.r)[match(gene_order2,rownames(M.r))]
gene_order2.new <- gene_order2.new[!is.na(gene_order2.new)]

# 用于标注显著性
M.r2 <- process_table(M.r2,gene_order2.new)
# 画出所有r
M.r3 <- M.r
M.r3 <- process_table(M.r3,gene_order2.new)
M.r3 <- M.r3[!M.r3$Pathway=="Anaerobic oxidation\nof methane (AOM)",]
M.r2 <- M.r2[!M.r2$Pathway=="Anaerobic oxidation\nof methane (AOM)",]

p <- ggplot(M.r3, aes(x = ENV, y = Gene)) +
  geom_square(aes(r0=r,fill=r),color=NA)+
  geom_square(data = M.r2,aes(x = ENV, y = Gene,r0=r), fill = "NA",colour = "#457033",size = 5, linetype = "dashed")+
  facet_grid(rows =  vars(Pathway), cols = vars(ENV_CLASS), scales = "free",space = "free", switch = "y")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1),
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
        legend.text = element_text(face = "bold.italic", family = "serif"),
        legend.title = element_text(face = "bold.italic", family = "serif",hjust=0.5),
        panel.background = element_rect(colour = "black", linewidth = 0.5, linetype = "solid",fill = NULL),
        panel.grid = element_blank())+
  guides(fill = guide_colorbar(title = "Spearman's r"))

p

# stat 
stat <- M.r2%>%dplyr::group_by(ENV_CLASS,ENV)%>%summarise(Positive=sum(r>0),Negative=sum(r<0),NS=sum(r==0))
stat$ENV_CLASS <- factor(stat$ENV_CLASS,levels = c("CLIM","Influent_ENV","Reactor_ENV","Performance"))
stat$ENV <- factor(stat$ENV,levels = c(group1_y,group2_y,group3_y,group4_y))
stat.longer <- stat%>%pivot_longer(c(Positive,Negative,NS),names_to = "Correlation",values_to = "Count")
stat.longer$Correlation <- factor(stat.longer$Correlation,levels = c("NS","Positive","Negative"))
# rowSums(stat[3:5])
stat.longer$Percent <- stat.longer$Count/112*100

p2 <- ggplot(stat.longer, aes(x = ENV, y = Percent, fill = Correlation)) +
  geom_bar(width = 0.8,stat = "identity",position='stack') +
  facet_grid(cols =  vars(ENV_CLASS), scales = "free",space = "free",switch = "x")+
  scale_x_discrete(position = "bottom") +
  labs(title = "MPD")+
  scale_color_manual(values=c("lightgrey","red","blue"))+
  scale_fill_manual(values=c("lightgrey","red","blue"))+
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

