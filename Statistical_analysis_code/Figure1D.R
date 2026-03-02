######Partial Mantel tests (geographic distance-corrected)#####
set.seed(123)
setwd("/AD_Methane_Data")
rm(list=ls())
library(dplyr)
taxon_16S <-  read.delim('AD_OTU_Table.txt', row.names = 1)
taxon_16S <- as.data.frame(t(taxon_16S)) 

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
#####match otu and env#####
samp.otu.16S = rownames(taxon_16S)
samp.env = rownames(env)
my.otu.16S = match(samp.otu.16S,samp.env)
env.16S = env[my.otu.16S,]

# install.packages("devtools")
# devtools::install_github("Hy4m/linkET", force = TRUE)

library(linkET)
library(ggplot2)
library(dplyr)
library(ggcorrplot)
library(vegan)
library(viridis)
library(ecodist)


mantel.prepare.data <- function(env_data,taxon,taxon_name,partial,sorensen=F){
  taxon.dis = vegdist(taxon, method="bray",binary=sorensen)
  geo <- env_data[c('Longitude', 'Latitude')]
  env_data <- decostand(env_data,method = "standardize",MARGIN = 2)
  env_data2 <- env_data[,-1]
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

#partial.mantel
m.16s.p <- mantel.prepare.data(env.16S,taxon_16S,"Bray-Curtis",partial="TRUE")
mantel <- rbind(m.16s.p)
mantel.tmp <- mantel


# u-unifrac
data_file_path <- "UniFrac_Matrix/Community_Unweighted_Unifrac_matrix.txt"
unifrac <- read.delim(data_file_path, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE,row.names = 1)

#####match otu and env#####(必做)
samp.otu.16S <- rownames(unifrac)
samp.env = rownames(env)
my.otu.16S = match(samp.otu.16S,samp.env)
env.16S = env[my.otu.16S,]

u.dist <- as.dist(unifrac)

mantel.prepare.data <- function(env_data,taxon,taxon_name,partial,distance_matrix,sorensen=F){
  if(is.null(taxon)){
    taxon.dis = distance_matrix
  }else{
    taxon.dis = vegdist(taxon, method="bray",binary=sorensen)
  }
  #计算单个环境因子距离
  geo <- env_data[c('Longitude', 'Latitude')]
  env_data <- decostand(env_data,method = "standardize",MARGIN = 2)
  env_data2 <- env_data[,-1]
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

#partial.mantel
m.16s.p.uu <- mantel.prepare.data(env.16S,NULL,"Unweighted UniFrac","TRUE",u.dist)
mantel <- rbind(mantel.tmp,m.16s.p.uu)

for(i in 3:ncol(mantel)){mantel[,i] <- as.numeric(mantel[,i])}
mantel <- mutate(mantel,
                 rd = cut(r, breaks = c(-Inf, 0.1, 0.2, Inf),
                          labels = c("< 0.1", "0.1 - 0.2", ">= 0.2")),
                 pd = cut(p, breaks = c(-Inf, 0.05, Inf),
                          labels = c("< 0.05",">= 0.05")))
mantel

env2 <- na.omit(env[-which(colnames(env)=="Longitude")])
correlate(env2,method = "pearson",engine = "Hmisc")
#筛选p值小于等于0.05
cor_p <- Hmisc::rcorr(as.matrix(env2),type = "pearson")
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
  mutate(value.x = as.character(value.x),
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

p
