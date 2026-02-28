# 1. library -------
library(multtest)
if(!require(Seurat))install.packages("Seurat")
library(DoubletFinder)
if(!require(dplyr))install.packages("dplyr")
#if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(clustree))install.packages("clustree")
if(!require(hdf5r))install.packages("hdf5r")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
library(ggplot2)
library(Matrix)
library(harmony)
library(future)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(scRNAtoolVis)
library(GSVA) # BiocManager::install('GSVA')
library(GSEABase)
library(escape) #BiocManager::install("escape")


# 2. load data-------
load('all.cell.annotation.counts.new.Rdata')

# 3. gene trans -------
a <- read.csv("human.inflammation.genes.symbol.csv",header = T)
library(org.Rn.eg.db)
library("biomaRt")
library(dplyr)
hgene <- a$ENSEMBL
human <- read.table("mart_export.txt",header = TRUE,sep="\t")
test <- human[human$Gene.stable.ID %in% hgene,]
test <- distinct(test[,c(1,5,6,8)])

# 4. inflammaton feature plot -------
Rat.infl.gene <- read.table("inflammation/Rat.inflammation.genes.symbol.txt",header=F)
Rat.infl.gene <- Rat.infl.gene$V1
sce.all.1 <- AddModuleScore(sce.all.1,features = list(Rat.infl.gene),name = 'Inflammation.scores')

library(scCustomize)
library(Nebulosa)

#BiocManager::install("Nebulosa")

FeaturePlot_scCustom(sce.all.1, features= 'Inflammation.scores1',order = TRUE,split.by = 'group',
                     pt.size = 0.1) 
ggsave("inflammation/FeaturePlot.inflammation.score.pdf")

# 5. inflammation score vlnplot -------
t <- sce.all.1@meta.data[c(1:3,4,19,20)]
t$group <- factor(t$group,levels=c("Control","CD","PB_IFX"))

p1 <- ggplot(t,aes(group,Inflammation.scores1,colour=group))+
  geom_violin(trim=T)+
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE)+
  theme_classic()+ 
  scale_color_manual(values=c("#9B989B","#BE0F20","#272B82"))+
  geom_signif(             
    comparisons = list(c("CD", "Control"),                               
                       c("CD", "PB_IFX"),                               
                       c("Control", "PB_IFX")),             
    map_signif_level=T,             
    tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0),         
    y_position = c(0.3,0.35,0.4),          
    size=1,      
    textsize = 4,         
    test = "wilcox.test")+
  ylab("Inflammation scores")+xlab("")+
  theme(axis.text.x = element_text(colour = 'black',size = 12,angle = 40,hjust = 1, vjust = 1),   
        axis.text.y = element_text(colour = 'black',size = 12),     
        axis.title = element_text(colour = 'black',size = 15),     
        axis.line = element_line(color = 'black', size = 1))

p2 <-ggplot(t,aes(group,Inflammation.scores1,colour=group))+
  geom_violin(trim=T)+
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE)+
  theme_classic()+ 
  scale_color_manual(values=c("#9B989B","#BE0F20","#272B82"))+
  stat_compare_means(method="wilcox.test",hide.ns = F,    
                     comparisons =list(c("CD", "Control"),                              
                                       c("CD", "PB_IFX")                           
                                                           ),     
                     label = "p.value",
                     bracket.size=0.6,size=4)+
  ylab("Inflammation scores")+xlab("")+
  theme(axis.text.x = element_text(colour = 'black',size = 12,angle = 40,hjust = 1, vjust = 1),   
        axis.text.y = element_text(colour = 'black',size = 12),     
        axis.title = element_text(colour = 'black',size = 15),     
        axis.line = element_line(color = 'black', size = 1))
p1+p2






