# 1. library -----
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(DOSE)  
library(presto)
if(!require(dplyr))install.packages("dplyr")
library(dplyr)

library(multtest)
if(!require(Seurat))install.packages("Seurat")
library(DoubletFinder)

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
library(ComplexHeatmap)
# 2. load data -------
load("immune/immune.annotated.Rdata")
sce<-imu[, imu$new_celltype %in% c( 'Neu_Csf1', 'Neu_Retnlg' )]

#3. cluster profile --------
type <- unique(sce.all.1$new_celltype)
for (i in 1:length(type)) {
  degs_data <- markers %>% filter(cluster==type[[i]]) %>% arrange(desc(avg_log2FC))
  degs_data <- markers %>% dplyr::filter(cluster ==c("Neu_Csf1","Neu_Retnlg")) %>% arrange(desc(avg_log2FC)) 
  #%>% head(100)
  genelist <- degs_data$avg_log2FC
  names(genelist) <- degs_data$gene
  genelist <- sort(genelist, decreasing = TRUE) 
  genelist_entrez <- mapIds(org.Rn.eg.db,            
                            keys = names(genelist),                  
                            column = "ENTREZID",         
                            keytype = "SYMBOL",                  
                            multiVals = "first")  

  entrez_list <- genelist[!is.na(genelist_entrez)]
  names(entrez_list) <- genelist_entrez[!is.na(genelist_entrez)]
  GO_ges <- gseGO(entrez_list,          
                  OrgDb = org.Rn.eg.db,        
                  ont = "ALL",               
                  minGSSize = 10,             
                  maxGSSize = 500,        
                  pvalueCutoff = 0.05) 
  
  GO_ges <- setReadable(GO_ges, OrgDb = org.Rn.eg.db,
                        keyType = "ENTREZID")
  GO_ges_result <- GO_ges@result
  filename = paste0(unique(type[i]),".markergene.gsea.go.csv")
  write.csv(GO_ges_result,filename)
}

kegg <- read.csv("allcellgsea/NEU.go.csv")
kegg_csf1 <- kegg %>% filter(Celltype == "Csf1")
level_csf1 <- kegg_csf1$Description[kegg_csf1$Celltype=='Csf1']
kegg_csf1$Description <- factor(kegg_csf1$Description,
                                levels = rev(level_csf1))
p1 <- ggplot(kegg_csf1,aes(x = NES,y = Description,fill=Celltype,color=pvalue)) + 
  geom_col(alpha=0.8) +
  scale_fill_manual(values =  "#CC4C02") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(color="black",size=10),
        axis.ticks.x=element_line(size = 1  ),
        axis.line.x=element_line(size = 1),
        # y轴设置
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(), 
        axis.text.y=element_text(size = 10),
        legend.title=element_blank(),
        legend.spacing.x=unit(0.2,'cm'),
        legend.key=element_blank(),
        legend.key.width=unit(0.5,'cm'),
        legend.key.height=unit(0.5,'cm'),
        plot.margin = margin(1,0.5,0.5,1,unit="cm"))
p1
kegg_Retnlg <- kegg %>% filter(Celltype == "Retnlg")
level_Retnlg <- kegg_Retnlg$Description[kegg_Retnlg$Celltype=='Retnlg']
kegg_Retnlg$Description <- factor(kegg_Retnlg$Description,
                                  levels = rev(level_Retnlg))
p2 <- ggplot(kegg_Retnlg,aes(x = NES,y = Description,fill=Celltype)) + 
  geom_col(alpha=0.8) +
  scale_fill_manual(values = "#7A5399") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(color="black",size=10),
        axis.ticks.x=element_line(size = 1  ),
        axis.line.x=element_line(size = 1),
        # y轴设置
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(), 
        axis.text.y=element_text(size = 10),
        legend.title=element_blank(),
        legend.spacing.x=unit(0.2,'cm'),
        legend.key=element_blank(),
        legend.key.width=unit(0.5,'cm'),
        legend.key.height=unit(0.5,'cm'),
        plot.margin = margin(1,0.5,0.5,1,unit="cm"))
p2
p1/p2


#4. 5 genes vlnpolt ------
genes <- c("Mmp8","Mmp9","Ccl3",'Itgam','Itgb2')
p1 <- VlnPlot(sce,  
              features = genes,
              ncol=5,
              # stack=T,  
              pt.size=0.5,  #flip=T,
              split.by = 'new_celltype',
              cols=c("#9970AB","#E08214") 
) +
  theme(legend.position = "none") 
p1

#5. NETs formation related genes heatmap ------
symbol <- cluster414.genes$feature
entrez <- bitr(symbol, 
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Rn.eg.db")

sce.genes <- wilcoxauc(sce,'seurat_clusters')
dplyr::count(sce.genes,group) 
cluster414.genes <- sce.genes %>% dplyr::filter(group ==c("4","14")) %>% 
  arrange(desc(auc)) %>% dplyr::select(feature,auc)

genelist <- cluster414.genes$auc
names(genelist) <- cluster414.genes$feature
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
head(genelist)

KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "rno",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

KEGG_ges <- setReadable(KEGG_ges,
                        OrgDb = org.Rn.eg.db,
                        keyType = "ENTREZID")

library(GseaVis)

#Neutrophil extracellular trap formation
gseaNb(object = KEGG_ges,
       geneSetID = 'Neutrophil extracellular trap formation')
gseaNb(object = KEGG_ges,
       geneSetID = 'rno04613',
       addPval = T,
       pvalX = 0.75,pvalY = 0.75,
       pCol = 'black',
       subPlot = 3,rmHt =T,
       pHjust = 0)

#heatmap
NETs.genes <- c(
  "Clec7a",'Selplg','Fcgr1a','Itgam','Itgb2','H3f3a',
  'Rac2','H2az1','Actg1','H3f3b','Actb',
  'Cyba','Ncf2','Ncf4',
  'Padi4','Elane',"Tlr4","Tlr2","Syk",  'Il1b',
  'Osm','Cxcl2','Cxcl1','Nos2','Ptgs2','Tnf',
  'Nfkb1','Nfkb2','Nfkbia'           
)
NETs.gene.exp <- AverageExpression(sce,features = NETs.genes,
                                   group.by = 'new_celltype',
                                   slot = 'data')

NETs.gene.exp <- as.data.frame(NETs.gene.exp$RNA)

NETs_exp <- t(scale(t(NETs.gene.exp),scale=T,center = T))

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)
Breaks = c(seq(min(NETs_exp), 0, 
               length.out=ceiling(paletteLength/2) + 1),
           seq(max(NETs_exp)/paletteLength, 
               max(NETs_exp), 
               length.out=floor(paletteLength/2)))

pheatmap(NETs_exp ,
         fontsize=14, 
         cluster_cols  = F,
         cluster_rows = F,
         fontsize_row = 10, 
         color=myColor, breaks = Breaks, 
         #angle_col = 90,
         treeheight_col = 0,  border_color = "white")


