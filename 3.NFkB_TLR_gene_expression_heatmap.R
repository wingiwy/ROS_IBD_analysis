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
library(progeny) #devtools::install_github("saezlab/progeny")
library(tidyr)
library(readr)
library(pheatmap)
#install.packages("devtools")
#devtools::install_github("junjunlab/jjAnno")

# 2. load data -------
load('all.cell.annotation.counts.new.Rdata')

#### 2.1 CD group -------
#提取CD组
sce = sce.all.1[, sce.all.1$group %in% c( 'CD' )]

#### 2.2 pathway score analysis -------
CellsClusters <- data.frame(Cell = names(Idents(sce)), 
                            CellType = as.character(Idents(sce)),
                            stringsAsFactors = FALSE)

sce <- progeny(sce,scale = FALSE,organism = "Mouse",
               top=500,perm=1,return_assay = T)

sce <- Seurat::ScaleData(sce, assay = "progeny") 

progeny_scores_df <- 
  as.data.frame(t(GetAssayData(sce, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

dim(progeny_scores_df)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))


summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 


paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))

df <- t(summarized_progeny_scores_df[,-1])
pheatmap(df,
         fontsize=14, 
         cluster_cols  = F,
         cluster_rows = F,
         fontsize_row = 10, 
         color=myColor, breaks = progenyBreaks, 
         angle_col = 90,
         treeheight_col = 0,  border_color = "grey")

#### 2.3 cytokine gene expression dotplot -------
genes <- c(
  'Nos2','Ptgs2',
  'Ccl2','Cxcl2','Cxcl1','Icam1','Vcam1',
  'Osm','Il18','Il6','Il1b','Il1a',
  'Tnf',  'Acod1','Rela','Myd88',
  'Tlr8',"Tlr7","Tlr6","Tlr5","Tlr4",'Tlr3','Tlr2',"Tlr1", 
  'Il1r2','Il1r1',  "Tnfrsf1b","Tnfrsf1a"
)
p<-DotPlot(sce_CD,features = genes, group.by ="new_celltype")
p$data$id <- factor(p$data$id,levels = c(
  'Stem','EECs_Neurod1','EECs_Lix1','Undiff_Entero','Entero_Rbp2','Entero_Mep1a',
  'Entero_Sult1a1','Entero_Duoxa2_s100a9','Goblet_Best2','Goblet_Cenpa',
  'Goblet_Gmds','Goblet_pre','Mucous_cells','Tuft','TA',
  'Mast','pDCs','cDC2_Zeb2','cDC1',
  'Macs_Mrc1_Cd163','Macs_Itgax','Macs_Spp1','Ly6c_Monocytes',
  'Neu_Csf1','Neu_Retnlg',
  'Cd8_Gzmm', 'Cd4_Cd8','Cd4_treg',
  'IgG_plasma_B','Memory_B','Plasmablast',
  'Endothelial','Lymphatics','Fibro_IL11_Cxcl14','Fibro_Smoc2_Ccl11',
  'Fibro_Mfap5_Gsn','Glial','Myofibroblasts','Cycling_stroma'
))

pCD <-  ggplot(p$data,aes(x = features.plot,y = id)) +  
  geom_point(aes(fill = avg.exp.scaled,size = pct.exp),
             color='black',shape=21) +  
  theme_bw(base_size = 14) +  xlab('') + ylab('') +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0,                     
                       limits = c(-2,3),                    
                       breaks = seq(-2,3,1),        
                       labels = seq(-2,3,1),   
                       name = 'Mean expression') + 
  scale_size(range = c(0,6),            
             limits = c(0, 100),          
             breaks = seq(20,100,20),       
             labels = seq(20,100,20)  ) + 
  theme( legend.position = 'top' , 
         panel.grid = element_blank(),    
         axis.text = element_text(color = 'black'),    
         axis.text.y = element_text( face = 'italic'),
         axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5 ) 
  ) + 
  coord_flip()+
  guides(  
    fill = guide_colorbar(    
      direction = "horizontal",   
      title.position = "top",    
      barwidth = unit(5, "cm")  
    ), 
    size = guide_legend(     
      title = "Percent of cells",      
      direction = "horizontal",      
      title.position = "top",     
      label.position = "right",     
      override.aes = list(      
        color = "black",      
        fill = "grey"     
      ) ) )
pCD 

#### 2.4 ROS-related gene expression vlnplot -------

sce_mye <- sce[, sce$new_celltype %in% c( 'Mast',
                                          'pDCs','cDC2_Zeb2','cDC1',
                                          'Macs_Mrc1_Cd163','Macs_Itgax','Macs_Spp1','Ly6c_Monocytes',
                                          'Neu_Csf1','Neu_Retnlg'
)]
ROS.genes <- c(
  'Ncf1','Ncf2','Ncf4',
  'Cyba','Cybb',  #Cybb=Nox2
  'Nox1', 'Nox4' ,   #NADPH酶家族，促氧化酶
  'Nfe2l2',"Hmox1",'Nqo1','Nlrp3','Casp1'
)

Idents(sce_mye) <- factor(Idents(sce_mye), levels= c('Mast','pDCs','cDC2_Zeb2','cDC1',
                                                     'Macs_Mrc1_Cd163','Macs_Itgax','Macs_Spp1','Ly6c_Monocytes',
                                                     'Neu_Csf1','Neu_Retnlg',"Fibro_IL11_Cxcl14","Cycling_stroma"))
p1 <- VlnPlot(sce_mye,  
              features = ROS.genes,  
              stack=T,  pt.size=0,  #flip=T,
              split.by = 'new_celltype',    
              cols=c("#006ADD","#C10000","#00A500","#6B5BC2","#FF8800","#13C4FF","#89194A","#00770F","#135988"),  ) +
  theme(legend.position = "none") 
p1

#2. ROS generation score ------
ros_geneset <- list(  ROS_Generation = c("Cybb","Nox4","Duox2","Nox1","Cyba", 
                                         "Ncf1","Ncf2","Xdh","Ptgs2","Nos2") )

sce <- AddModuleScore(sce, features = ros_geneset, name = "ROS_", assay = "RNA")
sce@meta.data$ROS_Score <- sce@meta.data$ROS_1 

#3. NLRP3 score ------
nlrp3_genes <-list( nlrp3= c("Nlrp3","Casp1","Il1b","Casp4","Gsdmd","Il18","Txnip","Slc7a10",
                             "Gsdme", "Pycard", "Aim2", "Nlrc4", "Nek7"))

sce <- AddModuleScore(sce, features = nlrp3_genes, name = "NLRP3_", assay = "RNA")
sce@meta.data$NLRP3_Score <- sce@meta.data$NLRP3_1 

#4. Response to ROS score ---------
test <- list(  ROS_genes = c("Abcc1",	"Atox1",	"Cat",	"Cdkn2d",	"Egln2",	"Ercc2",
                             "Fes",	"Ftl",	"G6pd",	"Gclc",	"Gclm",	"Glrx",	"Glrx2",'Hmox1',
                             "Gpx3",	"Gpx4",	"Gsr",	"Hbxip",	"Hhex",	"Hmox2",	"Ipcef1",
                             "Junb",	"Lsp1",	"Mbp",	"Mgst1",	"Mpo",	"Msra",	"Ndufa6",	
                             "Ndufb4",	"Ndufs2",	"Nqo1",	"Oxsr1",	"Pdlim1",	"Pfkp",	"Ppp2r4",	
                             "Prdx1",	"Prdx2",	"Prdx4",	"Prdx6",	"Prnp",	"Sbno2",	"Scaf4",
                             "Sels",	"Sod1",	"Sod2",	"Srxn1",	"Stk25",	"Txn",	"Txnrd1",	"Txnrd2"))

sce <- AddModuleScore(sce, features = test, name = "Response_", assay = "RNA")
sce@meta.data$Response_Score <- sce@meta.data$Response_1

#5. vlnplot ---------
plot_data <- sce@meta.data %>%
  dplyr::select(orig.ident, new_celltype,PathwayNFkB,PathwayTNFa,ROS_Score,
                NLRP3_Score,Response_Score
  )  

plot_data$new_celltype <- factor(plot_data$new_celltype,
                                 levels = c(
                                   'Stem','EECs_Neurod1','EECs_Lix1','Undiff_Entero','Entero_Rbp2','Entero_Mep1a',
                                   'Entero_Sult1a1','Entero_Duoxa2_s100a9','Goblet_Best2','Goblet_Cenpa',
                                   'Goblet_Gmds','Goblet_pre','Mucous_cells','Tuft','TA',
                                   'Mast','pDCs','cDC2_Zeb2','cDC1',
                                   'Macs_Mrc1_Cd163','Macs_Itgax','Macs_Spp1','Ly6c_Monocytes',
                                   'Neu_Csf1','Neu_Retnlg',
                                   'Cd8_Gzmm','Cd4_Cd8','Cd4_treg',
                                   'IgG_plasma_B','Memory_B','Plasmablast',
                                   'Endothelial','Lymphatics','Fibro_IL11_Cxcl14','Fibro_Smoc2_Ccl11',
                                   'Fibro_Mfap5_Gsn','Glial','Myofibroblasts','Cycling_stroma'
                                 ))

###5.1 ROS generation score -----
p1<- ggplot(plot_data, aes(x = new_celltype, y = ROS_Score, fill = new_celltype)) +
  geom_violin(alpha = 0.6, scale = "width", width = 0.8, trim = TRUE) +
  geom_boxplot(aes(),width = 0.3, alpha = 0.8, outlier.shape = NA) +
  labs(x="",y="ROS generation score")+
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12)
  ) 

###5.2 NF-kB Score vlnplot -------
p3<- ggplot(plot_data, aes(x = new_celltype, y = PathwayNFkB, fill = new_celltype)) +
  geom_violin(alpha = 0.6, scale = "width", width = 0.8, trim = TRUE) +
  geom_boxplot(aes(),width = 0.3, alpha = 0.8, outlier.shape = NA) +
  labs(x="",y="NF-κB activity score")+
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 10),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12)
  ) 
ggplot(plot_data, aes(x = new_celltype, y = PathwayTNFa, fill = new_celltype)) +
  geom_violin(alpha = 0.6, scale = "width", width = 0.8, trim = TRUE) +
  geom_boxplot(aes(color = new_celltype),width = 0.3, alpha = 0.8, outlier.shape = NA) +
  labs(x="",y="PathwayTNFa activity score")+
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 10),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12)
  ) 

###5.3 NLRP3 Score vlnplot -------
p2<- ggplot(plot_data, aes(x = new_celltype, y = NLRP3_Score, fill = new_celltype)) +
  geom_violin(alpha = 0.6, scale = "width", width = 0.8, trim = TRUE) +
  geom_boxplot(aes(),width = 0.3, alpha = 0.8, outlier.shape = NA) +
  labs(x="",y="Pyroptosis score")+
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 10),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12)
  ) 

p4<- ggplot(plot_data, aes(x = new_celltype, y = Response_Score, fill = new_celltype)) +
  geom_violin(alpha = 0.6, scale = "width", width = 0.8, trim = TRUE) +
  geom_boxplot(aes(),width = 0.3, alpha = 0.8, outlier.shape = NA) +
  labs(x="",y="Response to ROS score")+
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 10),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12)
  ) 
p1/p2/p4
#6. NF-kB. cor ROS score ------
expression_data <-  subset(plot_data, new_celltype == "Neu_Retnlg")
p1 <- ggplot(expression_data, aes(x = ROS_Score, y = PathwayNFkB)) +
  geom_point(alpha = 0.5,size=0.8) +
  geom_smooth(method = "lm",formula = y ~ x, 
              color = "blue", fill = "grey70", alpha = 0.8, 
              se = TRUE, level = 0.95) +
  stat_cor(method = "pearson", label.sep = "\n") +
  theme_classic() +theme( panel.border = element_rect(color = "black", fill = NA, size = 0.8))+
  labs(x = "ROS score", y = "NF-kB activity score", 
       title = "Retnlg Neu")
expression_data <-  subset(plot_data, new_celltype == "Neu_Csf1")

p2<- ggplot(expression_data, aes(x = ROS_Score, y = PathwayNFkB)) +
  geom_point(alpha = 0.5,size=0.8) +
  geom_smooth(method = "lm",formula = y ~ x, 
              color = "blue", fill = "grey70", alpha = 0.8, 
              se = TRUE, level = 0.95) +
  stat_cor(method = "spearman", label.sep = "\n") +
  theme_classic() +theme( panel.border = element_rect(color = "black", fill = NA, size = 0.8))+
  labs(x = "ROS score", y = "NF-kB activity score", 
       title = "Csf1 Neu")

expression_data <-  subset(plot_data, new_celltype == "Ly6c_Monocytes")

p3<- ggplot(expression_data, aes(x = ROS_Score, y = PathwayNFkB)) +
  geom_point(alpha = 0.5,size=0.8) +
  geom_smooth(method = "lm",formula = y ~ x, 
              color = "blue", fill = "grey70", alpha = 0.8, 
              se = TRUE, level = 0.95) +
  stat_cor(method = "spearman", label.sep = "\n") +
  theme_classic() +theme( panel.border = element_rect(color = "black", fill = NA, size = 0.8))+
  labs(x = "ROS score", y = "NF-kB activity score", 
       title = "Ly6c Monocytes")

expression_data <-  subset(plot_data, new_celltype == "Macs_Spp1")

p4 <- ggplot(expression_data, aes(x = ROS_Score, y = PathwayNFkB)) +
  geom_point(alpha = 0.5,size=0.8) +
  geom_smooth(method = "lm",formula = y ~ x, 
              color = "blue", fill = "grey70", alpha = 0.8, 
              se = TRUE, level = 0.95) +
  stat_cor(method = "spearman", label.sep = "\n") +
  theme_classic() +theme( panel.border = element_rect(color = "black", fill = NA, size = 0.8))+
  labs(x = "ROS score", y = "NF-kB activity score", 
       title = "Macs_Spp1")

expression_data <-  subset(plot_data, new_celltype == "Macs_Itgax")

p5 <- ggplot(expression_data, aes(x = ROS_Score, y = PathwayNFkB)) +
  geom_point(alpha = 0.5,size=0.8) +
  geom_smooth(method = "lm",formula = y ~ x, 
              color = "blue", fill = "grey70", alpha = 0.8, 
              se = TRUE, level = 0.95) +
  stat_cor(method = "spearman", label.sep = "\n") +
  theme_classic() +theme( panel.border = element_rect(color = "black", fill = NA, size = 0.8))+
  labs(x = "ROS score", y = "NF-kB activity score", 
       title = "Macs_Itgax")

(p1+p2+p3+p4+p5)


