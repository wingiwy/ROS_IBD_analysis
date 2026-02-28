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

# 2. working dir -------
#### 2.1 dir -------
rm(list = ls()); gc () 
dirs = dir(pattern="CD_|Control|PB")  
f = "dat.Rdata"
if(!file.exists(f)){scelist = list()}

tmp = list.dirs('/work/project/Rat_CD_model/single_cell_sequencing')[-1] 

ct = Read10X(tmp)  
sce.all=CreateSeuratObject(counts = ct  ,   # 创建Seurat对象
                           min.cells = 3,
                           min.features = 200)

new.cluster.ids <-  basename(tmp)
names(new.cluster.ids) <- levels(sce.all$orig.ident)
pbmc <- RenameIdents(sce.all, new.cluster.ids)
sce.all$orig.ident = Idents(pbmc)  

aa <- as.character(sce.all@meta.data$orig.ident) 
names(aa) <- rownames(sce.all@meta.data) 
sce.all@active.ident <- as.factor(aa) 

sce.all$group = sce.all$orig.ident

# 3. QC -------
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all,pattern = "^Mt-")

sce.all[["percent.rb"]] <- PercentageFeatureSet(sce.all,pattern = "^Rbs|Rpl")

sce.all[["percent.hb"]] <- PercentageFeatureSet(sce.all,pattern = "^Hb[^(p)]")

sce.all[["lograti"]] <- log10(sce.all$nFeature_RNA) / log10(sce.all$nCount_RNA)

sce.all <- subset(sce.all, subset = nFeature_RNA < 6000 &  percent.mt < 20 & lograti > 0.8 )

counts <- GetAssayData(object = sce.all, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
sce.all <- CreateSeuratObject(filtered_counts, meta.data = sce.all@meta.data)

sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize", scale.factor = 10000)
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(sce.all), 20)  

sce.all <- ScaleData(sce.all)
sce.all <- RunPCA(sce.all, npcs=30, features = VariableFeatures(object = sce.all), 
                  ndims.print=1:20, nfeatures.print=5, seed.use=114514)
VizDimLoadings(sce.all, dims = 1:2, reduction = "pca") 
DimPlot(sce.all, reduction = "pca") + NoLegend()
DimHeatmap(sce.all, dims = 1:15, cells = 500, balanced = TRUE) 

sce.all <- JackStraw(sce.all, num.replicate = 100)  
sce.all <- ScoreJackStraw(sce.all, dims = 1:20)  
JackStrawPlot(sce.all, dims = 1:20)  
ElbowPlot(sce.all)  
ggsave("sce.all.elbowplot.pdf",height=6,width=6)

sce.all <- RunHarmony(sce.all,group.by.vars = "orig.ident")
sce.all <- RunUMAP(sce.all, reduction = "harmony", dims = 1:20)
DimPlot(sce.all,reduction = "umap",label=F ) 
sce.all <- RunTSNE(sce.all, dims = 1:20, reduction = "harmony")
DimPlot(sce.all,reduction = "tsne",label=F )

sce.all <- FindNeighbors(sce.all, dims = 1:20)
sce.all <- FindClusters(sce.all, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1))
clustree(sce.all@meta.data, prefix = "RNA_snn_res.")

colnames(sce.all@meta.data)
apply(sce.all@meta.data[,grep("RNA_snn",colnames(sce.all@meta.data))],2,table) 

sce.all.markers <- FindAllMarkers(sce.all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 4. annotation ----------

sce.all.0.01 <- SetIdent(sce.all,value="RNA_snn_res.0.01")
table(sce.all.0.01@active.ident)
sce.all.0.01.cluster.markers <-  FindAllMarkers(sce.all.0.01, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10pos <- sce.all.0.01.cluster.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC) 
top10negtive <- sce.all.0.01.cluster.markers %>% group_by(cluster) %>% top_n(n = -10, wt = avg_log2FC)
top_bottom_gene <-data.frame(top10pos,top10negtive)

genes_to_check_3_main_cluster = c(
  "Epcam","Krt8","Krt18",        #Epithelial
  "Cdh5","Col1a1","Col1a2","Col6a2","Vwf",   #Stromal
  'RT1-Db1',"Ptprc","Cd3d","Cd3g","Cd3e",'Cd8a','Trdc', 'ENSRNOG00000071219', 
  "Cd79a","Cd79b","Cd14","Cd68","Cd83","Csf1r","Fcer1g"   #Immune
)
p = DotPlot(sce.all.0.01, features = unique(genes_to_check_3_main_cluster),
            assay='RNA'  )  + coord_flip()

tsne = DimPlot(sce.all.0.01, reduction = "tsne",
               pt.size = 0.8,
               group.by = "RNA_snn_res.0.01",label = T,label.box = T)

umap = DimPlot(sce.all.0.01, reduction='umap',
               pt.size = 0.8,
               group.by = "RNA_snn_res.0.01",label = T,label.box = T)
##### 
celltype=data.frame(ClusterID=0:9,
                    celltype= 0:9) 

celltype[celltype$ClusterID %in% c( 2,6),2]='Stromal'
celltype[celltype$ClusterID %in% c( 0,1,5,9 ),2]='Epithelial'
celltype[celltype$ClusterID %in% c(3,4,7,8),2]='Immune'

sce.all.0.01@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all.0.01@meta.data[which(sce.all.0.01@meta.data$RNA_snn_res.0.01 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all.0.01@meta.data$celltype)


th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 

celltype_tsne = DimPlot(sce.all.0.01, reduction = "tsne",
                        cols=c("#D68747","#BEABD2","#754D40"),
                        pt.size = 1,
                        group.by = "celltype",label = T,label.size=7)

ggsave(plot=celltype_tsne,filename="3_main_cluster/tsne.3mainclus.name.pdf",height=8,width=9)

celltype_umap = DimPlot(sce.all.0.01, reduction = "umap",
                        cols=c("#D68747","#BEABD2","#754D40"),
                        pt.size = 0.01,
                        group.by = "celltype",label = T,label.size=7)


####umap plot
umap = sce.all.0.01@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = sce.all.0.01@meta.data$celltype) 
library(ggpubr)
umap$cell_type <- factor(umap$cell_type,levels=c("Stromal","Immune","Epithelial"))
p1 <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +
  geom_point(size = 0.5 , shape=16, stroke=0 )  +
  scale_color_manual(values = c("#D68747","#BEABD2","#754D40"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(), 
    legend.key=element_rect(fill='white'), 
    legend.text = element_text(size=9),
    legend.key.size=unit(1.2,'line') ,
    legend.position = 'right') +  
  guides(color = guide_legend(override.aes = list(size=5,alpha=1))) +
  geom_segment(aes(x = min(umap$umap_1) , y = min(umap$umap_2) ,
                   xend = min(umap$umap_1) +3, yend = min(umap$umap_2) ),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$umap_1)  , y = min(umap$umap_2)  ,
                   xend = min(umap$umap_1) , yend = min(umap$umap_2) + 3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$umap_1) +1.5, y = min(umap$umap_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$umap_1) -1, y = min(umap$umap_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 

umap_group = sce.all.0.01@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(group = sce.all.0.01@meta.data$group) 

sce.all.0.01@meta.data$group <- 
  
  p2 <- ggplot(umap_group,aes(x= umap_1 , y = umap_2 ,color = group)) +
  geom_point(size = 0.5 , shape=16, stroke=0 )  +
  scale_color_manual(values = c("#C14231","#aaaaaa","#417ab5"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(), 
    legend.key=element_rect(fill='white'), 
    legend.text = element_text(size=9), 
    legend.key.size=unit(1.2,'line') ,
    legend.position = 'right') +  
  guides(color = guide_legend(override.aes = list(size=5,alpha=1))) + 
  geom_segment(aes(x = min(umap$umap_1) , y = min(umap$umap_2) ,
                   xend = min(umap$umap_1) +3, yend = min(umap$umap_2) ),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$umap_1)  , y = min(umap$umap_2)  ,
                   xend = min(umap$umap_1) , yend = min(umap$umap_2) + 3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$umap_1) +1.5, y = min(umap$umap_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$umap_1) -1, y = min(umap$umap_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
p2

#####marker gene heatmap
mks = sce.all.0.01.cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(subset(sce.all.0.01,downsample=100),mks$gene,size=3)+
  scale_fill_gradientn(colors=colorRampPalette(c("#aaaaaa","#ffffff", "#C14231"))(100))

########subcluster annotation #########

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
#install.packages("devtools")
#devtools::install_github("junjunlab/jjAnno")
library(enrichplot)
library(ggstatsplot)
library(Vennerable) 
library(ggplot2)
library(pheatmap)
library(cowplot)
library(ggrepel)
library(clusterProfiler) #BiocManager::install("org.Rn.eg.db")

library(org.Rn.eg.db)
library(DOSE)
library(pathview)


sce = sce.all.1[, sce.all.1$RNA_snn_res.1 %in% c( '2','4','10','15',
                                                  '18','24','27','31')]
Idents(sce) = sce$orig.ident   
sce_counts <- LayerData(sce,assay = "RNA",layer = "counts")
metadata <- sce.all.1@meta.data[c(1:8,19)]

sce=CreateSeuratObject(counts = sce_counts,
                       meta.data = metadata)

sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce))
DimHeatmap(sce, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2))
clustree(sce@meta.data, prefix = "RNA_snn_res.")

sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)
sce <- RunUMAP(object = sce, dims = 1:15, do.fast = TRUE)

seuratObj <- RunHarmony(sce, "orig.ident")
seuratObj <- RunUMAP(seuratObj, dims = 1:8, reduction = "harmony")
sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:5)
sce <- FindClusters(sce, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1,1.2))
sce.all.0.01 <- SetIdent(sce.all,value="RNA_snn_res.0.01")

sce.cluster.markers <-  FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sce.cluster.markers,file="immune.all.markers.csv")
top10pos <- sce.cluster.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top10negtive <- sce.cluster.markers %>% group_by(cluster) %>% top_n(n = -10, wt = avg_log2FC)
top_bottom_gene <-data.frame(top10pos,top10negtive)

p1 <- DimPlot(sce,reduction = 'umap',pt.size = .1,label = T,group.by = "group",
              label.size = 4,raster=FALSE)
p2 <- DimPlot(sce,reduction = 'umap',pt.size = .1,label = T,
              label.size = 4,raster=FALSE)
p1+p2


genes_to_check_all_immune_subcluster <- c(
  #  "RT1-Da",'RT1-Db1','RT1-Db2','RT1-DOa','RT1-DOb',           
  'RT1-DMa','RT1-DMb',"RT1-Bb", "RT1-Ba",                     
  "Cpa3","Cma1",                                               
  'Cd79a', 'Ms4a1',                                           
  "Ighm",'Vpreb3',"Mki67",                                    
  'Mzb1','Ighg1',                                              
  'Bank1','Pax5',                                            
  "Flt3",'Itgax','Itgam',                                     
  'Zeb2',                                                      
  'Cd24','Siglech','Tcf4',                                         
  'Ly6c',"Ccr2",                                               
  'Cd68','C1qb','C1qa', "Csf1r",                       
  "Mrc1","Cd163",                                              
  "Mmp12",  "Il10",    	                                         
  'Spp1',"Tnf", "Il1b",                                        
  'Nos2',  'Cx3cr1',                                           
  'S100a8', 'S100a9', 'Retnlg', 'Mmp8',                        
  'Cxcr2', 'Csf1',                                             
  'Cd3d','Cd3e', 'ENSRNOG00000071219',                             
  'Cd8a','Cd8b',
  'Foxp3','Il2ra','Ctla4','Tnfrsf18','Tnfrsf4',          
  #'Ctsw','Klrd1','Ccl5',"Xcl1",
  'Ikzf2','Il10','Tgfb1','Tnf', 'Stat5a',                 
  'Gzmm','Pdcd1',
  'Cd3e', 'ENSRNOG00000071219',                                
  'Cd8a','Cd8b','Lck','Rora','Il17a',
  'Ctsw','Gzmb','Klrd1','Ccl5',"Xcl1",
  'Slc9a9','Inpp4b','Ikzf2',                            
  'Foxp3','Il2ra','Ctla4','Tnfrsf18','Tnfrsf4'          
)

p <- DotPlot(sce, features = unique(genes_to_check_all_immune_subcluster),
             assay='RNA'  )  +
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45,size=7))
p
p$data$id <- factor(p$data$id,levels = c(
  '11',
  '13','8','15','10',
  '16','5','3',
  '7','9','6','17','12',
  '14','4',
  '1','0','2'
))


p_immune <-  ggplot(p$data,aes(x = features.plot,y = id)) +  
  geom_point(aes(fill = avg.exp.scaled,size = pct.exp),
             color='black',shape=21) +  
  theme_bw(base_size = 14) +  xlab('') + ylab('') +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5,                     
                       limits = c(-2,3),                    
                       breaks = seq(-2,3,1),        
                       labels = seq(-2,3, 1),   
                       name = 'Mean\nexpression') + 
  scale_size(range = c(0,6),            
             limits = c(0, 100),          
             breaks = seq(20,100,20),       
             labels = seq(20,100,20)  ) + 
  theme(panel.grid = element_blank(),    
        axis.text = element_text(color = 'black'),    
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = 'italic') 
  ) +   
  guides(  
    fill = guide_colorbar(    
      direction = "vertical",   
      title.position = "top",    
      barwidth = unit(0.5, "cm")  
    ), 
    size = guide_legend(     
      title = "Percent\nof cells",      
      direction = "vertical",      
      title.position = "top",     
      label.position = "right",     
      override.aes = list(      
        color = "black",      
        fill = "grey"     
      ) ) )
p_immune

VlnPlot(sce,features=c("Nfe2l2",'Ptgs2','Nos2','Bcl2','Nfkb1','Nfkbia','Tnf','Il6'), pt.size = 0) 
VlnPlot(sce,features=c('Nfkb1','Nfkbia'), group.by = "celltype", pt.size = 0)


########
celltype=data.frame(ClusterID=0:17,
                    celltype= 0:17) 

celltype[celltype$ClusterID %in% c( 4 ),2]='Neu_Retnlg'
celltype[celltype$ClusterID %in% c( 14 ),2]='Neu_Csf1'
celltype[celltype$ClusterID %in% c( 11 ),2]='Mast'
celltype[celltype$ClusterID %in% c( 8 ),2]='IgG_plasma_B'
celltype[celltype$ClusterID %in% c( 13 ),2]='Memory_B'
celltype[celltype$ClusterID %in% c( 10,15 ),2]='Plasmablast'
celltype[celltype$ClusterID %in% c( 0 ),2]='T_Cd8_Fopx3_Il2ra'
celltype[celltype$ClusterID %in% c( 1 ),2]='T_Cd8_Foxp3_Gzmm'
celltype[celltype$ClusterID %in% c( 2 ),2]='Treg_Cd4_Ikzf2'
celltype[celltype$ClusterID %in% c( 3 ),2]='cDC1'
celltype[celltype$ClusterID %in% c( 5 ),2]='cDC2_Zeb2'
celltype[celltype$ClusterID %in% c( 16 ),2]='pDCs'
celltype[celltype$ClusterID %in% c( 7,9 ),2]='Macs_Mrc1_Cd163'
celltype[celltype$ClusterID %in% c( 6 ),2]='Macs_Itgax'
celltype[celltype$ClusterID %in% c( 17 ),2]='Macs_Spp1'
celltype[celltype$ClusterID %in% c( 12 ),2]='Ly6c_Monocytes'

sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)
DimPlot(sce,reduction = 'umap',pt.size = .1,label = T,group.by = "celltype",
        label.size = 2.5,raster=FALSE)

p <- DotPlot(sce, features = unique(genes_to_check_all_immune_subcluster),
             group.by = "celltype",assay='RNA'  )
p$data$id <- factor(p$data$id,levels = c(
  'Mast','IgG_plasma_B','Memory_B','Plasmablast',
  'pDCs','cDC2_Zeb2','cDC1',
  'Ly6c_Monocytes','Macs_Spp1','Macs_Itgax','Macs_Mrc1_Cd163',
  'Neu_Csf1','Neu_Retnlg',
  'Treg_Cd4_Ikzf2','T_Cd8_Foxp3_Gzmm','T_Cd8_Fopx3_Il2ra'
))

p_immune <-  ggplot(p$data,aes(x = features.plot,y = id)) +  
  geom_point(aes(fill = avg.exp.scaled,size = pct.exp),
             color='black',shape=21) +  
  theme_bw(base_size = 14) +  xlab('') + ylab('') +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5,                     
                       limits = c(-2,3),                    
                       breaks = seq(-2,3,1),        
                       labels = seq(-2,3, 1),   
                       name = 'Mean\nexpression') + 
  scale_size(range = c(0,6),            
             limits = c(0, 100),          
             breaks = seq(20,100,20),       
             labels = seq(20,100,20)  ) + 
  theme(panel.grid = element_blank(),    
        axis.text = element_text(color = 'black'),    
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = 'italic') 
  ) +   
  guides(  
    fill = guide_colorbar(    
      direction = "vertical",   
      title.position = "top",    
      barwidth = unit(0.5, "cm")  
    ), 
    size = guide_legend(     
      title = "Percent\nof cells",      
      direction = "vertical",      
      title.position = "top",     
      label.position = "right",     
      override.aes = list(      
        color = "black",      
        fill = "grey"     
      ) ) )
p_immune


library(ggpubr)

umap = sce@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = sce@meta.data$celltype) 
mycolors_immune <- c(
  "#9970AB","#E08214","#ee6aa7","#A6DBA0",
  "#5E1415","#417ab5","#F1B6DA",
  "#1B7837","#762A83","#B32B34","#684A11", 
  "#CC4C02","#7A5399",
  "#542788","#FEE0B6","#AEC42F"
)

umap$cell_type <- factor(umap$cell_type,
                         levels=c('Mast','IgG_plasma_B','Memory_B','Plasmablast',
                                  'pDCs','cDC2_Zeb2','cDC1',
                                  'Macs_Mrc1_Cd163','Macs_Itgax','Macs_Spp1','Macs_Nos2',
                                  'Neu_Csf1','Neu_Retnlg',
                                  'Treg_Cd4_Ikzf2','T_Cd8_Foxp3_Gzmm','T_Cd8_Fopx3_Il2ra'	))
p1 <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +
  geom_point(size = 0.5 , shape=16, stroke=0 )  +
  #  theme_pubr()+
  scale_color_manual(values = mycolors_immune)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(), 
    legend.key=element_rect(fill='white'), 
    legend.text = element_text(size=9), 
    legend.key.size=unit(1.2,'line') ,
    legend.position = 'right') +  
  guides(color = guide_legend(override.aes = list(size=5,alpha=1)))  
p1


library(readxl)
library(ggplot2)
library(ggh4x)
library(tidyr)
library(reshape2)
library(gplots)

phe=sce@meta.data

colnames(phe)
tb=table(phe$celltype, phe$orig.ident)
head(tb)

balloonplot(tb)
bar_data<-as.data.frame(tb)
bar_per<-bar_data%>%
  group_by(Var2)%>%
  mutate(sum(Freq))%>%
  mutate(percent=Freq/`sum(Freq)`)

tb2=table(phe$celltype, phe$group)
head(tb2)
balloonplot(tb2)
bar_data2<-as.data.frame(tb2)

bar_per2<-bar_data2%>%
  group_by(Var2)%>%
  mutate(sum(Freq))%>%
  mutate(percent=Freq/`sum(Freq)`)
head(bar_per2)
write.csv(bar_per2,file="immune_subcluster_by_group_percent.csv")

bar_per2$Var2 <- factor(bar_per2$Var2,levels=c("Control","CD","PB_IFX"))
#排序
bar_per2$Var1 <- factor(bar_per2$Var1,
                        levels=c('Mast','IgG_plasma_B','Memory_B','Plasmablast',
                                 'pDCs','cDC2_Zeb2','cDC1',
                                 'Macs_Mrc1_Cd163','Macs_Itgax','Macs_Spp1','Macs_Nos2',
                                 'Neu_Csf1','Neu_Retnlg',
                                 'Treg_Cd4_Ikzf2','T_Cd8_Foxp3_Gzmm','T_Cd8_Fopx3_Il2ra'))

p1 <- ggplot(bar_per2,aes(x=Var1,y=percent,color=Var2, fill=Var2))+
  geom_bar(stat="identity",position = position_dodge(0.9),width=0.75,alpha=0.8)+
  scale_y_continuous(expand = c(0,0),
                     limits = c(-0.01,0.4))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size = 10,
                                   face="bold" ,angle=45, vjust=1,hjust=1),
        axis.text.y = element_text(colour = "black", size = 8),
        legend.position="top"
  )+
  labs(x="", y="Fraction of total cells")+
  scale_color_manual(values=c("#9B989B","#BE0F20","#272B82"))+
  scale_fill_manual(values=c("#9B989B","#BE0F20","#272B82"))
p1


##### 
Idents(sce) = sce$group
table(Idents(sce))
sce.1  = sce[, sce$celltype %in% c( 'Mast','IgG_plasma_B','Memory_B','Plasmablast',
                                    'pDCs','cDC2_Zeb2','cDC1',
                                    'Macs_Mrc1_Cd163','Macs_Itgax','Macs_Spp1','Macs_Nos2',
                                    'Neu_Retnlg',
                                    'Treg_Cd4_Ikzf2','T_Cd8_Foxp3_Gzmm','T_Cd8_Fopx3_Il2ra')]

degs = lapply(unique(sce.1$celltype), function(x){
  FindMarkers(sce.1[,sce.1$celltype==x],ident.1 = 'PB_IFX',
              ident.2 = 'CD')
})
names(degs) <- unique(sce.1$celltype)

for (i in 1:length(degs)){
  dat <- degs[[i]]
  dat$cluster <- rep(names(degs)[[i]],times=nrow(dat))
  dat$gene <- rownames(dat)
  filename1= paste0(unique(dat$cluster),".clusterdeg.csv")
  write.csv(dat,file=filename1)
}

##### 

cell <- c( 'Mast','B_Tnfrsf13c','B_Bank1','B_Mki67',
           'pDCs','cDC2_Zeb2','cDC1',
           'Macs_Mrc1_Cd163','Macs_Itgax','Macs_Spp1','Macs_Nos2',
           'Neu_Retnlg','Treg','T')

for (i in 1:1) {
  filename1= paste0(unique(cell[i]),".clusterdeg.csv")
  DEG = read.csv(filename1,row.names = 1)
  gene_up <- rownames(DEG[with(DEG,avg_log2FC>0.1 & p_val_adj<0.05),])
  gene_down <- rownames(DEG[with(DEG,avg_log2FC< -0.1 & p_val_adj<0.05),])
  gene_up_entrez <- as.character(na.omit(bitr(gene_up, 
                                              fromType="SYMBOL", 
                                              toType="ENTREZID", 
                                              OrgDb="org.Rn.eg.db")[,2])) 
  
  kegg_enrich_up_results <- enrichKEGG(gene  = gene_up_entrez,
                                       organism  = "rno" , 
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff = 0.2
  )
  kegg_enrich_up_results <- DOSE::setReadable(kegg_enrich_up_results, 
                                              OrgDb="org.Rn.eg.db", 
                                              keyType='ENTREZID')
  filename2= paste0(unique(cell[i]),".KEGG.up.list.csv")
  gene_down_entrez <- as.character(na.omit(bitr(gene_down, 
                                                fromType="SYMBOL", 
                                                toType="ENTREZID",
                                                OrgDb="org.Rn.eg.db")[,2])) 
  
  kegg_enrich_down_results <- enrichKEGG(gene  = gene_down_entrez,
                                         organism  = "rno" , 
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.2
  )
  kegg_enrich_down_results <- DOSE::setReadable(kegg_enrich_down_results, 
                                                OrgDb="org.Rn.eg.db", 
                                                keyType='ENTREZID')
}


####
sce = sce.all.1[, sce.all.1$RNA_snn_res.1 %in% c( '1','8','11','17',
                                                  '20','22','26')]
Idents(sce) = sce$orig.ident    
sce_counts <- LayerData(sce,assay = "RNA",layer = "counts")
metadata <- sce.all.1@meta.data[c(1:8,19)]

sce=CreateSeuratObject(counts = sce_counts,
                       meta.data = metadata)

sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce))
DimHeatmap(sce, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1))
clustree(sce@meta.data, prefix = "RNA_snn_res.")
setwd("stromal")
ggsave("clustree.pdf")
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)
sce <- RunUMAP(object = sce, dims = 1:15, do.fast = TRUE)

seuratObj <- RunHarmony(sce, "orig.ident")
seuratObj <- RunUMAP(seuratObj, dims = 1:8, reduction = "harmony")
sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:5)
sce <- FindClusters(sce, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1))

save(sce, file='stromal.counts.Rdata')
load("stromal.counts.Rdata")

sce.cluster.markers <-  FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sce.cluster.markers,file="stromal.all.markers.csv")
top10pos <- sce.cluster.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC) 
top10negtive <- sce.cluster.markers %>% group_by(cluster) %>% top_n(n = -10, wt = avg_log2FC)
top_bottom_gene <-data.frame(top10pos,top10negtive)
write.csv(top_bottom_gene,file="stromal.top_bottom_gene.markers.csv")

p1 <- DimPlot(sce,reduction = 'umap',pt.size = .1,label = T,group.by = "group",
              label.size = 4,raster=FALSE)
p2 <- DimPlot(sce,reduction = 'umap',pt.size = .1,label = T,
              label.size = 4,raster=FALSE)
p1+p2

######
genes_to_check_all_stroma_subcluster <- c(
  'Vwf','Podxl', 'RGD1565355', 'Plvap',              
  'Ackr1','Sele', 'Il33',                           
  'Ccl21', 'Mmrn1',                                  
  'Pdgfra','Dpt','Adamdec1',                         
  'Fap','Il11', 'Mmp3', 'Cxcl14',                    
  'Ccl11', 'Kcnn3','Smoc2',                         
  'Mfap5', 'Gsn',                                    
  'Plp1','Sox10','S100b', 'Nrxn1',                   
  'Acta2','Actg2','Tagln',                           
  'Des',                                             
  'Myh11',                                       
  'Mki67','Hmgb2','Ube2c','Pttg1','Top2a','Nusap1'   
)

p <- DotPlot(sce, features = unique(genes_to_check_all_stroma_subcluster),
             assay='RNA'  )  +
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45,size=7))

p$data$id <- factor(p$data$id,levels = c(
  '1','8',
  '7',
  '0','2','3','4','11','16','5','6','10',
  '13',
  '9','12','15',
  '14'
))


p_stroma <-  ggplot(p$data,aes(x = features.plot,y = id)) +  
  geom_point(aes(fill = avg.exp.scaled,size = pct.exp),
             color='black',shape=21) +  
  theme_bw(base_size = 14) +  xlab('') + ylab('') +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5,                     
                       limits = c(-2,3),                    
                       breaks = seq(-2,3,1),        
                       labels = seq(-2,3, 1),   
                       name = 'Mean\nexpression') + 
  scale_size(range = c(0,6),            
             limits = c(0, 100),          
             breaks = seq(20,100,20),       
             labels = seq(20,100,20)  ) + 
  theme(panel.grid = element_blank(),    
        axis.text = element_text(color = 'black'),    
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = 'italic') 
  ) +   
  guides(  
    fill = guide_colorbar(    
      direction = "vertical",   
      title.position = "top",    
      barwidth = unit(0.5, "cm")  
    ), 
    size = guide_legend(     
      title = "Percent\nof cells",      
      direction = "vertical",      
      title.position = "top",     
      label.position = "right",     
      override.aes = list(      
        color = "black",      
        fill = "grey"     
      ) ) )
p_stroma

#####
celltype=data.frame(ClusterID=0:16,
                    celltype= 0:16) 

celltype[celltype$ClusterID %in% c( 1,8 ),2]='Endothelial'
celltype[celltype$ClusterID %in% c( 7 ),2]='Lymphatics'
celltype[celltype$ClusterID %in% c( 6,10 ),2]='Fibro_IL11_Cxcl14'
celltype[celltype$ClusterID %in% c( 5 ),2]='Fibro_Smoc2_Ccl11'
celltype[celltype$ClusterID %in% c( 0,2,3,4,11,16 ),2]='Fibro_Mfap5_Gsn'
celltype[celltype$ClusterID %in% c( 13 ),2]='Glial'
celltype[celltype$ClusterID %in% c( 9,12,15 ),2]='Myofibroblasts'
celltype[celltype$ClusterID %in% c( 14 ),2]='Cycling_stroma'


sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)
save(sce,file="stromal.annotated.Rdata")
load("stromal.annotated.Rdata")

#####
p <- DotPlot(sce, features = unique(genes_to_check_all_stroma_subcluster),
             group.by = "celltype",assay='RNA'  )
p$data$id <- factor(p$data$id,levels = c(
  'Endothelial','Lymphatics','Fibro_IL11_Cxcl14','Fibro_Smoc2_Ccl11',
  'Fibro_Mfap5_Gsn','Glial','Myofibroblasts','Cycling_stroma'
))


p_stroma <-  ggplot(p$data,aes(x = features.plot,y = id)) +  
  geom_point(aes(fill = avg.exp.scaled,size = pct.exp),
             color='black',shape=21) +  
  theme_bw(base_size = 14) +  xlab('') + ylab('') +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5,                     
                       limits = c(-2,3),                    
                       breaks = seq(-2,3,1),        
                       labels = seq(-2,3, 1),   
                       name = 'Mean\nexpression') + 
  scale_size(range = c(0,6),            
             limits = c(0, 100),          
             breaks = seq(20,100,20),       
             labels = seq(20,100,20)  ) + 
  theme(panel.grid = element_blank(),    
        axis.text = element_text(color = 'black'),    
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = 'italic') 
  ) +   
  guides(  
    fill = guide_colorbar(    
      direction = "vertical",   
      title.position = "top",    
      barwidth = unit(0.5, "cm")  
    ), 
    size = guide_legend(     
      title = "Percent\nof cells",      
      direction = "vertical",      
      title.position = "top",     
      label.position = "right",     
      override.aes = list(      
        color = "black",      
        fill = "grey"     
      ) ) )
p_stroma
ggsave('stromal.dotplot.pdf')

##### 
library(ggpubr)

umap = sce@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = sce@meta.data$celltype) 
mycolors_stromal <- c(
  "#FE9929","#684A11", "#B32B34", "#762A83","#993404","#F4E726" , "#754555","#5AAE61" 
)

umap$cell_type <- factor(umap$cell_type,
                         levels=c( 'Endothelial','Lymphatics',
                                   'Fibro_IL11_Cxcl14','Fibro_Smoc2_Ccl11',
                                   'Fibro_Mfap5_Gsn','Glial','Myofibroblasts',
                                   'Cycling_stroma'	))
p1 <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +
  geom_point(size = 0.5 , shape=16, stroke=0 )  +
  #  theme_pubr()+
  scale_color_manual(values = mycolors_stromal)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(), 
    legend.key=element_rect(fill='white'), 
    legend.text = element_text(size=9), 
    legend.key.size=unit(1.2,'line') ,
    legend.position = 'right') +  
  guides(color = guide_legend(override.aes = list(size=5,alpha=1)))   
p1

##### 
library(readxl)
library(ggplot2)
library(ggh4x)
library(tidyr)
library(reshape2)
library(gplots)

phe=sce@meta.data

colnames(phe)
tb=table(phe$celltype, phe$orig.ident)
head(tb)

balloonplot(tb)
bar_data<-as.data.frame(tb)
bar_per<-bar_data%>%
  group_by(Var2)%>%
  mutate(sum(Freq))%>%
  mutate(percent=Freq/`sum(Freq)`)
head(bar_per)
write.csv(bar_per,file="stromal_subcluster_by_sample_percent.csv")
tb2=table(phe$celltype, phe$group)
head(tb2)
balloonplot(tb2)
bar_data2<-as.data.frame(tb2)

bar_per2<-bar_data2%>%
  group_by(Var2)%>%
  mutate(sum(Freq))%>%
  mutate(percent=Freq/`sum(Freq)`)
head(bar_per2)
write.csv(bar_per2,file="stromal_subcluster_by_group_percent.csv")

bar_per2$Var2 <- factor(bar_per2$Var2,levels=c("Control","CD","PB_IFX"))

bar_per2$Var1 <- factor(bar_per2$Var1,
                        levels=c('Endothelial','Lymphatics','Fibro_IL11_Cxcl14','Fibro_Smoc2_Ccl11',
                                 'Fibro_Mfap5_Gsn','Glial','Myofibroblasts','Cycling_stroma'))

p1 <- ggplot(bar_per2,aes(x=Var1,y=percent,color=Var2, fill=Var2))+
  geom_bar(stat="identity",position = position_dodge(0.9),width=0.75,alpha=0.8)+
  scale_y_continuous(expand = c(0,0),
                     limits = c(-0.01,0.8))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size = 10,
                                   face="bold" ,angle=45, vjust=1,hjust=1),
        axis.text.y = element_text(colour = "black", size = 8),
        legend.position="top"
  )+
  labs(x="", y="Fraction of total cells")+
  scale_color_manual(values=c("#9B989B","#BE0F20","#272B82"))+
  scale_fill_manual(values=c("#9B989B","#BE0F20","#272B82"))
p1

#####
Idents(sce) = sce$group
table(Idents(sce))
sce.1  = sce[, sce$celltype %in% c( 'Endothelial','Lymphatics','Fibro_IL11_Cxcl14','Fibro_Smoc2_Ccl11',
                                    'Fibro_Mfap5_Gsn','Glial','Myofibroblasts')]

degs = lapply(unique(sce.1$celltype), function(x){
  FindMarkers(sce.1[,sce.1$celltype==x],ident.1 = 'PB_IFX',
              ident.2 = 'CD')
})
names(degs) <- unique(sce.1$celltype)
for (i in 1:length(degs)){
  dat <- degs[[i]]
  dat$cluster <- rep(names(degs)[[i]],times=nrow(dat))
  dat$gene <- rownames(dat)
  filename1= paste0(unique(dat$cluster),".clusterdeg.csv")
  write.csv(dat,file=filename1)
}

up <- data.frame()
down <- data.frame()
for (i in 1:length(degs)){
  dat <- degs[[i]]
  dat$cell_type <- rep(names(degs)[[i]],times=nrow(dat))
  dat$gene <- rownames(dat)
  up_regulate_df <- dat[dat$avg_log2FC>0 & dat$p_val_adj<0.05,c('avg_log2FC','cell_type','gene')]
  
  down_regulate_df <- dat[dat$avg_log2FC<0 & dat$p_val_adj<0.05,c('avg_log2FC','cell_type','gene')]
  up <- rbind(up,up_regulate_df)
  down <- rbind(down,down_regulate_df)
}

#####

cell <- c( 'Endothelial','Fibro_IL11_Cxcl14','Fibro_Smoc2_Ccl11',
           'Fibro_Mfap5_Gsn','Myofibroblasts')

for (i in 1:5) {
  filename1= paste0(unique(cell[i]),".clusterdeg.csv")
  DEG = read.csv(filename1,row.names = 1)
  gene_up <- rownames(DEG[with(DEG,avg_log2FC>0.1 & p_val_adj<0.05),])
  gene_down <- rownames(DEG[with(DEG,avg_log2FC< -0.1 & p_val_adj<0.05),])
  gene_up_entrez <- as.character(na.omit(bitr(gene_up,
                                              fromType="SYMBOL", 
                                              toType="ENTREZID", 
                                              OrgDb="org.Rn.eg.db")[,2])) 
  
  kegg_enrich_up_results <- enrichKEGG(gene  = gene_up_entrez,
                                       organism  = "rno" , 
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff = 0.2
  )
  kegg_enrich_up_results <- DOSE::setReadable(kegg_enrich_up_results, 
                                              OrgDb="org.Rn.eg.db", 
                                              keyType='ENTREZID')
  filename2= paste0(unique(cell[i]),".KEGG.up.list.csv")
  write.csv(kegg_enrich_up_results@result,filename2)
  
  gene_down_entrez <- as.character(na.omit(bitr(gene_down, 
                                                fromType="SYMBOL", 
                                                toType="ENTREZID", 
                                                OrgDb="org.Rn.eg.db")[,2])) 
  kegg_enrich_down_results <- enrichKEGG(gene  = gene_down_entrez,
                                         organism  = "rno" , 
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.2
  )
  kegg_enrich_down_results <- DOSE::setReadable(kegg_enrich_down_results, 
                                                OrgDb="org.Rn.eg.db", 
                                                keyType='ENTREZID')
  filename3= paste0(unique(cell[i]),".KEGG.down.list.csv")
  write.csv(kegg_enrich_down_results@result,filename3)
}


###
sce = sce.all.1[, sce.all.1$RNA_snn_res.1 %in% c( '0','5','6','7','25',
                                                  '3','32','12','21','23',
                                                  '16','19','9','30','34','33','13','14',
                                                  '28','29')]
Idents(sce) = sce$orig.ident     
sce_counts <- LayerData(sce,assay = "RNA",layer = "counts")
metadata <- sce.all.1@meta.data[c(1:8,19)]

sce=CreateSeuratObject(counts = sce_counts,
                       meta.data = metadata)

sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce))
DimHeatmap(sce, dims = 19:30, cells = 100, balanced = TRUE)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1))
clustree(sce@meta.data, prefix = "RNA_snn_res.")
setwd("epithelium")
ggsave("clustree.pdf")
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)
sce <- RunUMAP(object = sce, dims = 1:15, do.fast = TRUE)

seuratObj <- RunHarmony(sce, "orig.ident")
seuratObj <- RunUMAP(seuratObj, dims = 1:8, reduction = "harmony")
sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:5)
sce <- FindClusters(sce, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1))

sce.cluster.markers <-  FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sce.cluster.markers,file="epithelium.all.markers.csv")
top10pos <- sce.cluster.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC) 
top10negtive <- sce.cluster.markers %>% group_by(cluster) %>% top_n(n = -10, wt = avg_log2FC)
top_bottom_gene <-data.frame(top10pos,top10negtive)
write.csv(top_bottom_gene,file="epithelium.top_bottom_gene.markers.csv")

p1 <- DimPlot(sce,reduction = 'umap',pt.size = .1,label = T,group.by = "group",
              label.size = 4,raster=FALSE)
p2 <- DimPlot(sce,reduction = 'umap',pt.size = .1,label = T,
              label.size = 4,raster=FALSE)
p1+p2

###### 
genes_to_check_epithelium_subcluster = c(
  "Epcam",
  'Lgr5', 'Olfm4', "Satb2",'Rgmb',                       
  'Chga',"Chgb",'Scg5','Scgn',                           
  'Fev', 'Neurod1',                                      
  'Lix1',                                                
  'Car1','Car2', 'Selenop','Selenbp1',                   
  "Fabp2", "Rbp2", 'Alpi',                               
  'Mep1a',                                               
  'Sult1a1',                                             
  'Duoxa2','Tff1', 'S100a9',                            
  'Muc2', "Spdef",                                                   
  'Spink4', 'Itln1', 'Clca1',                           
  'Elapor1','Gmds',                                     
  'Wfdc2','Best2',                                       
  'Cenpa',                                               
  'Scgb1a1','Sftpd',                                     
  'Msmb', 'Pbsn',                                       
  'Sh2d6', 'Irag2',"Trpm5",                              
  'Nusap1',"Ube2c","Hmgb2","Top2a","Mki67"              
)

p <- DotPlot(sce, features = unique(genes_to_check_epithelium_subcluster),
             assay='RNA'  )  +
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45,size=7))

p$data$id <- factor(p$data$id,levels = c(
  '6',
  '19','20',
  '12','7','15','9','13','14','21',
  '0','1','4','10','17',
  '18','3',
  '16',
  '2','11',
  '5','8'
))


p_epithelium <-  ggplot(p$data,aes(x = features.plot,y = id)) +  
  geom_point(aes(fill = avg.exp.scaled,size = pct.exp),
             color='black',shape=21) +  
  theme_bw(base_size = 14) +  xlab('') + ylab('') +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5,                     
                       limits = c(-2,3),                    
                       breaks = seq(-2,3,1),        
                       labels = seq(-2,3, 1),   
                       name = 'Mean\nexpression') + 
  scale_size(range = c(0,6),            
             limits = c(0, 100),          
             breaks = seq(20,100,20),       
             labels = seq(20,100,20)  ) + 
  theme(panel.grid = element_blank(),    
        axis.text = element_text(color = 'black'),    
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = 'italic') 
  ) +   
  guides(  
    fill = guide_colorbar(    
      direction = "vertical",   
      title.position = "top",    
      barwidth = unit(0.5, "cm")  
    ), 
    size = guide_legend(     
      title = "Percent\nof cells",      
      direction = "vertical",      
      title.position = "top",     
      label.position = "right",     
      override.aes = list(      
        color = "black",      
        fill = "grey"     
      ) ) )
p_epithelium

###### 

celltype=data.frame(ClusterID=0:21,
                    celltype= 0:21) 

celltype[celltype$ClusterID %in% c( 6 ),2]='Stem'
celltype[celltype$ClusterID %in% c( 19 ),2]='EECs_Neurod1'
celltype[celltype$ClusterID %in% c( 20 ),2]='EECs_Lix1'
celltype[celltype$ClusterID %in% c( 7,12 ),2]='Entero_Rbp2'
celltype[celltype$ClusterID %in% c( 15 ),2]='Undiff_Entero'
celltype[celltype$ClusterID %in% c( 13 ),2]='Entero_Mep1a'
celltype[celltype$ClusterID %in% c( 9 ),2]='Entero_Sult1a1'
celltype[celltype$ClusterID %in% c( 14,21 ),2]='Entero_Duoxa2_s100a9'
celltype[celltype$ClusterID %in% c( 0,1,4 ),2]='Goblet_Best2'
celltype[celltype$ClusterID %in% c( 10 ),2]='Goblet_Cenpa'
celltype[celltype$ClusterID %in% c( 17 ),2]='Goblet_Gmds'
celltype[celltype$ClusterID %in% c( 3,18 ),2]='Goblet_pre'
celltype[celltype$ClusterID %in% c( 16 ),2]='Mucous_cells'
celltype[celltype$ClusterID %in% c( 2,11 ),2]='Tuft'
celltype[celltype$ClusterID %in% c( 5,8 ),2]='TA'

sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)

##### 
p <- DotPlot(sce, features = unique(genes_to_check_epithelium_subcluster),
             group.by = "celltype",assay='RNA'  )
p$data$id <- factor(p$data$id,levels = c(
  'Stem','EECs_Neurod1','EECs_Lix1','Undiff_Entero','Entero_Rbp2','Entero_Mep1a',
  'Entero_Sult1a1','Entero_Duoxa2_s100a9','Goblet_Best2','Goblet_Cenpa',
  'Goblet_Gmds','Goblet_pre','Mucous_cells','Tuft','TA'
))


p_epithelim <-  ggplot(p$data,aes(x = features.plot,y = id)) +  
  geom_point(aes(fill = avg.exp.scaled,size = pct.exp),
             color='black',shape=21) +  
  theme_bw(base_size = 14) +  xlab('') + ylab('') +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5,                     
                       limits = c(-2,3),                    
                       breaks = seq(-2,3,1),        
                       labels = seq(-2,3, 1),   
                       name = 'Mean\nexpression') + 
  scale_size(range = c(0,6),            
             limits = c(0, 100),          
             breaks = seq(20,100,20),       
             labels = seq(20,100,20)  ) + 
  theme(panel.grid = element_blank(),    
        axis.text = element_text(color = 'black'),    
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = 'italic') 
  ) +   
  guides(  
    fill = guide_colorbar(    
      direction = "vertical",   
      title.position = "top",    
      barwidth = unit(0.5, "cm")  
    ), 
    size = guide_legend(     
      title = "Percent\nof cells",      
      direction = "vertical",      
      title.position = "top",     
      label.position = "right",     
      override.aes = list(      
        color = "black",      
        fill = "grey"     
      ) ) )
p_epithelim

###### 

library(ggpubr)

umap = sce@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = sce@meta.data$celltype) 
mycolors_epithelium <- c(
  "#4D9221","#543005","#B7794C","#CC4C02" ,
  "#D79D9D","#c4e6c3","#36877a","#4E1365",
  "#9A5C18","#A29F9D",
  "#8189C5","#003C30","#EDCF46","#E1DF97","#576B26"
)
"#754D40"
umap$cell_type <- factor(umap$cell_type,
                         levels=c( 'Stem','EECs_Neurod1','EECs_Lix1','Undiff_Entero','Entero_Rbp2','Entero_Mep1a',
                                   'Entero_Sult1a1','Entero_Duoxa2_s100a9','Goblet_Best2','Goblet_Cenpa',
                                   'Goblet_Gmds','Goblet_pre','Mucous_cells','Tuft','TA'	))
p1 <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +
  geom_point(size = 0.5 , shape=16, stroke=0 )  +
  scale_color_manual(values = mycolors_epithelium )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(), 
    legend.key=element_rect(fill='white'), 
    legend.text = element_text(size=9), 
    legend.key.size=unit(1.2,'line') ,
    legend.position = 'right') +  
  guides(color = guide_legend(override.aes = list(size=5,alpha=1))) 
p1

##### 
library(readxl)
library(ggplot2)
library(ggh4x)
library(tidyr)
library(reshape2)
library(gplots)

phe=sce@meta.data

colnames(phe)
tb=table(phe$celltype, phe$orig.ident)
head(tb)

balloonplot(tb)
bar_data<-as.data.frame(tb)
bar_per<-bar_data%>%
  group_by(Var2)%>%
  mutate(sum(Freq))%>%
  mutate(percent=Freq/`sum(Freq)`)
head(bar_per)

tb2=table(phe$celltype, phe$group)
head(tb2)
balloonplot(tb2)
bar_data2<-as.data.frame(tb2)

bar_per2<-bar_data2%>%
  group_by(Var2)%>%
  mutate(sum(Freq))%>%
  mutate(percent=Freq/`sum(Freq)`)
head(bar_per2)
write.csv(bar_per2,file="epithelium_subcluster_by_group_percent.csv")

bar_per2$Var2 <- factor(bar_per2$Var2,levels=c("Control","CD","PB_IFX"))
bar_per2$Var1 <- factor(bar_per2$Var1,
                        levels=c('Stem','EECs_Neurod1','EECs_Lix1','Undiff_Entero',
                                 'Entero_Rbp2','Entero_Mep1a',
                                 'Entero_Sult1a1','Entero_Duoxa2_s100a9',
                                 'Goblet_Best2','Goblet_Cenpa',
                                 'Goblet_Gmds','Goblet_pre','Mucous_cells','Tuft',
                                 'TA'))

p1 <- ggplot(bar_per2,aes(x=Var1,y=percent,color=Var2, fill=Var2))+
  geom_bar(stat="identity",position = position_dodge(0.9),width=0.75,alpha=0.8)+
  scale_y_continuous(expand = c(0,0),
                     limits = c(-0.01,0.4))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size = 10,
                                   face="bold" ,angle=45, vjust=1,hjust=1),
        axis.text.y = element_text(colour = "black", size = 8),
        legend.position="top"
  )+
  labs(x="", y="Fraction of total cells")+
  scale_color_manual(values=c("#9B989B","#BE0F20","#272B82"))+
  scale_fill_manual(values=c("#9B989B","#BE0F20","#272B82"))
p1


#####
library(enrichplot)
library(ggstatsplot)
library(Vennerable) 
library(ggplot2)
library(pheatmap)
library(cowplot)
library(ggrepel)


Idents(sce) = sce$group
table(Idents(sce))
sce.1  = sce[, sce$celltype %in% c( 'Stem','EECs_Neurod1','EECs_Lix1',
                                    'Entero_Mep1a',
                                    'Entero_Sult1a1','Entero_Duoxa2_s100a9',
                                    'Goblet_Best2','Goblet_Cenpa',
                                    'Goblet_Gmds','Goblet_pre','Mucous_cells','Tuft',
                                    'TA')]

degs = lapply(unique(sce.1$celltype), function(x){
  FindMarkers(sce.1[,sce.1$celltype==x],ident.1 = 'PB_IFX',
              ident.2 = 'CD')
})
names(degs) <- unique(sce.1$celltype)

for (i in 1:length(degs)){
  dat <- degs[[i]]
  dat$cluster <- rep(names(degs)[[i]],times=nrow(dat))
  dat$gene <- rownames(dat)
  filename1= paste0(unique(dat$cluster),".clusterdeg.csv")
  write.csv(dat,file=filename1)
}

#####
cell <- c( 'Stem','EECs_Neurod1','EECs_Lix1',
           'Entero_Mep1a',
           'Entero_Duoxa2_s100a9',
           'Goblet_Best2','Goblet_Cenpa',
           'Goblet_Gmds','Goblet_pre','Mucous_cells','Tuft',
           'TA')
for (i in 1:12) {
  filename1= paste0(unique(cell[i]),".clusterdeg.csv")
  DEG = read.csv(filename1,row.names = 1)
  gene_up <- rownames(DEG[with(DEG,avg_log2FC>0.1 & p_val_adj<0.05),])
  gene_down <- rownames(DEG[with(DEG,avg_log2FC< -0.1 & p_val_adj<0.05),])
  gene_up_entrez <- as.character(na.omit(bitr(gene_up, 
                                              fromType="SYMBOL", 
                                              toType="ENTREZID", 
                                              OrgDb="org.Rn.eg.db")[,2])) #"org.Hs.eg.db" "org.Mm.eg.db"
  
  kegg_enrich_up_results <- enrichKEGG(gene  = gene_up_entrez,
                                       organism  = "rno" , 
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff = 0.2
  )
  kegg_enrich_up_results <- DOSE::setReadable(kegg_enrich_up_results, 
                                              OrgDb="org.Rn.eg.db", 
                                              keyType='ENTREZID')
  filename2= paste0(unique(cell[i]),".KEGG.up.list.csv")
  write.csv(kegg_enrich_up_results@result,filename2)
  
  gene_down_entrez <- as.character(na.omit(bitr(gene_down, 
                                                fromType="SYMBOL", 
                                                toType="ENTREZID", 
                                                OrgDb="org.Rn.eg.db")[,2])) #"org.Hs.eg.db" "org.Mm.eg.db"
  
  kegg_enrich_down_results <- enrichKEGG(gene  = gene_down_entrez,
                                         organism  = "rno" , 
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.2
  )
  kegg_enrich_down_results <- DOSE::setReadable(kegg_enrich_down_results, 
                                                OrgDb="org.Rn.eg.db", 
                                                keyType='ENTREZID')
  filename3= paste0(unique(cell[i]),".KEGG.down.list.csv")
  write.csv(kegg_enrich_down_results@result,filename3)
}


# 5. re-annotation -------
load('intermediate.cell.annotation.counts.Rdata')
load("immune/immune.annotated.Rdata")
imu <- imu
load("epithelium/epithelium.annotated.Rdata")
epi <- sce
load("stromal/stromal.annotated.Rdata")
str <- sce

Idents(imu)='new_celltype'
Idents(sce.all.1)='celltype'
Idents(sce.all.1, cells = colnames(imu)) <- Idents(imu)
sce.all.1$new_celltype <- Idents(sce.all.1)


Idents(epi)='celltype'
Idents(sce.all.1, cells = colnames(epi)) <- Idents(epi)
sce.all.1$new_celltype <- Idents(sce.all.1)


Idents(str)='celltype'
Idents(sce.all.1, cells = colnames(str)) <- Idents(str)
sce.all.1$new_celltype <- Idents(sce.all.1)
save(sce.all.1, file='all.cell.annotation.counts.new.Rdata')

DimPlot(sce.all.1, reduction='umap',
        pt.size = 0.1,
        group.by = "new_celltype",label = T,
        label.size=3,repel=T ) + NoLegend()
