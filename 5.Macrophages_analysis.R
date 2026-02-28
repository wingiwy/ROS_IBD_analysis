# 1. library -----
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(presto)
#if(!require(dplyr))install.packages("dplyr")
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
sce_macs <-  imu[, imu$new_celltype %in% c( 'Macs_Itgax',"Macs_Mrc1_Cd163","Ly6c_Monocytes","Macs_Spp1" )]
#3. clusterprofile -----
sce.markers<-FindAllMarkers(object=sce_macs,only.pos=TRUE,
                            min.pct=0.25,
                            thresh.use=0.25)

markers<-sce.markers|>group_by(cluster)|>
  filter(p_val_adj<0.005)|>
  ungroup()

library(clusterProfiler)
gid<-bitr(unique(markers$gene),'SYMBOL','ENTREZID',OrgDb='org.Rn.eg.db')
markers<-full_join(markers,gid,by=c('gene'='SYMBOL'))
#Itgax
m <- markers[grep("Macs_Itgax" , markers$cluster),]
gene <- m$gene
gene_entrez <- as.character(na.omit(bitr(gene, 
                                         fromType="SYMBOL", 
                                         toType="ENTREZID", # 
                                         OrgDb="org.Rn.eg.db")[,2])) #"org.Hs.eg.db" "org.Mm.eg.db"

#GO
ego <- enrichGO(
  gene          = m$ENTREZID,
  OrgDb         = "org.Rn.eg.db",
  keyType       = "ENTREZID",
  ont           = "ALL",         
  pAdjustMethod = "BH",         
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE           
)
go_results <- ego@result %>%
  filter(p.adjust < 0.05) %>%    
  arrange(p.adjust)             

write.csv(go_results, "macs_function/Itgax.go.csv", row.names = FALSE)

#Mrc1
m <- markers[grep("Macs_Mrc1_Cd163" , markers$cluster),]
gene <- m$gene
gene_entrez <- as.character(na.omit(bitr(gene, 
                                         fromType="SYMBOL", 
                                         toType="ENTREZID", 
                                         OrgDb="org.Rn.eg.db")[,2])) #"org.Hs.eg.db" "org.Mm.eg.db"

#GO分析
ego <- enrichGO(
  gene          = m$ENTREZID,
  OrgDb         = "org.Rn.eg.db",
  keyType       = "ENTREZID",
  ont           = "ALL",          
  pAdjustMethod = "BH",          
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE           
)
go_results <- ego@result %>%
  filter(p.adjust < 0.05) %>%    
  arrange(p.adjust)              

write.csv(go_results, "macs_function/Mrc1.go.csv", row.names = FALSE)

kegg <- read.csv("macs_function/pathway.csv")
kegg$'-log10pvalue' <- -log10(kegg$pvalue)

level_up <- kegg$Description[kegg$Cluster=='Itgax']
level_down <- kegg$Description[kegg$Cluster=='Mrc1']
level <- c(level_up)
level <- c(level_down)

kegg_Up <- kegg[kegg$Cluster=='Itgax']
kegg$Description <- factor(kegg$Description,
                           levels = rev(level))

mytheme <- theme(
  legend.position = 'none',
  #axis.text.y = element_blank(),
  #axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'grey60',size = 1.1),
  axis.text = element_text(size = 12)
)

up <- kegg[which(kegg$Cluster=='Itgax'),]
down <- kegg[which(kegg$Cluster=='Mrc1'),]
mycol <- c("#272B82","#BE0F20")
p1 <- ggplot(kegg,aes(x = `-log10pvalue`,y = Description,fill = Cluster)) + 
  geom_col(alpha=0.8) +
  theme_classic() + mytheme +
  geom_text(data = up,
            aes(x = -0.2, y = Description, label = Description),
            size = 3.5,
            hjust = 1)+
  geom_text(data = down,
            aes(x = 0.2, y = Description, label = Description),
            size = 3.5,
            hjust = 0) +
  scale_fill_manual(values = mycol)
p1


#4. pseudotime analysis-------
library(Seurat)
library(monocle)
#library(SeuratWrappers)
library(dplyr)

scRNA = subset(sce.all.1,idents = c("Macs_Mrc1_Cd163","Macs_Itgax","Macs_Spp1","Macs_Nos2"))
table(scRNA$new_celltype)
t <- scRNA@assays$RNA$counts
gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  row.names = row.names(ct)
)

pd <- new("AnnotatedDataFrame",
          data=scRNA@meta.data)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)

sc_cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds

sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- estimateDispersions(sc_cds)


fdif = "diff_test_res.Rdata"
if(!file.exists(fdif)){
  diff_test_res <- differentialGeneTest(sc_cds,
                                        fullModelFormulaStr = " ~ new_celltype", 
                                        #reducedModelFormulaStr = " ~ orig.ident", 
                                        relative_expr=TRUE, cores=4)
  save(diff_test_res,file = fdif)
}
load(fdif)


ordering_genes <- VariableFeatures(scRNA) 

head(ordering_genes)
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)
plot_ordering_genes(sc_cds)

sc_cds <- reduceDimension(sc_cds,residualModelFormulaStr = "~orig.ident")

sc_cds <- reduceDimension(sc_cds) 

sc_cds <- orderCells(sc_cds) 

#绘图
library(ggsci)
p1 = plot_cell_trajectory(sc_cds)+ scale_color_nejm()
p2 = plot_cell_trajectory(sc_cds, color_by = 'Pseudotime') 
p3 = plot_cell_trajectory(sc_cds, color_by = 'new_celltype')  + scale_color_npg()
library(patchwork)
p2+p1/p3


cols_1 <- c("#1B7837","#762A83","#B32B34","#684A11")
plot_cell_trajectory(sc_cds, color_by = 'new_celltype') + 
  scale_color_manual(values = cols_1)

plot_cell_trajectory(sc_cds, color_by = "Pseudotime")+ scale_color_gradientn(
  colors = c("#2E3192","#dddddc", "#FF0000"))


gene_to_cluster = diff_test_res %>% arrange(qval) %>% head(50) %>% pull(gene_short_name);head(gene_to_cluster)
gene_to_cluster = diff_test_res %>% arrange(qval) %>% head(100) %>% pull(gene_short_name);head(gene_to_cluster)


plot_pseudotime_heatmap(sc_cds[gene_to_cluster,],
                        num_clusters = nlevels(Idents(scRNA)), 
                        #num_clusters = 4,
                        show_rownames = TRUE,
                        cores = 1,return_heatmap = TRUE,
                        hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))


plot_genes_in_pseudotime(
  sc_cds[c("Spp1","Mrc1",'Cd163','Nos2','Mmp12','Ly6c',"Gpnmb","Folr2"), ],
  color_by = "new_celltype",
  ncol = 2
) +scale_color_manual(values = cols_1)

#5. macrophage polarization analysis -------
sce_macs <-  sce[, sce$celltype %in% c( 'Macs_Itgax',"Macs_Mrc1_Cd163","Macs_Nos2","Macs_Spp1" )]
library(tidyverse)
sce_macs <- sce_macs %>% NormalizeData %>% FindVariableFeatures %>% ScaleData

gl = list(
  M1 =  c( 'Tnf', 'Il1b', 'Il6','Nos2',"Ly6c","Spp1","Trem2"),
  M2 =  c( 'Il10', 'Cd163', 'Mrc1',"Arg1","Apoe","Folr2",'Selenop','Dap','Stab1')
)

sce_macs =  AddModuleScore(object = sce_macs,features = gl )
colnames(sce_macs@meta.data)
plot(sce_macs$Cluster1,sce_macs$Cluster2 ,col= sce_macs$celltype)

cols_1 <- c("#1B7837","#762A83","#B32B34","#684A11")

p <- ggplot(sce_macs@meta.data, aes(x=sce_macs$Cluster1,y=sce_macs$Cluster2,color=sce_macs$celltype))+
  geom_point(size=0.7)+
  scale_color_manual(values=c("#762A83","#1B7837","#684A11","#B32B34"))+
  scale_x_continuous(breaks = seq(-0.8,1.6,0.4),
                     limits = c(-0.8,1.6))+
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 1) +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 1) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 12),   
        axis.text.y = element_text(colour = 'black',size = 12),     
        axis.title = element_text(colour = 'black',size = 14)     
  )+
  labs(x="M1 Polarization",y="M2 Polarization")
p

#6. DEG gene KEGG analysis-----
kegg <- read.csv("DEGS/pseudobulk/CD_vs_PB_IFX_Macs_Itgax/keggpathway-down.csv")
kegg$'-log10pvalue' <-  -log10(kegg$pvalue)
kegg$Description <- factor(kegg$Description, levels = rev(kegg$Description))
p1 <- ggplot(kegg,aes(x = `-log10pvalue`,y = Description,fill = "#762A83")) + 
  geom_col(alpha=0.8) +
  theme_classic()
p1

kegg <- read.csv("DEGS/pseudobulk/CD_vs_PB_IFX_Macs_Mrc1_Cd163/pathway-down.csv")
kegg$'-log10pvalue' <-  -log10(kegg$pvalue)
kegg$Description <- factor(kegg$Description, levels = rev(kegg$Description))
p2 <- ggplot(kegg,aes(x = `-log10pvalue`,y = Description,fill = "#1B7837")) + 
  geom_col(alpha=0.8) +
  theme_classic()
p2


#7. heatmap -------
genes <- c(
  'Fcgr1a','Fcgr2a-ps1', 'Fcgr2b','Fcgr3a', 
  'Ifngr1','Ifngr2',                
  'Adgre1',       
  'Rac1','Cdc42',     
  'Ly6c','Ccr2','Plac8', 'Il6',      
  'Lipa','Treml4',       
  'Il10',       
  'Arg1',   
  'Mmp9','Vegfa','Vegfb','Pdgfa', 'Pdgfc',       
  'Folr2',
  'Timd4',"Mertk",'Pparg' 
)
gene.exp <- AverageExpression(sce,features = genes,
                              group.by = 'celltype',
                              slot = 'data')
gene.exp <- as.data.frame(gene.exp$RNA)

exp <- t(scale(t(gene.exp),scale=T,center = T))
exp <- exp[,c(6,7,4,5)]

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)
Breaks = c(seq(min(exp), 0, 
               length.out=ceiling(paletteLength/2) + 1),
           seq(max(exp)/paletteLength, 
               max(exp), 
               length.out=floor(paletteLength/2)))
pheatmap(exp,
         fontsize=14, 
         cluster_cols  = F,
         cluster_rows = F,
         fontsize_row = 10, 
         color=myColor, breaks = Breaks, 
         #angle_col = 90,
         treeheight_col = 0,  border_color = "white")


#8. gene expression heatmap -----
matu.genes <- c(
  'Ly6c','Ccr2','Sell','Trem1',"Itgam", "Il6",'Csf3',              
  'RT1-DMa','RT1-DMb',"RT1-Bb",
  "RT1-Ba",'RT1-Db1','Stab1','Ccl2',
  'C1qa','C1qb','C1qc', 'Adgre1',"Fcgr1a","Cx3cr1",
  'Mrc1','Cd163','Retnla','Il10','Tgfb1',   
  'Mmp12','Tlr4','Mmp9','Tnf','Il1b','Il1a','Nos2','Osm',
  'Spp1','Sod2','Trem2','Gpnmb','Hbegf'             
)
NETs.gene.exp <- AverageExpression(sce_macs,features = matu.genes,
                                   group.by = 'celltype',
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

pheatmap(NETs_exp,
         fontsize=14, 
         cluster_cols  = F,
         cluster_rows = F,
         fontsize_row = 10, 
         color=myColor, breaks = Breaks, 
         #angle_col = 90,
         treeheight_col = 0,  border_color = "white")

#8.KEGG pathway barplot -----
kegg <- read.csv("macs_function/pathway.csv")

dt <- kegg[1:5,]
dt$'-log10pvalue' <-  -log10(dt$pvalue)
dt$Description <- factor(dt$Description,
                         levels = rev(dt$Description))

mycol <- c("#272B82")
p1 <- ggplot(dt,aes(x = `-log10pvalue`,y = Description,fill = Cluster)) + 
  geom_col(alpha=0.8) +
  theme_classic() +
  scale_fill_manual(values = mycol)
p1

dt2 <- kegg[6:10,]
dt2$'-log10pvalue' <-  -log10(dt2$pvalue)
dt2$Description <- factor(dt2$Description,
                          levels = rev(dt2$Description))
mycol <- c("#272B82")
p2 <- ggplot(dt2,aes(x = `-log10pvalue`,y = Description,fill = Cluster)) + 
  geom_col(alpha=0.8) +
  theme_classic() +
  scale_fill_manual(values = mycol)
p2

p1+p2



