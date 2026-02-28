#1. library ----
library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)
#2. load data ------
plan("multisession", workers = 4)
load("../all.cell.annotation.counts.Rdata")
table(sce.all.1$new_celltype)
celltypes_to_plot <- c(
  'EECs_Lix1',
  'Entero_Sult1a1','Entero_Duoxa2_s100a9','Goblet_Best2',
  'pDCs','cDC2_Zeb2','cDC1',
  'Macs_Mrc1_Cd163','Macs_Itgax','Macs_Spp1','Ly6c_Monocytes',
  'Neu_Csf1','Neu_Retnlg',
  'Cd8_Gzmm','Cd4_Cd8','Cd4_treg',
  'IgG_plasma_B','Memory_B','Plasmablast',
  'Endothelial','Lymphatics','Fibro_IL11_Cxcl14','Fibro_Smoc2_Ccl11',
  'Fibro_Mfap5_Gsn','Myofibroblasts','Cycling_stroma'
)

sce_CD <- sce.all.1[, sce.all.1$group %in% c( 'CD' )]
sce_CD <- sce_CD[, sce_CD$new_celltype %in% celltypes_to_plot]


data.input <- sce_CD[["RNA"]]$data
meta <- sce_CD@meta.data[,c("group","new_celltype")]
colnames(meta) <- c("group","celltypes")
table(meta$celltypes)
meta$celltypes = droplevels(meta$celltypes, 
                            exclude = setdiff(levels(meta$celltypes),unique(meta$celltypes)))

cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "celltypes")
levels(cellchat@idents)


CellChatDB <- CellChatDB.mouse
dplyr::glimpse(CellChatDB$interaction)
cellchat@DB <- CellChatDB
table(CellChatDB.mouse$interaction$annotation)

#3. cellchat2 analysis------
cellchat <- subsetData(cellchat)  
cellchat <- identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)  

cellchat <- computeCommunProb(cellchat,type="triMean") 

cellchat <- filterCommunication(cellchat, min.cells = 5)  

cellchat <- computeCommunProbPathway(cellchat)  
cellchat <- aggregateNet(cellchat)             
df.net <- subsetCommunication(cellchat)     

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(4,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

sce_PB <- sce.all.1[, sce.all.1$group %in% c( 'PB_IFX' )]
sce_PB <- sce_PB[, sce_PB$new_celltype %in% celltypes_to_plot]
data.input <- sce_PB[["RNA"]]$data
meta <- sce_PB@meta.data[,c("group","new_celltype")]
colnames(meta) <- c("group","celltypes")
table(meta$celltypes)
meta$celltypes = droplevels(meta$celltypes, 
                            exclude = setdiff(levels(meta$celltypes),unique(meta$celltypes)))

cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "celltypes")
levels(cellchat@idents)

CellChatDB <- CellChatDB.mouse
dplyr::glimpse(CellChatDB$interaction)
cellchat@DB <- CellChatDB
table(CellChatDB.mouse$interaction$annotation)

cellchat <- subsetData(cellchat)  
cellchat <- identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)  

cellchat <- computeCommunProb(cellchat,type="triMean") 
cellchat <- filterCommunication(cellchat, min.cells = 5)  

cellchat <- computeCommunProbPathway(cellchat)  
cellchat <- aggregateNet(cellchat)             
df.net <- subsetCommunication(cellchat)     


CDcellchat <-  readRDS("CD.cellchat.rds")
PBcellchat <- readRDS("PB.cellchat.rds")

object.list <- list(CD = CDcellchat,PB = PBcellchat )
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)



a<- df.net$CD  #查看CD组的pathway
b<- df.net$PB  #查看PB组的pathway

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
p

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2






library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(circlize)
library(grDevices)
library(cowplot)

setwd("/Users/wangying/work/project/Rat_CD_model/single_cell_sequencing/res_1/cell-cell/")

data = read.csv("group_prioritization.csv",row.names = 1)

df = data[,c(8,3,4,5,6,2,17)]

###CD
df_CD = df[df$group == "CD",]
df_CD = df_CD[1:50,]
df_CD = dplyr::mutate(df_CD, prioritization_rank = row_number())
group_oi = "CD"

prioritized_tbl_oi = df_CD %>% dplyr::ungroup() # if grouped: things will be messed up downstream

mycolors2 <- c(  
  "#fb8e41",    
  "#e31b1e",    
  "#ad5f2c",    
  "#f6e36d",   
  "#fad1e0",    
  "#36877a",    
  "#d3bad9",   
  "#b982bc",    
  "#caeac2",   
  "#bbe173",   
  "#65c2a4",   
  "#a8cee0",   
  "#89a6af",   
  "#1b79af",   
  "#543005"    
)

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), 
                          prioritized_tbl_oi$receiver %>% unique()) %>% sort()
colors_sender = mycolors2 %>% magrittr::set_names(senders_receivers)
colors_receiver = mycolors2 %>%  magrittr::set_names(senders_receivers)


prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() 
# Link each cell type to a color
grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
groups_oi = prioritized_tbl_oi$group %>% unique()

title = group_oi
circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
df = circos_links

ligand.uni = unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i = df[df$ligand == ligand.uni[i], ]
  sender.uni = unique(df.i$sender)
  for (j in 1:length(sender.uni)) {
    df.i.j = df.i[df.i$sender == sender.uni[j], ]
    df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
  }
}
receptor.uni = unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i = df[df$receptor == receptor.uni[i], ]
  receiver.uni = unique(df.i$receiver)
  for (j in 1:length(receiver.uni)) {
    df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
    df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
  }
}

intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))

while(length(intersecting_ligands_receptors) > 0){
  df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
  df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
  df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
  df = dplyr::bind_rows(df_unique, df_duplicated)
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
}

circos_links = df

# Link ligands/Receptors to the colors of senders/receivers
circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency

# Define order of the ligands and receptors and the gaps
ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
}) %>% unlist()

receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
}) %>% unlist()

order = c(ligand_order,receptor_order)

width_same_cell_same_ligand_type = 0.275
width_different_cell = 3
width_ligand_receptor = 9
width_same_cell_same_receptor_type = 0.275

#######
sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
sender_gaps = sender_gaps[-length(sender_gaps)]
# 
receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
receiver_gaps = receiver_gaps[-length(receiver_gaps)]
# 
gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
# 
# # print(length(gaps))
# # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
  warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
}
#########

links_circle$weight[links_circle$weight == 0] = 0.01
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle,
             directional = 1,
             order=order,
             link.sort = TRUE,
             link.decreasing = TRUE,
             grid.col = grid_col,
             # transparency = transparency,
             diffHeight = 0.0075,
             direction.type = c("diffHeight", "arrows"),
             link.visible = links_circle$weight > 0.01,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.175),
             grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
             reduce = 0,
             scale = TRUE)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
}, bg.border = NA) #
title(title)
p_circos = recordPlot()
#return(p_circos)


# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Receiver")
ComplexHeatmap::draw(legend, just = c("left", "bottom"))
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Sender")
ComplexHeatmap::draw(legend, just = c("left", "top"))

p_legend = grDevices::recordPlot()


#p_circos 
#p_legend

data = read.csv("group_prioritization.csv",row.names = 1)
cell <- c("Macs_Spp1","Macs_Mrc1_Cd163","Ly6c_Monocytes","Macs_Itgax")
signal <- c("Il1a_Il1rap","Il1b_Il1rap","Tnf_Tnfrsf1a","Tnf_Tnfrsf1b","Ccl3_Ccr5","Ccl2_Ccr2",
            "Cxcl1_Cxcr2","Cxcl2_Cxcr2","Icam1_Itgam","S100a9_Itgb2",
            "Anxa2_Tlr2","Hp_Tlr4","Bgn_Tlr4","Tgfb1_Tgfbr1",
            "Tgfb1_Itgav","Eng_Bmpr2","Spp1_Itgav","Thbs1_Itgb1",
            "Fn1_Itga5","Col1a1_Itga2","Mmp2_Sdc2","Mmp9_Lrp1",
            "Timp2_Mmp2","Timp1_Mmp9","Timp1_Mmp2",
            "Spp1_Itgav","Il1a_Il1r1","Il1b_Il1r1","Spp1_Itgb6","Spp1_Itga8",
            "Cd86_Ctla4","Cd86_Cd28","Sirpa_Cd69","Icosl_Cd28","Icosl_Ctla4",
            "Cxcl12_Cxcr4","Ccl7_Cxcr3","Cxcl14_Cxcr4","Osm_Osmr"
)
exclude_strs <- c("Mucous_cells", "Stem","IgG_plasma_B","Lymphatics","Mast",
                  "Fibro_Mfap5_Gsn","TA","Glial"
)
df_filter <- data[
  data$sender %in% cell & 
    grepl(paste(signal, collapse = "|"), data$lr_interaction) &
    !(data$receiver %in% exclude_strs), 
]

df_CD = df_filter[df_filter$group == "CD",]
unique(df_CD$receiver)
df = df_CD[,c(8,3,4,5,6,2,17)]

df_CD = dplyr::mutate(df_CD, prioritization_rank = row_number())
group_oi = "CD"

prioritized_tbl_oi = df_CD %>% dplyr::ungroup() # if grouped: things will be messed up downstream


mycolors2 <- c(  
  "#fb8e41",     
  "#FEE0B6",  
  "#F1B6DA",    
  "#417ab5",    
  "#65c2a4",    
  "#caeac2",    
  "#f6e36d",    
  "#ad5f2c",    
  "#a8cee0",    
  "#684A11",    
  "#762A83",   
  "#1B7837",    
  "#B32B34",    
  "#d3bad9",    
  "#7A5399",    
  "#5E1415"    
)

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), 
                          prioritized_tbl_oi$receiver %>% unique()) %>% sort()
colors_sender = mycolors2 %>% magrittr::set_names(senders_receivers)
colors_receiver = mycolors2 %>%  magrittr::set_names(senders_receivers)


prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() 
# Link each cell type to a color
grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
groups_oi = prioritized_tbl_oi$group %>% unique()

title = group_oi
circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
df = circos_links

ligand.uni = unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i = df[df$ligand == ligand.uni[i], ]
  sender.uni = unique(df.i$sender)
  for (j in 1:length(sender.uni)) {
    df.i.j = df.i[df.i$sender == sender.uni[j], ]
    df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
  }
}
receptor.uni = unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i = df[df$receptor == receptor.uni[i], ]
  receiver.uni = unique(df.i$receiver)
  for (j in 1:length(receiver.uni)) {
    df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
    df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
  }
}

intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))

while(length(intersecting_ligands_receptors) > 0){
  df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
  df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
  df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
  df = dplyr::bind_rows(df_unique, df_duplicated)
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
}

circos_links = df

# Link ligands/Receptors to the colors of senders/receivers
circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency

# Define order of the ligands and receptors and the gaps
ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
}) %>% unlist()

receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
}) %>% unlist()

order = c(ligand_order,receptor_order)

width_same_cell_same_ligand_type = 0.275
width_different_cell = 3
width_ligand_receptor = 9
width_same_cell_same_receptor_type = 0.275

#######
sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
sender_gaps = sender_gaps[-length(sender_gaps)]
# 
receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
receiver_gaps = receiver_gaps[-length(receiver_gaps)]
# 
gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
# 
# # print(length(gaps))
# # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
  warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
}
#########

links_circle$weight[links_circle$weight == 0] = 0.01
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle,
             directional = 1,
             order=order,
             link.sort = TRUE,
             link.decreasing = TRUE,
             grid.col = grid_col,
             # transparency = transparency,
             diffHeight = 0.0075,
             direction.type = c("diffHeight", "arrows"),
             link.visible = links_circle$weight > 0.01,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.175),
             grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
             reduce = 0,
             scale = TRUE)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
}, bg.border = NA) #
title(title)
p_circos = recordPlot()
#return(p_circos)


# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Receiver")
ComplexHeatmap::draw(legend, just = c("left", "bottom"))
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Sender")
ComplexHeatmap::draw(legend, just = c("left", "top"))

p_legend = grDevices::recordPlot()




data = read.csv("group_prioritization.csv",row.names = 1)
cell <- c("Macs_Spp1","Macs_Mrc1_Cd163","Ly6c_Monocytes","Macs_Itgax")
signal <- c("Il1a_Il1rap","Il1b_Il1rap","Tnf_Tnfrsf1a","Ccl3_Ccr5","Ccl2_Ccr2",
            "Cxcl1_Cxcr2","Cxcl2_Cxcr2","Icam1_Itgam","S100a9_Itgb2",
            "Anxa2_Tlr2","Hp_Tlr4","Bgn_Tlr4","Tgfb1_Tgfbr1",
            "Tgfb1_Itgav","Eng_Bmpr2","Spp1_Itgav","Thbs1_Itgb1",
            "Fn1_Itga5","Col1a1_Itga2","Mmp2_Sdc2","Mmp9_Lrp1","Timp2_Mmp2",
            "Spp1_Itgav","Il1a_Il1r1","Il1b_Il1r1","Spp1_Itgb6","Spp1_Itga8",
            "Cd86_Ctla4","Cd86_Cd28","Sirpa_Cd69","Icosl_Cd28","Icosl_Ctla4",
            "Cxcl12_Cxcr4","Ccl7_Cxcr3","Cxcl14_Cxcr4"
)
exclude_strs <- c("Mucous_cells", "Stem","IgG_plasma_B","Lymphatics","Mast",
                  "Fibro_Mfap5_Gsn","TA","Glial","Cd4_Cd8","Fibro_Smoc2_Ccl11"
)
df_filter <- data[
  data$receiver %in% cell & 
    grepl(paste(signal, collapse = "|"), data$lr_interaction) &
    !(data$sender %in% exclude_strs), 
]


df_CD = df_filter[df_filter$group == "CD",]
unique(df_CD$sender)
df = df_CD[,c(8,3,4,5,6,2,17)]

df_CD = dplyr::mutate(df_CD, prioritization_rank = row_number())
group_oi = "CD"

prioritized_tbl_oi = df_CD %>% dplyr::ungroup() # if grouped: things will be messed up downstream



mycolors2 <- c(  
  "#FEE0B6",  
  # "#fb8e41",    
  "#F1B6DA",    
  "#417ab5",    
  "#FE9929", 
  #"#65c2a4",    
  "#caeac2",    
  "#1b79af",    
  "#D79D9D", 
  "#ad5f2c",    
  "#003C30", 
  "#684A11",   
  "#762A83",    
  "#1B7837",   
  "#B32B34",   
  "#d3bad9",    
  "#7A5399",    
  "#A6DBA0", 
  "#E1DF97" 
)

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), 
                          prioritized_tbl_oi$receiver %>% unique()) %>% sort()
colors_sender = mycolors2 %>% magrittr::set_names(senders_receivers)
colors_receiver = mycolors2 %>%  magrittr::set_names(senders_receivers)


prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() 
# Link each cell type to a color
grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
groups_oi = prioritized_tbl_oi$group %>% unique()

title = group_oi
circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
df = circos_links

ligand.uni = unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i = df[df$ligand == ligand.uni[i], ]
  sender.uni = unique(df.i$sender)
  for (j in 1:length(sender.uni)) {
    df.i.j = df.i[df.i$sender == sender.uni[j], ]
    df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
  }
}
receptor.uni = unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i = df[df$receptor == receptor.uni[i], ]
  receiver.uni = unique(df.i$receiver)
  for (j in 1:length(receiver.uni)) {
    df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
    df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
  }
}

intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))

while(length(intersecting_ligands_receptors) > 0){
  df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
  df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
  df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
  df = dplyr::bind_rows(df_unique, df_duplicated)
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
}

circos_links = df

# Link ligands/Receptors to the colors of senders/receivers
circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency

# Define order of the ligands and receptors and the gaps
ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
}) %>% unlist()

receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
}) %>% unlist()

order = c(ligand_order,receptor_order)

width_same_cell_same_ligand_type = 0.275
width_different_cell = 3
width_ligand_receptor = 9
width_same_cell_same_receptor_type = 0.275

#######
sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
sender_gaps = sender_gaps[-length(sender_gaps)]
# 
receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
receiver_gaps = receiver_gaps[-length(receiver_gaps)]
# 
gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
# 
# # print(length(gaps))
# # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
  warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
}
#########

links_circle$weight[links_circle$weight == 0] = 0.01
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle,
             directional = 1,
             order=order,
             link.sort = TRUE,
             link.decreasing = TRUE,
             grid.col = grid_col,
             # transparency = transparency,
             diffHeight = 0.0075,
             direction.type = c("diffHeight", "arrows"),
             link.visible = links_circle$weight > 0.01,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.175),
             grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
             reduce = 0,
             scale = TRUE)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
}, bg.border = NA) #
title(title)
p_circos = recordPlot()
#return(p_circos)


# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Receiver")
ComplexHeatmap::draw(legend, just = c("left", "bottom"))
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Sender")
ComplexHeatmap::draw(legend, just = c("left", "top"))

p_legend = grDevices::recordPlot()



####3. Treg and T

cell <- c("Cd4_treg","Cd8_Gzmm")
signal <- c("Col4a1_Itga1",  "Col17a1_Itga1" ,"Col4a1_Itgb1",  "Col5a1_Itga1"  ,"Col17a1_Itgb1",
            "Col4a1_Itgav",  "Col5a1_Itgb1",  "Col1a2_Cd44" ,  "Col15a1_Itga1",
            "Col5a1_Sdc1" ,  "Col14a1_Cd44" , "Col1a1_Cd44"  , "Col6a2_Itga1" , "Col3a1_Itga1",
            "Col6a3_Itga1",  "Col12a1_Itga1" ,"Col5a3_Itgb1",  "Col1a1_Itga1" ,
            "Col15a1_Itgb1", "Col6a2_Itgb1",  "Col6a1_Itga1"  ,"Col5a2_Itga1",  
            "Col1a2_Itga1",  "Col14a1_Itga1", "Col1a2_Itgb1",  "Col4a1_Itgb8"  ,"Col3a1_Itgb1" ,
            "Prnp_Tnfrsf25" ,"Timp2_Mmp2","Omg_Tnfrsf1b","Timp2_Cd44" , "Timp2_Itgb1","Wnt5a_Lrp6",
            "Thbs1_Itga4", "Thbs4_Itgb1", "Thbs4_Cd47",  "Thbs2_Itga4", "Thbs4_Sdc1",  "Thbs1_Itgb1",
            "Mmp14_Cd44", "Mmp14_Sdc1", "Mmp14_Mmp2", "Mmp9_Itgb2", "Mmp9_Cd44" ,"Furin_Adam19",
            "C3_Itgb2","Nxph3_Nrxn1","App_Ncstn" ,"Lrrc4b_Ptprs" , "Calca_Ramp1","Plat_Itgb2",
            "Ptn_Ptprs","Fbn1_Itgb1","Agrn_Itgb1","Lum_Itgb1" ,"Lgals3_Lag3" ,"Hp_Itgb2" ,
            "Ccn2_Lrp6"  ,"Lum_Itgb1" ,"Lrfn1_Ptprs" ,
            "C3_Itgax","F13a1_Itga4" ,"Lyz2_Itgal" ,"Ncam1_Nptn","Hbegf_Cd6","Siglec1_Spn","Siglech_Kir3dl1" ,
            "Lamb3_Itgb1", "Lamb3_Itga6", "Lama3_Itgb1" ,"Lama3_Itga6", "Lama1_Nt5e",  "Lamc2_Itgb1"
)
exclude_strs <- c("Mucous_cells", "Stem","IgG_plasma_B","Lymphatics","Mast",
                  "TA","Glial","Cd4_Cd8","Fibro_Smoc2_Ccl11","Goblet_pre","Tuft",
                  "Memory_B",	"EECs_Neurod1"
)
df_filter <- data[
  data$receiver %in% cell &
    #grepl(paste(signal, collapse = "|"), data$lr_interaction) &
    !(data$lr_interaction %in% signal) &
    !(data$sender %in% exclude_strs), 
]

df_CD = df_filter[df_filter$group == "CD",]

df = df_CD[,c(8,3,4,5,6,2,17)]

df_CD = dplyr::mutate(df_CD, prioritization_rank = row_number())
group_oi = "CD"

prioritized_tbl_oi = df_CD %>% dplyr::ungroup() # if grouped: things will be messed up downstream

mycolors2 <- c(  
  "#fb8e41",     
  "#FEE0B6",  
  "#F1B6DA",   
  "#417ab5",    
  "#65c2a4",    
  "red", 
  "#caeac2",    
  "#1b79af",   
  "#D79D9D", 
  "#ad5f2c",    
  "#a8cee0",    
  "#684A11",    
  "#762A83",    
  "#1B7837",    
  "#B32B34",    
  "#d3bad9",    
  "#7A5399",    
  "#5E1415"     
)

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), 
                          prioritized_tbl_oi$receiver %>% unique()) %>% sort()
colors_sender = mycolors2 %>% magrittr::set_names(senders_receivers)
colors_receiver = mycolors2 %>%  magrittr::set_names(senders_receivers)


prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() 
# Link each cell type to a color
grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
groups_oi = prioritized_tbl_oi$group %>% unique()

title = group_oi
circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
df = circos_links

ligand.uni = unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i = df[df$ligand == ligand.uni[i], ]
  sender.uni = unique(df.i$sender)
  for (j in 1:length(sender.uni)) {
    df.i.j = df.i[df.i$sender == sender.uni[j], ]
    df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
  }
}
receptor.uni = unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i = df[df$receptor == receptor.uni[i], ]
  receiver.uni = unique(df.i$receiver)
  for (j in 1:length(receiver.uni)) {
    df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
    df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
  }
}

intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))

while(length(intersecting_ligands_receptors) > 0){
  df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
  df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
  df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
  df = dplyr::bind_rows(df_unique, df_duplicated)
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
}

circos_links = df

# Link ligands/Receptors to the colors of senders/receivers
circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency

# Define order of the ligands and receptors and the gaps
ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
}) %>% unlist()

receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
}) %>% unlist()

order = c(ligand_order,receptor_order)

width_same_cell_same_ligand_type = 0.275
width_different_cell = 3
width_ligand_receptor = 9
width_same_cell_same_receptor_type = 0.275

#######
sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
sender_gaps = sender_gaps[-length(sender_gaps)]
# 
receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
receiver_gaps = receiver_gaps[-length(receiver_gaps)]
# 
gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
# 
# # print(length(gaps))
# # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
  warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
}
#########

links_circle$weight[links_circle$weight == 0] = 0.01
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle,
             directional = 1,
             order=order,
             link.sort = TRUE,
             link.decreasing = TRUE,
             grid.col = grid_col,
             # transparency = transparency,
             diffHeight = 0.0075,
             direction.type = c("diffHeight", "arrows"),
             link.visible = links_circle$weight > 0.01,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.175),
             grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
             reduce = 0,
             scale = TRUE)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
}, bg.border = NA) #
title(title)
p_circos = recordPlot()
#return(p_circos)


# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Receiver")
ComplexHeatmap::draw(legend, just = c("left", "bottom"))
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Sender")
ComplexHeatmap::draw(legend, just = c("left", "top"))

p_legend = grDevices::recordPlot()


###
group_prioritization_tbl = multinichenet_output$prioritization_tables$group_prioritization_tbl
ligand <- c("Tgfb1","Igf1","Vegfa",   "Bdnf",
            "Il12a",   "Il1b",    "Il10",    "Igf1",
            "Bmp2",    "Bmp4",    "Il1a",    "Inhba",   "Tgfb3",   "Tgfb2",
            "Tgfb1",   "Osm",     "Cxcl12",  "Tnf",     "Il6",
            "Ifng",    "Bmp6","Tnfsf11", "Il21",    "Il15" )
test <- group_prioritization_tbl[
  group_prioritization_tbl$receiver %in% ("Fibro_IL11_Cxcl14") &
    group_prioritization_tbl$ligand %in% ligand &
    group_prioritization_tbl$prioritization_score > 0.8,
  #grepl(paste(signal, collapse = "|"), data$lr_interaction) &
]



#PBIFX
data = read.csv("group_prioritization.csv",row.names = 1)
cell <- c("Macs_Spp1","Macs_Mrc1_Cd163","Ly6c_Monocytes","Macs_Itgax")
signal <- c("Il1a_Il1rap","Il1b_Il1rap","Tnf_Tnfrsf1a","Ccl3_Ccr5","Ccl2_Ccr2",
            "Cxcl1_Cxcr2","Cxcl2_Cxcr2","Icam1_Itgam","S100a9_Itgb2",
            "Anxa2_Tlr2","Hp_Tlr4","Bgn_Tlr4","Tgfb1_Tgfbr1",
            "Tgfb1_Itgav","Eng_Bmpr2","Spp1_Itgav","Thbs1_Itgb1",
            "Fn1_Itga5","Col1a1_Itga2","Mmp2_Sdc2","Mmp9_Lrp1","Timp2_Mmp2",
            "Spp1_Itgav","Il1a_Il1r1","Il1b_Il1r1","Spp1_Itgb6","Spp1_Itga8",
            "Cd86_Ctla4","Cd86_Cd28","Sirpa_Cd69","Icosl_Cd28","Icosl_Ctla4",
            "Cxcl12_Cxcr4","Ccl7_Cxcr3","Cxcl14_Cxcr4","Osm_Osmr"
)
exclude_strs <- c("Mucous_cells", "Stem","IgG_plasma_B","Lymphatics","Mast",
                  "Fibro_Mfap5_Gsn","TA","Glial","Goblet_Best2"
)
df_filter <- data[
  data$sender %in% cell & 
    grepl(paste(signal, collapse = "|"), data$lr_interaction) &
    !(data$receiver %in% exclude_strs), 
]

df_PB = df_filter[df_filter$group == "PB_IFX",]
#df_PB = df_PB[1:50,]
#df_PB = data[data$group == "PB_IFX",]
unique(df_PB$receiver)
df = df_PB[,c(8,3,4,5,6,2,17)]

df_PB = dplyr::mutate(df_PB, prioritization_rank = row_number())
group_oi = "PB_IFX"

prioritized_tbl_oi = df_PB %>% dplyr::ungroup() # if grouped: things will be messed up downstream


mycolors2 <- c(  
  "#fb8e41",     
  "#FEE0B6",  
  "#417ab5",    
  "#D79D9D",   
  "#762A83",   
  "#1B7837",    
  "#ee6aa7",  
  "#A6DBA0"  
  )

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), 
                          prioritized_tbl_oi$receiver %>% unique()) %>% sort()
colors_sender = mycolors2 %>% magrittr::set_names(senders_receivers)
colors_receiver = mycolors2 %>%  magrittr::set_names(senders_receivers)


prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() 
# Link each cell type to a color
grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
groups_oi = prioritized_tbl_oi$group %>% unique()

title = group_oi
circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
df = circos_links

ligand.uni = unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i = df[df$ligand == ligand.uni[i], ]
  sender.uni = unique(df.i$sender)
  for (j in 1:length(sender.uni)) {
    df.i.j = df.i[df.i$sender == sender.uni[j], ]
    df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
  }
}
receptor.uni = unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i = df[df$receptor == receptor.uni[i], ]
  receiver.uni = unique(df.i$receiver)
  for (j in 1:length(receiver.uni)) {
    df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
    df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
  }
}

intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))

while(length(intersecting_ligands_receptors) > 0){
  df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
  df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
  df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
  df = dplyr::bind_rows(df_unique, df_duplicated)
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
}

circos_links = df

# Link ligands/Receptors to the colors of senders/receivers
circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency

# Define order of the ligands and receptors and the gaps
ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
}) %>% unlist()

receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
}) %>% unlist()

order = c(ligand_order,receptor_order)

width_same_cell_same_ligand_type = 0.275
width_different_cell = 3
width_ligand_receptor = 9
width_same_cell_same_receptor_type = 0.275

#######
sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
sender_gaps = sender_gaps[-length(sender_gaps)]
# 
receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
receiver_gaps = receiver_gaps[-length(receiver_gaps)]
# 
gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
# 
# # print(length(gaps))
# # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
  warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
}
#########

links_circle$weight[links_circle$weight == 0] = 0.01
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle,
             directional = 1,
             order=order,
             link.sort = TRUE,
             link.decreasing = TRUE,
             grid.col = grid_col,
             # transparency = transparency,
             diffHeight = 0.0075,
             direction.type = c("diffHeight", "arrows"),
             link.visible = links_circle$weight > 0.01,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.175),
             grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
             reduce = 0,
             scale = TRUE)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
}, bg.border = NA) #
title(title)
p_circos = recordPlot()
#return(p_circos)


# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Receiver")
ComplexHeatmap::draw(legend, just = c("left", "bottom"))
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Sender")
ComplexHeatmap::draw(legend, just = c("left", "top"))

p_legend = grDevices::recordPlot()

####PB@IFX
data = read.csv("group_prioritization.csv",row.names = 1)
cell <- c("Macs_Spp1","Macs_Mrc1_Cd163","Ly6c_Monocytes","Macs_Itgax")
signal <- c("Il1a_Il1rap","Il1b_Il1rap","Tnf_Tnfrsf1a","Ccl3_Ccr5","Ccl2_Ccr2",
            "Cxcl1_Cxcr2","Cxcl2_Cxcr2","Icam1_Itgam","S100a9_Itgb2",
            "Anxa2_Tlr2","Hp_Tlr4","Bgn_Tlr4","Tgfb1_Tgfbr1",
            "Tgfb1_Itgav","Eng_Bmpr2","Spp1_Itgav","Thbs1_Itgb1",
            "Fn1_Itga5","Col1a1_Itga2","Mmp2_Sdc2","Mmp9_Lrp1","Timp2_Mmp2",
            "Spp1_Itgav","Il1a_Il1r1","Il1b_Il1r1","Spp1_Itgb6","Spp1_Itga8",
            "Cd86_Ctla4","Cd86_Cd28","Sirpa_Cd69","Icosl_Cd28","Icosl_Ctla4",
            "Cxcl12_Cxcr4","Ccl7_Cxcr3","Cxcl14_Cxcr4"
)
exclude_strs <- c("Mucous_cells", "Stem","IgG_plasma_B","Lymphatics","Mast",
                  "Fibro_Mfap5_Gsn","TA","Glial","Goblet_Best2","Fibro_Smoc2_Ccl11"
)
df_filter <- data[
  data$receiver %in% cell & 
    grepl(paste(signal, collapse = "|"), data$lr_interaction) &
    !(data$sender %in% exclude_strs), 
]

df_PB = df_filter[df_filter$group == "PB_IFX",]
unique(df_PB$receiver)
df = df_PB[,c(8,3,4,5,6,2,17)]

df_PB = dplyr::mutate(df_PB, prioritization_rank = row_number())
group_oi = "PB_IFX"

prioritized_tbl_oi = df_PB %>% dplyr::ungroup() # if grouped: things will be messed up downstream


mycolors2 <- c(  
  #"#fb8e41",     
  "#FEE0B6",  
  "#F1B6DA",    
  "#417ab5",    
  "red", 
  "#ad5f2c",    
  #"#D79D9D",    
  "#762A83",   
  "#1B7837",    
  "#d3bad9",    
  "#A6DBA0"  
)

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), 
                          prioritized_tbl_oi$receiver %>% unique()) %>% sort()
colors_sender = mycolors2 %>% magrittr::set_names(senders_receivers)
colors_receiver = mycolors2 %>%  magrittr::set_names(senders_receivers)


prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() 
# Link each cell type to a color
grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
groups_oi = prioritized_tbl_oi$group %>% unique()

title = group_oi
circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
df = circos_links

ligand.uni = unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i = df[df$ligand == ligand.uni[i], ]
  sender.uni = unique(df.i$sender)
  for (j in 1:length(sender.uni)) {
    df.i.j = df.i[df.i$sender == sender.uni[j], ]
    df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
  }
}
receptor.uni = unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i = df[df$receptor == receptor.uni[i], ]
  receiver.uni = unique(df.i$receiver)
  for (j in 1:length(receiver.uni)) {
    df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
    df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
  }
}

intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))

while(length(intersecting_ligands_receptors) > 0){
  df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
  df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
  df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
  df = dplyr::bind_rows(df_unique, df_duplicated)
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
}

circos_links = df

# Link ligands/Receptors to the colors of senders/receivers
circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency

# Define order of the ligands and receptors and the gaps
ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
}) %>% unlist()

receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
}) %>% unlist()

order = c(ligand_order,receptor_order)

width_same_cell_same_ligand_type = 0.275
width_different_cell = 3
width_ligand_receptor = 9
width_same_cell_same_receptor_type = 0.275

#######
sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
sender_gaps = sender_gaps[-length(sender_gaps)]
# 
receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
receiver_gaps = receiver_gaps[-length(receiver_gaps)]
# 
gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
# 
# # print(length(gaps))
# # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
  warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
}
#########

links_circle$weight[links_circle$weight == 0] = 0.01
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle,
             directional = 1,
             order=order,
             link.sort = TRUE,
             link.decreasing = TRUE,
             grid.col = grid_col,
             # transparency = transparency,
             diffHeight = 0.0075,
             direction.type = c("diffHeight", "arrows"),
             link.visible = links_circle$weight > 0.01,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.175),
             grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
             reduce = 0,
             scale = TRUE)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
}, bg.border = NA) #
title(title)
p_circos = recordPlot()
#return(p_circos)


# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Receiver")
ComplexHeatmap::draw(legend, just = c("left", "bottom"))
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Sender")
ComplexHeatmap::draw(legend, just = c("left", "top"))

p_legend = grDevices::recordPlot()


#####PB@IFX T
data = read.csv("group_prioritization.csv",row.names = 1)
cell <- c("Cd4_treg","Cd8_Gzmm")
signal <- c("Col4a1_Itga1",  "Col17a1_Itga1" ,"Col4a1_Itgb1",  "Col5a1_Itga1"  ,"Col17a1_Itgb1",
            "Col4a1_Itgav",  "Col5a1_Itgb1",  "Col1a2_Cd44" ,  "Col15a1_Itga1",
            "Col5a1_Sdc1" ,  "Col14a1_Cd44" , "Col1a1_Cd44"  , "Col6a2_Itga1" , "Col3a1_Itga1",
            "Col6a3_Itga1",  "Col12a1_Itga1" ,"Col5a3_Itgb1",  "Col1a1_Itga1" ,
            "Col15a1_Itgb1", "Col6a2_Itgb1",  "Col6a1_Itga1"  ,"Col5a2_Itga1",  
            "Col1a2_Itga1",  "Col14a1_Itga1", "Col1a2_Itgb1",  "Col4a1_Itgb8"  ,"Col3a1_Itgb1" ,
            "Prnp_Tnfrsf25" ,"Timp2_Mmp2","Omg_Tnfrsf1b","Timp2_Cd44" , "Timp2_Itgb1","Wnt5a_Lrp6",
            "Thbs1_Itga4", "Thbs4_Itgb1", "Thbs4_Cd47",  "Thbs2_Itga4", "Thbs4_Sdc1",  "Thbs1_Itgb1",
            "Mmp14_Cd44", "Mmp14_Sdc1", "Mmp14_Mmp2", "Mmp9_Itgb2", "Mmp9_Cd44" ,"Furin_Adam19",
            "C3_Itgb2","Nxph3_Nrxn1","App_Ncstn" ,"Lrrc4b_Ptprs" , "Calca_Ramp1","Plat_Itgb2",
            "Ptn_Ptprs","Fbn1_Itgb1","Agrn_Itgb1","Lum_Itgb1" ,"Lgals3_Lag3" ,"Hp_Itgb2" ,
            "Ccn2_Lrp6"  ,"Lum_Itgb1" ,"Lrfn1_Ptprs" ,"Hbegf_Cd9","Col4a5_Cd93",
            "Col4a6_Cd93",
            "C3_Itgax","F13a1_Itga4" ,"Lyz2_Itgal" ,"Ncam1_Nptn","Hbegf_Cd6","Siglec1_Spn","Siglech_Kir3dl1" ,
            "Lamb3_Itgb1", "Lamb3_Itga6", "Lama3_Itgb1" ,"Lama3_Itga6", "Lama1_Nt5e",  "Lamc2_Itgb1"
)
exclude_strs <- c("Mucous_cells", "Stem","IgG_plasma_B","Lymphatics","Mast",
                  "TA","Glial","Cd4_Cd8","Fibro_Smoc2_Ccl11","Goblet_pre","Tuft",
                  "Memory_B",	"EECs_Neurod1","Goblet_Best2","Goblet_Cenpa",
                  "Plasmablast"
)
df_filter <- data[
  data$receiver %in% cell &
    #grepl(paste(signal, collapse = "|"), data$lr_interaction) &
    !(data$lr_interaction %in% signal) &
    !(data$sender %in% exclude_strs), 
]

df_PB = df_filter[df_filter$group == "PB_IFX",]
unique(df_PB$receiver)
df = df_PB[,c(8,3,4,5,6,2,17)]

df_PB = dplyr::mutate(df_PB, prioritization_rank = row_number())
group_oi = "PB_IFX"

prioritized_tbl_oi = df_PB %>% dplyr::ungroup() # if grouped: things will be messed up downstream


mycolors2 <- c(  
  "#fb8e41",     
  "#FEE0B6",  
  "#F1B6DA",    
  "#417ab5",    
  "#65c2a4",   
  "red", 
  "#caeac2",    
  "#1b79af",    
  "#D79D9D", 
  "#ad5f2c",   
   "#a8cee0",    
  #"#684A11",    
  "#762A83",    
  "#1B7837",   
  "#d3bad9"   
)

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), 
                          prioritized_tbl_oi$receiver %>% unique()) %>% sort()
colors_sender = mycolors2 %>% magrittr::set_names(senders_receivers)
colors_receiver = mycolors2 %>%  magrittr::set_names(senders_receivers)


prioritized_tbl_oi = prioritized_tbl_oi %>% dplyr::ungroup() 
# Link each cell type to a color
grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
groups_oi = prioritized_tbl_oi$group %>% unique()

title = group_oi
circos_links_oi = prioritized_tbl_oi %>% dplyr::filter(group == group_oi)
circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)
df = circos_links

ligand.uni = unique(df$ligand)
for (i in 1:length(ligand.uni)) {
  df.i = df[df$ligand == ligand.uni[i], ]
  sender.uni = unique(df.i$sender)
  for (j in 1:length(sender.uni)) {
    df.i.j = df.i[df.i$sender == sender.uni[j], ]
    df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
    df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
  }
}
receptor.uni = unique(df$receptor)
for (i in 1:length(receptor.uni)) {
  df.i = df[df$receptor == receptor.uni[i], ]
  receiver.uni = unique(df.i$receiver)
  for (j in 1:length(receiver.uni)) {
    df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
    df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
    df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
  }
}

intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))

while(length(intersecting_ligands_receptors) > 0){
  df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
  df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
  df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
  df = dplyr::bind_rows(df_unique, df_duplicated)
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
}

circos_links = df

# Link ligands/Receptors to the colors of senders/receivers
circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency

# Define order of the ligands and receptors and the gaps
ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
}) %>% unlist()

receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
}) %>% unlist()

order = c(ligand_order,receptor_order)

width_same_cell_same_ligand_type = 0.275
width_different_cell = 3
width_ligand_receptor = 9
width_same_cell_same_receptor_type = 0.275

#######
sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
  sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
sender_gaps = sender_gaps[-length(sender_gaps)]
# 
receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
  sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
  gap = width_different_cell
  return(c(sector,gap))
}) %>% unlist()
receiver_gaps = receiver_gaps[-length(receiver_gaps)]
# 
gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
# 
# # print(length(gaps))
# # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
  warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
}
#########

links_circle$weight[links_circle$weight == 0] = 0.01
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle,
             directional = 1,
             order=order,
             link.sort = TRUE,
             link.decreasing = TRUE,
             grid.col = grid_col,
             # transparency = transparency,
             diffHeight = 0.0075,
             direction.type = c("diffHeight", "arrows"),
             link.visible = links_circle$weight > 0.01,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.175),
             grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
             reduce = 0,
             scale = TRUE)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
}, bg.border = NA) #
title(title)
p_circos = recordPlot()
#return(p_circos)


# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Receiver")
ComplexHeatmap::draw(legend, just = c("left", "bottom"))
legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                type = "grid",
                                legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                title_position = "topleft",
                                title = "Sender")
ComplexHeatmap::draw(legend, just = c("left", "top"))

p_legend = grDevices::recordPlot()




