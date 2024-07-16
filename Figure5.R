#--------------------------------------------------------------
# filename : Figure5.R
# Date : 2023-09-01
# contributor : Minghao Dang, PhD
# function: 
# R version: R/4.2.0
#--------------------------------------------------------------

source('global_settings.R')

###################### Fig 5A MGUS06s(SMM-01) B cell trajectory ######################
monocle_cds<-readRDS('SMM-01_Bcell.monocle3.rds')

p1<-plot_cells(monocle_cds, color_cells_by = 'cell.type',
               label_groups_by_cluster=F,
               label_cell_groups=T,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=5,
               group_label_size=5,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 1,alpha=1)+ theme(legend.position = "none")

p2<-plot_cells(monocle_cds, color_cells_by = 'pseudotime',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 1,alpha=1)+
  theme(legend.position = "none")

p3<-plot_cells(monocle_cds, color_cells_by = 'infercnv.score.sd.normal_CD138_ref',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 1,alpha=1)+
  scale_color_gradientn(colors=rev(c("#67001f","#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","white","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061","#053061")),na.value = "#f0f0f0")+
  theme(legend.position = "none")

p4<-plot_cells(monocle_cds, color_cells_by = 'BCR_clonotype_raw',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 1,alpha=1)+ theme(legend.position = "none")

p5<-plot_cells(monocle_cds, color_cells_by = 'BCR',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 1,alpha=1)+ theme(legend.position = "none")

aa<-read.table('./CytoTRACE/CytoTRACE_plot_table.txt',sep='\t',header = T,row.names = 1)
monocle_cds$CytoTRACE<-aa[colnames(monocle_cds),]$CytoTRACE

temp <- RColorBrewer::brewer.pal(11, "Spectral")
temp[6] <- "gold"
rbPal <- colorRampPalette(temp)(length(CytoTRACE_value))

p6<-plot_cells(monocle_cds, color_cells_by = 'CytoTRACE',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_roots = F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 1,alpha=1)+
  scale_color_gradientn(colors=rev(temp),na.value = "#f0f0f0")+
  theme(legend.position = "none")

do.call(ggarrange,c(list(p1,p2,p3,p4,p5,p6),ncol = 1,nrow = 1)) -> combined.gp
pdf(paste0('SMM-01_Bcell_monocle3.',Sys.Date(),'.pdf'),height=5,width = 5.5)
print(combined.gp)
dev.off()

###################### Fig 5B MGUS06s(SMM-01) CD138+ inferCNV visualization ######################
library('data.table')
fdat = fread('./SMM-01/infercnv.observations.txt', sep = ' ')
dat_obs = data.frame(t(fdat[,-1]),check.names = F)
colnames(dat_obs)=fdat$V1
dat_obs = data.frame('barcode'=rownames(dat_obs),dat_obs,check.names = F)

gene_order<-readRDS('SMM-01.infercnv.gene_order.rds')

meta.data<-read.table('SMM-01.CD138_pos.filter.meta.data.txt',sep='\t',header=T)

meta.data<-meta.data[,c('barcode','cluster')]

dat_obs<-dat_obs[meta.data$barcode,]
dat_obs<-left_join(meta.data[,c('barcode','cluster')],dat_obs,by=c('barcode'='barcode'))

avg.cnv.df<-dat_obs%>%
  select(-barcode)%>%
  group_by(cluster)%>%
  summarize_all(mean, na.rm=TRUE)%>%
  as.data.frame(check.names=F)


cluster.cnv.score<-data.frame('cluster'=avg.cnv.df$cluster,
                              'cnv.score'=apply(avg.cnv.df[,-1],1,sd))
rownames(cluster.cnv.score)<-cluster.cnv.score$cluster

rownames(avg.cnv.df)<-avg.cnv.df$cluster
avg.cnv.df<-avg.cnv.df[,-1]

anno.col<-gene_order[,'chr',drop=FALSE]
anno.col$chr<-factor(anno.col$chr,levels=paste0('chr',1:22))

mat_colors <- list(chr = c(rep(c("lightgrey", "darkgrey"), 11)))
names(mat_colors$chr) <- unique(anno.col$chr)

cols_gain <- c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","white")
colfunc_gain <- colorRampPalette(cols_gain)

cols_loss <- c("white","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061")
colfunc_loss <- colorRampPalette(cols_loss)

cnv.cols <- c(colfunc_gain(20),'white',colfunc_loss(round(20*(1-min(avg.cnv.df))/(max(avg.cnv.df)-1),0)))


library('ComplexHeatmap')
top_ann = HeatmapAnnotation( df=anno.col,
                             which='column',
                             annotation_name_side = "left",
                             col=mat_colors,
                             show_legend = F
)


avg.cnv.df<-avg.cnv.df[c('SMM-01_C5','SMM-01_C2','SMM-01_C1','SMM-01_C3','SMM-01_C8','SMM-01_C6','SMM-01_C7','SMM-01_C0','SMM-01_C4'),]
pdf(paste0('SMM01.inferCNV.cluster.level.heatmap.',Sys.Date(),'.pdf'),height=2.5,width=15)
ht_list = Heatmap(avg.cnv.df, name = "inferCNV",
                  col=rev(cnv.cols),
                  cluster_rows = F,
                  cluster_columns = F,
                  show_row_names=T,
                  row_names_side = "left",
                  show_column_names = F,
                  #column_names_side = "bottom",
                  show_heatmap_legend = T,
                  
                  top_annotation = top_ann,
                  #left_annotation = left_ann,
                  #right_annotation = rowAnnotation(cnv.score = anno_points(cluster.cnv.score[rownames(avg.cnv.df),'cnv.score'])),
                  row_split = factor(rownames(avg.cnv.df),levels=rownames(avg.cnv.df)),
                  row_title = NULL,
                  column_split = anno.col$chr,
                  column_title = NULL,
                  row_gap = unit(0, "mm"), 
                  column_gap = unit(0, "mm"), 
                  border = TRUE,
                  
                  heatmap_legend_param = list(
                    #at = c(-2, 0, 2),
                    #labels = c("low", "zero", "high"),
                    #title = "Some values",
                    title_gp = gpar(fontsize = 12, fontface = "bold"),
                    labels_gp = gpar(fontsize = 12),
                    legend_height = unit(5, "cm"),
                    title_position ='topleft'
                  )) 
draw(ht_list)#, annotation_legend_list = lgd_list, annotation_legend_side = "top")
dev.off()

###################### Fig 5B MGUS06s(SMM-01) CD138+ inferCNV phylogenetic tree ######################
library('data.table')
library('ggtree')
library('tibble')
library('tidyr')
library('ape')
library('dplyr')
library('ggplot2')

fdat = fread('./SMM-01/infercnv.observations.txt', sep = ' ')

dat_obs = data.frame(t(fdat[,-1]),check.names = F)
colnames(dat_obs)=fdat$V1

meta.data<-read.table('SMM-01.CD138_pos.filter.meta.data.txt',sep='\t',header=T)
meta.data<-meta.data[,c('barcode','cluster')]

dat_obs = data.frame('barcode'=rownames(dat_obs),dat_obs,check.names = F)
dat_obs<-dat_obs[meta.data$barcode,]
dat_obs<-left_join(meta.data[,c('barcode','cluster')],dat_obs,by=c('barcode'='barcode'))

avg.cnv.df<-dat_obs%>%
  select(-barcode)%>%
  group_by(cluster)%>%
  summarize_all(mean, na.rm=TRUE)%>%
  as.data.frame(check.names=F)

rownames(avg.cnv.df)<-avg.cnv.df$cluster
avg.cnv.df<-avg.cnv.df[,-1]


run_me_tree <- function(consensus_df,
                        clusters,
                        ploidy_VAL,
                        rotate_nodes = NULL,
                        plot = TRUE) {
  
  consensus_int <- ploidy_VAL*consensus_df
  
  #adding a neutral state, will use as root
  consensus_int[nrow(consensus_int) + 1, ] <- round(ploidy_VAL)
  consensus_int[nrow(consensus_int) + 1, ] <- round(ploidy_VAL)
  
  tree <- ape::fastme.bal(dist(consensus_int, method = "manhattan"))
  
  tree <-
    root.phylo(tree,
               outgroup = which(tree$tip.label == Ntip(tree)),
               resolve.root = T)
  
  tree <-
    drop.tip(tree, tip = as.character(c(
      nrow(consensus_int), nrow(consensus_int) - 1
    )))
  
  tree <- ladderize(tree)
  
  # getting order
  is_tip <- tree$edge[, 2] <= length(tree$tip.label)
  ordered_tips_index <- tree$edge[is_tip, 2]
  tree_tips_order <- tree$tip.label[ordered_tips_index] %>% rev()
  
  # adding superclones information to frequencies_df
  clones_df <- clusters %>%
    distinct(cluster) %>%
    mutate(taxa = cluster)
  
  clones_df <- clones_df[c("taxa","cluster")]
  
  if (!is.null(rotate_nodes)) {
    tree <-
      phytools::rotateNodes(tree, rotate_nodes)
    is_tip <-
      tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips_index <- tree$edge[is_tip, 2]
    tree_tips_order <-
      tree$tip.label[ordered_tips_index] %>% rev()
    tree_tips_order <-
      tree_tips_order[str_detect(tree_tips_order, "c")]
  }
  
  p <- ggtree(tree, ladderize = F, size = 2) +
    geom_treescale()
  
  p <- p %<+% clones_df +
    geom_tippoint(aes(color = cluster),
                  size = 10) 
  
  if(plot) print(p)
  
  results <- list(tree = tree,
                  cs_plot = p,
                  cs_tree_order = tree_tips_order)
  
  return(results)
  
}

consensus_df<-avg.cnv.df;
clusters<-meta.data;
ploidy_VAL<-2;

pdf('SMM-01.phylo_tree.pdf')
run_me_tree(avg.cnv.df,meta.data,2)
dev.off()

###################### Fig 5B MGUS06s(SMM-01) B cell DEG heatmap ######################
SMM01_Bcell<-readRDS('SMM-01.Bcell.filter.scs.rds')

Idents(SMM01_Bcell)<-factor(SMM01_Bcell$cell.type,levels=c('Pre B cells','Proliferating B cells','B cells','SMM-01_C5','SMM-01_C2','SMM-01_C1','SMM-01_C3',
                                                           'SMM-01_C8','SMM-01_C6','SMM-01_C7','SMM-01_C0','SMM-01_C4'))
SMM01_Bcell$cluster<-Idents(SMM01_Bcell)

load('SMM-01.Bcell.cluster.markers.RData')
cluster.markers.rm.IG<-cluster.markers[!grepl('^IG.*V',cluster.markers$gene),]
top10 <- cluster.markers.rm.IG %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10<-unique(top10$gene)

# cell level
SMM01_Bcell.use<-subset(SMM01_Bcell,downsample=100)
aa<-GetAssayData(SMM01_Bcell.use,slot='scale.data')
aa<-aa[top10,]
aa.anno<-FetchData(SMM01_Bcell.use,vars = 'cluster')
aa.anno<-aa.anno%>%arrange(cluster)
aa<-aa[,rownames(aa.anno)]


diagnosis_ann_v = HeatmapAnnotation( 'cluster' = aa.anno$cluster,
                                     col = list(cluster = c('Pre B cells'='#f8766d',
                                                            'Proliferating B cells'='#de8c00',
                                                            'B cells'='#b79f00',
                                                            'SMM-01_C5'='#619cff',
                                                            'SMM-01_C2'='#00c08b',
                                                            'SMM-01_C1'='#00ba38',
                                                            'SMM-01_C3'='#00bfc4',
                                                            'SMM-01_C8'='#ff64b0',
                                                            'SMM-01_C6'='#c77cff',
                                                            'SMM-01_C7'='#f564e3',
                                                            'SMM-01_C0'='#7cae00',
                                                            'SMM-01_C4'='#00b4f0')),
                                     height=unit(2,'cm'),
                                     which='column')

gene=unique(c('DNTT','VPREB1',
              'CD9','CD99','TYMS','STMN1','IGHG1','TNFRSF17','XBP1','RRM2','TUBB4B','SPINK2','TNFRSF13B',
              'TCL1A','CD24','CD37','CD79B','SDC1','MYC','CD79A','CD27','FCRLA',
              'JCHAIN','CCND2','KLF6','CDKN1A','LAMP5','FRZB','IGHM',
              'ASS1','FOS','JNU','DUSP1','ITGB1','S100A6',
              'KLF2','ITGB7','MAFB')) 

gene_ann_v = HeatmapAnnotation(foo = anno_mark(at = match(gene,rownames(aa)), 
                                               labels = gene, side='left', 
                                               labels_gp=gpar(col = "black", fontsize = 15)),
                               which='row')

pdf('SMM-01.Bcell.seurat_clusters.TOP10.DEG.heatmap.pdf',height=9,width=16)
Heatmap(aa, name = "Z-score",
        col=circlize::colorRamp2(c(-2, 0, 2), c("blue", "black", "yellow")),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names=F,
        row_names_gp = gpar(fontsize = 6),
        #row_names_side = "bottom",
        show_column_names = F,
        column_names_side = "bottom",
        column_split=aa.anno$cluster,
        column_gap = unit(1, "mm"),
        show_heatmap_legend = T,
        use_raster=T,
        left_annotation = gene_ann_v,
        #left_annotation = left_ann,
        top_annotation = diagnosis_ann_v) 
#draw(ht_list, annotation_legend_list = lgd_list, annotation_legend_side = "top")
dev.off()

###################### Fig 5D MGUS06s(SMM-01) featureplot ######################
monocle_cds<-readRDS('SMM-01_Bcell.monocle3.rds')
gene<-c('CDR1', 'NSD2', 'ITGB7', 'CCND2', 'SPP1', 'LAMP5', 'FRZB', 'CCND1', 'XBP1',
        'KDM3A','KLF2','IRF4','CD27','JCHAIN',
        'EGFL7','ITGB3','ASS1',
        'PEBP1','GAS6','CIRBP','PERP','RRBP1','SEL1L');

pdf(paste0('SMM-01_Bcell_plot_cells.',Sys.Date(),'.pdf'),height=5,width=5.5)
for(i in gene){ #
  
  plotx<-plot_cells_log2(monocle_cds,genes=i,
                         label_groups_by_cluster=FALSE,
                         label_leaves=F,
                         label_branch_points=F,
                         graph_label_size=0,
                         group_label_size=0,
                         trajectory_graph_segment_size = 1,
                         cell_size = 1,alpha=1) +
    scale_colour_gradient2('log(Expression)',low = "#fffcfc", mid = "white",
                           high = "red", midpoint = 0, space = "Lab",#limits=c(0,3.1),
                           na.value = "#f0f0f0", guide = "colourbar", aesthetics = "colour") +
    theme(strip.background = element_blank(),strip.text.x = element_blank(),legend.title = element_blank()) +
    labs(title=i)
  
  plotx$data=plotx$data[order(plotx$data$value,na.last = F),]
  print(plotx)
}
dev.off()


gene2<-c('CKS1B','CTHRC1','LIME1','SDF2L1','LY6E','DUSP1','CD99',#'PTPRCAP',
         'SLAMF1','ITM2C','LAMP5','SULF2',
         'FGFR3','NSD2','CCND1','MAF','MAFB','CCND2',
         'MYC','CDKN2C',
         'CD38','SDC1','TNFRSF17','GPRC5D','FCRL5',
         'ITGB7', 'SPP1', 'FRZB','CADM1','LAPTM5','WFDC2','TXNIP','KLF2','IRF4','ASS1','CD47',
         'HLA-A','HLA-B','HLA-C')
pdf(paste0('SMM-01_Bcell_plot_cells2.',Sys.Date(),'.pdf'),height=5,width=5.5)
for(i in gene2){ #
  
  plotx<-plot_cells_log2(monocle_cds,genes=i,
                         label_groups_by_cluster=FALSE,
                         label_leaves=F,
                         label_branch_points=F,
                         graph_label_size=0,
                         group_label_size=0,
                         trajectory_graph_segment_size = 1,
                         cell_size = 1,alpha=1) +
    scale_colour_gradient2('log(Expression)',low = "#fffcfc", mid = "white",
                           high = "red", midpoint = 0, space = "Lab",#limits=c(0,3.1),
                           na.value = "#f0f0f0", guide = "colourbar", aesthetics = "colour") +
    theme(strip.background = element_blank(),strip.text.x = element_blank(),legend.title = element_blank()) +
    labs(title=i)
  
  plotx$data=plotx$data[order(plotx$data$value,na.last = F),]
  print(plotx)
}
dev.off()
###################### Fig 5E SMM07s(NDMM-05) B cell trajectory ######################
monocle_cds<-readRDS('NDMM-05_Bcell.monocle3.rds')

p1<-plot_cells(monocle_cds, color_cells_by = 'cell.type',
               label_groups_by_cluster=F,
               label_cell_groups=T,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=5,
               group_label_size=5,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 2,alpha=1)

p2<-plot_cells(monocle_cds, color_cells_by = 'pseudotime',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               label_roots = F,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 2,alpha=1)+
  theme(legend.position = "none")

p3<-plot_cells(monocle_cds, color_cells_by = 'infercnv.score.sd.normal_CD138_ref',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               label_roots = F,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 2,alpha=1)+
  scale_color_gradientn(colors=rev(c("#67001f","#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","white","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061","#053061")),na.value = "#f0f0f0")+
  theme(legend.position = "none")

p4<-plot_cells(monocle_cds, color_cells_by = 'BCR_clonotype_raw',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               label_roots = F,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 2,alpha=1)+ theme(legend.position = "none")

p5<-plot_cells(monocle_cds, color_cells_by = 'BCR',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               label_roots = F,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 2,alpha=1)+ theme(legend.position = "none")

aa<-read.table('./CytoTRACE/CytoTRACE_plot_table.txt',sep='\t',header = T,row.names = 1)
monocle_cds$CytoTRACE<-aa[colnames(monocle_cds),]$CytoTRACE
temp <- RColorBrewer::brewer.pal(11, "Spectral")
temp[6] <- "gold"
rbPal <- colorRampPalette(temp)(length(CytoTRACE_value))

p6<-plot_cells(monocle_cds, color_cells_by = 'CytoTRACE',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_roots = F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 1,alpha=1)+
  scale_color_gradientn(colors=rev(temp),na.value = "#f0f0f0")+
  theme(legend.position = "none")

do.call(ggarrange,c(list(p1,p2,p3,p4,p5,p6),ncol = 1,nrow = 1)) -> combined.gp
pdf(paste0('NDMM-05_Bcell_monocle3.',Sys.Date(),'.pdf'),height=5,width = 5.5)
print(combined.gp)
dev.off()

###################### Fig 5F SMM07s(NDMM-05) CD138+ inferCNV visualization ######################
library('data.table')
fdat = fread('./NDMM-05/infercnv.observations.txt', sep = ' ')
dat_obs = data.frame(t(fdat[,-1]),check.names = F)
colnames(dat_obs)=fdat$V1
dat_obs = data.frame('barcode'=rownames(dat_obs),dat_obs,check.names = F)

gene_order<-readRDS('NDMM-05.infercnv.gene_order.rds')

meta.data<-read.table('NDMM-05.CD138_pos.filter.meta.data.txt',sep='\t',header=T)

meta.data<-meta.data[,c('barcode','cluster')]

dat_obs<-dat_obs[meta.data$barcode,]
dat_obs<-left_join(meta.data[,c('barcode','cluster')],dat_obs,by=c('barcode'='barcode'))

avg.cnv.df<-dat_obs%>%
  select(-barcode)%>%
  group_by(cluster)%>%
  summarize_all(mean, na.rm=TRUE)%>%
  as.data.frame(check.names=F)


cluster.cnv.score<-data.frame('cluster'=avg.cnv.df$cluster,
                              'cnv.score'=apply(avg.cnv.df[,-1],1,sd))
rownames(cluster.cnv.score)<-cluster.cnv.score$cluster

rownames(avg.cnv.df)<-avg.cnv.df$cluster
avg.cnv.df<-avg.cnv.df[,-1]

anno.col<-gene_order[,'chr',drop=FALSE]
anno.col$chr<-factor(anno.col$chr,levels=paste0('chr',1:22))

mat_colors <- list(chr = c(rep(c("lightgrey", "darkgrey"), 11)))
names(mat_colors$chr) <- unique(anno.col$chr)

cols_gain <- c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","white")
colfunc_gain <- colorRampPalette(cols_gain)

cols_loss <- c("white","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061")
colfunc_loss <- colorRampPalette(cols_loss)

cnv.cols <- c(colfunc_gain(20),'white',colfunc_loss(round(20*(1-min(avg.cnv.df))/(max(avg.cnv.df)-1),0)))


library('ComplexHeatmap')
top_ann = HeatmapAnnotation( df=anno.col,
                             which='column',
                             annotation_name_side = "left",
                             col=mat_colors,
                             show_legend = F
)

avg.cnv.df<-avg.cnv.df[c('NDMM-05_C2','NDMM-05_C1','NDMM-05_C0','NDMM-05_C3'),]

pdf(paste0('NDMM-05.inferCNV.cluster.level.heatmap.',Sys.Date(),'.pdf'),height=2.5,width=15)
ht_list = Heatmap(avg.cnv.df, name = "inferCNV",
                  col=rev(cnv.cols),
                  cluster_rows = F,
                  cluster_columns = F,
                  show_row_names=T,
                  row_names_side = "left",
                  show_column_names = F,
                  #column_names_side = "bottom",
                  show_heatmap_legend = T,
                  
                  top_annotation = top_ann,
                  #left_annotation = left_ann,
                  #right_annotation = rowAnnotation(cnv.score = anno_points(cluster.cnv.score[rownames(avg.cnv.df),'cnv.score'])),
                  row_split = factor(rownames(avg.cnv.df),levels=rownames(avg.cnv.df)),
                  row_title = NULL,
                  column_split = anno.col$chr,
                  column_title = NULL,
                  row_gap = unit(0, "mm"), 
                  column_gap = unit(0, "mm"), 
                  border = TRUE,
                  
                  heatmap_legend_param = list(
                    #at = c(-2, 0, 2),
                    #labels = c("low", "zero", "high"),
                    #title = "Some values",
                    title_gp = gpar(fontsize = 12, fontface = "bold"),
                    labels_gp = gpar(fontsize = 12),
                    legend_height = unit(5, "cm"),
                    title_position ='topleft'
                  )) 
draw(ht_list)#, annotation_legend_list = lgd_list, annotation_legend_side = "top")
dev.off()

###################### Fig 5F SMM07s(NDMM-05) CD138+ inferCNV phylogenetic tree ######################
library('data.table')
library('ggtree')
library('tibble')
library('tidyr')
library('ape')
library('dplyr')
library('ggplot2')

fdat = fread('./NDMM-05/infercnv.observations.txt', sep = ' ')

dat_obs = data.frame(t(fdat[,-1]),check.names = F)
colnames(dat_obs)=fdat$V1

meta.data<-read.table('NDMM-05.CD138_pos.filter.meta.data.txt',sep='\t',header=T)
meta.data<-meta.data[,c('barcode','cluster')]

dat_obs = data.frame('barcode'=rownames(dat_obs),dat_obs,check.names = F)
dat_obs<-dat_obs[meta.data$barcode,]
dat_obs<-left_join(meta.data[,c('barcode','cluster')],dat_obs,by=c('barcode'='barcode'))

avg.cnv.df<-dat_obs%>%
  select(-barcode)%>%
  group_by(cluster)%>%
  summarize_all(mean, na.rm=TRUE)%>%
  as.data.frame(check.names=F)

rownames(avg.cnv.df)<-avg.cnv.df$cluster
avg.cnv.df<-avg.cnv.df[,-1]


run_me_tree <- function(consensus_df,
                        clusters,
                        ploidy_VAL,
                        rotate_nodes = NULL,
                        plot = TRUE) {
  
  consensus_int <- ploidy_VAL*consensus_df
  
  #adding a neutral state, will use as root
  consensus_int[nrow(consensus_int) + 1, ] <- round(ploidy_VAL)
  
  
  tree <- ape::fastme.bal(dist(consensus_int, method = "manhattan"))
  
  tree <-
    root.phylo(tree,
               outgroup = which(tree$tip.label == Ntip(tree)),
               resolve.root = T)
  
  tree <-
    drop.tip(tree, tip = as.character(c(
      nrow(consensus_int), nrow(consensus_int) - 1
    )))
  
  tree <- ladderize(tree)
  
  # getting order
  is_tip <- tree$edge[, 2] <= length(tree$tip.label)
  ordered_tips_index <- tree$edge[is_tip, 2]
  tree_tips_order <- tree$tip.label[ordered_tips_index] %>% rev()
  
  # adding superclones information to frequencies_df
  clones_df <- clusters %>%
    distinct(cluster) %>%
    mutate(taxa = cluster)
  
  clones_df <- clones_df[c("taxa","cluster")]
  
  if (!is.null(rotate_nodes)) {
    tree <-
      phytools::rotateNodes(tree, rotate_nodes)
    is_tip <-
      tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips_index <- tree$edge[is_tip, 2]
    tree_tips_order <-
      tree$tip.label[ordered_tips_index] %>% rev()
    tree_tips_order <-
      tree_tips_order[str_detect(tree_tips_order, "c")]
  }
  
  p <- ggtree(tree, ladderize = F, size = 2) +
    geom_treescale()
  
  p <- p %<+% clones_df +
    geom_tippoint(aes(color = cluster),
                  size = 10) 
  
  if(plot) print(p)
  
  results <- list(tree = tree,
                  cs_plot = p,
                  cs_tree_order = tree_tips_order)
  
  return(results)
  
}

pdf('NDMM-05.phylo_tree.pdf')
run_me_tree(avg.cnv.df,meta.data,2)
dev.off()

###################### Fig 5G SMM07s(NDMM-05) B cell DEG heatmap ######################
NDMM05_Bcell<-readRDS('NDMM-05.Bcell.filter.scs.rds')

Idents(NDMM05_Bcell)<-factor(NDMM05_Bcell$cell.type,levels=c('Pre B cells','Proliferating B cells','B cells','NDMM-05_C2','NDMM-05_C1','NDMM-05_C0','NDMM-05_C3'))
NDMM05_Bcell$cluster<-Idents(NDMM05_Bcell)

load('NDMM-05.Bcell.cluster.markers.RData')
cluster.markers.rm.IG<-cluster.markers[!grepl('^IG.*V',cluster.markers$gene),]
top10 <- cluster.markers.rm.IG %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10<-unique(top10$gene)

# cell level
aa<-GetAssayData(NDMM05_Bcell,slot='scale.data')
aa<-aa[top10,]
aa.anno<-FetchData(NDMM05_Bcell,vars = 'cluster')
aa.anno<-aa.anno%>%arrange(cluster)
aa<-aa[,rownames(aa.anno)]


diagnosis_ann_v = HeatmapAnnotation( 'cluster' = aa.anno$cluster,
                                     col = list(cluster = c('Pre B cells'='#a58aff',
                                                            'Proliferating B cells'='#fb61d7',
                                                            'B cells'='#f8766d',
                                                            'NDMM-05_C2'='#00c094',
                                                            'NDMM-05_C1'='#53b400',
                                                            'NDMM-05_C0'='#c49a00',
                                                            'NDMM-05_C3'='#00b6eb')),
                                     height=unit(2,'cm'),
                                     which='column')

gene=unique(c('DNTT','VPREB1',
              'CD9','CD99','TYMS','STMN1','IGHG1','TNFRSF17','XBP1','RRM2','TUBB4B','SPINK2','TNFRSF13B',
              'TCL1A','CD24','CD37','CD79B','SDC1','MYC','CD79A','CD27','FCRLA',
              'JCHAIN','CCND2','KLF6','CDKN1A','LAMP5','FRZB','IGHM','IRF4',
              'ASS1','FOS','JNU','DUSP1','ITGB1','S100A6','PRDX4','MZB1','NEB',
              'KLF2','ITGB7','MAFB')) 

gene_ann_v = HeatmapAnnotation(foo = anno_mark(at = match(gene,rownames(aa)), 
                                               labels = gene, side='left', 
                                               labels_gp=gpar(col = "black", fontsize = 18)),
                               which='row')

pdf('NDMM-05.Bcell.seurat_clusters.TOP10.DEG.heatmap.pdf',height=7,width=16)
Heatmap(aa, name = "Z-score",
        col=circlize::colorRamp2(c(-2, 0, 2), c("blue", "black", "yellow")),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names=F,
        row_names_gp = gpar(fontsize = 6),
        #row_names_side = "bottom",
        show_column_names = F,
        column_names_side = "bottom",
        column_split=aa.anno$cluster,
        column_gap = unit(1, "mm"),
        show_heatmap_legend = T,
        use_raster=T,
        left_annotation = gene_ann_v,
        #left_annotation = left_ann,
        top_annotation = diagnosis_ann_v) 
#draw(ht_list, annotation_legend_list = lgd_list, annotation_legend_side = "top")
dev.off()

###################### Fig 5H SMM07s(NDMM-05) featureplot ######################
monocle_cds<-readRDS('NDMM-05_Bcell.monocle3.rds')
gene<-c('CDR1', 'NSD2', 'ITGB7', 'CCND2', 'SPP1', 'LAMP5', 'FRZB', 'CCND1', 'XBP1',
        'KDM3A','KLF2','IRF4',
        'EGFL7','ITGB3','ASS1',
        'PEBP1','GAS6','CIRBP','PERP','RRBP1','SEL1L');

pdf(paste0('NDMM-05_Bcell_plot_cells.',Sys.Date(),'.pdf'),height=5,width=5.5)
for(i in gene){ #
  
  plotx<-plot_cells_log2(monocle_cds,genes=i,
                         label_groups_by_cluster=FALSE,
                         label_leaves=F,
                         label_branch_points=F,
                         graph_label_size=0,
                         group_label_size=0,
                         trajectory_graph_segment_size = 1,
                         cell_size = 1,alpha=1) +
    scale_colour_gradient2('log(Expression)',low = "#fffcfc", mid = "white",
                           high = "red", midpoint = 0, space = "Lab",#limits=c(0,3.1),
                           na.value = "#f0f0f0", guide = "colourbar", aesthetics = "colour") +
    theme(strip.background = element_blank(),strip.text.x = element_blank(),legend.title = element_blank()) +
    labs(title=i)
  plotx$data=plotx$data[order(plotx$data$value,na.last = F),]
  print(plotx)
}
dev.off()


gene2<-c('CKS1B','CTHRC1','LIME1','SDF2L1','LY6E','DUSP1','CD99',#'PTPRCAP',
         'SLAMF1','ITM2C','LAMP5','SULF2',
         'FGFR3','NSD2','CCND1','MAF','MAFB','CCND2',
         'MYC','CDKN2C',
         'CD38','SDC1','TNFRSF17','GPRC5D','FCRL5',
         'ITGB7', 'SPP1', 'FRZB','CADM1','LAPTM5','WFDC2','TXNIP','KLF2','IRF4','ASS1','CD47',
         'HLA-A','HLA-B','HLA-C')
pdf(paste0('NDMM-05_Bcell_plot_cells2.',Sys.Date(),'.pdf'),height=5,width=5.5)
for(i in gene2){ #
  
  plotx<-plot_cells_log2(monocle_cds,genes=i,
                         label_groups_by_cluster=FALSE,
                         label_leaves=F,
                         label_branch_points=F,
                         graph_label_size=0,
                         group_label_size=0,
                         trajectory_graph_segment_size = 1,
                         cell_size = 1,alpha=1) +
    scale_colour_gradient2('log(Expression)',low = "#fffcfc", mid = "white",
                           high = "red", midpoint = 0, space = "Lab",#limits=c(0,3.1),
                           na.value = "#f0f0f0", guide = "colourbar", aesthetics = "colour") +
    theme(strip.background = element_blank(),strip.text.x = element_blank(),legend.title = element_blank()) +
    labs(title=i)
  plotx$data=plotx$data[order(plotx$data$value,na.last = F),]
  print(plotx)
}
dev.off()