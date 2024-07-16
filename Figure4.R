#--------------------------------------------------------------
# filename : Figure4.R
# Date : 2023-09-01
# contributor : Minghao Dang, PhD
# function: 
# R version: R/4.2.0
#--------------------------------------------------------------

source('global_settings.R')

###################### Fig 4A MGUS01s B cell trajectory ######################
monocle_cds<-readRDS('M028_Bcell.monocle3.rds')
monocle_cds@int_colData@listData$reducedDims@listData$UMAP[,1]=-monocle_cds@int_colData@listData$reducedDims@listData$UMAP[,1]
monocle_cds <- learn_graph(monocle_cds,
                           close_loop = F,
                           learn_graph_control = list(minimal_branch_len=8,rann.k=NULL))
library('monocle3')

p1<-plot_cells(monocle_cds, color_cells_by = 'cell.type',
               label_groups_by_cluster=F,
               label_cell_groups=T,
               label_leaves=F,
               label_roots = F,
               label_branch_points=F,
               graph_label_size=10,
               group_label_size=10,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 1,alpha=1) + scale_fill_manual(values =(c("CD34+" = "#f8766d", 
                                                                     'MKI67+'='#cd9600',
                                                                     'Mature'='#7cae00',
                                                                     'CD138+_C3'='#00be67',
                                                                     'CD138+_C5'='#00bfc4',
                                                                     'CD138+_C4'='#00a9ff',
                                                                     'CD138+_C2'='#c77cff',
                                                                     'CD138+_C1'='#ff61cc')))
#p1<-AugmentPlot(p1)
ggsave(paste0('M028_Bcell_monocle3.cell.type.',Sys.Date(),'.pdf'),p1,width=5.5,height=5)


`%notin%`<-Negate(`%in%`)
monocle_cds$BCR<-monocle_cds$BCR_clonotype_raw
monocle_cds$BCR[monocle_cds$BCR_clonotype_raw%notin%paste0('clonotype',1,'_M028')]<-NA
p1<-plot_cells(monocle_cds, color_cells_by = 'BCR',
               label_groups_by_cluster=FALSE,
               label_cell_groups=FALSE,
               label_leaves=F,
               label_roots = F,
               label_branch_points=F,
               graph_label_size=3,
               group_label_size=3,
               trajectory_graph_color = "black",
               trajectory_graph_segment_size = 1.5,
               cell_size = 2,alpha=1)+ theme(legend.position = "none")
#p1<-AugmentPlot(p1)
ggsave(paste0('M028_Bcell_monocle3.top1.clonotype.',Sys.Date(),'.pdf'),p1,width=5.5,height=5)


p1<-plot_cells(monocle_cds, color_cells_by = 'infercnv.score.sd.normal_CD138_ref',
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
  scale_color_gradientn(colors=rev(c("#67001f","#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","white","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061","#053061")),na.value = "#f0f0f0")+
  theme(legend.position = "none")
#scale_color_gradientn(colors=viridis(20),na.value = "#f0f0f0")
#p1<-AugmentPlot(p1)
ggsave(paste0('M028_Bcell_monocle3.CNV_score.',Sys.Date(),'.pdf'),p1,width=5.5,height=5)


monocle_cds<-order_cells(monocle_cds,root_cells = c('TAGCCGGTCCTTAATC-1_M028'))
p1<-plot_cells(monocle_cds, color_cells_by = 'pseudotime',
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
  theme(legend.position = "none")
#scale_color_gradientn(colors=viridis(20),na.value = "#f0f0f0")
#p1<-AugmentPlot(p1)
ggsave(paste0('M028_Bcell_monocle3.pesudotime.',Sys.Date(),'.pdf'),p1,width=5.5,height=5)


aa<-read.table('./CytoTRACE/CytoTRACE_plot_table.txt',sep='\t',header = T,row.names = 1)
monocle_cds$CytoTRACE<-aa[colnames(monocle_cds),]$CytoTRACE

temp <- RColorBrewer::brewer.pal(11, "Spectral")
temp[6] <- "gold"
rbPal <- colorRampPalette(temp)(length(CytoTRACE_value))

p1<-plot_cells(monocle_cds, color_cells_by = 'CytoTRACE',
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
ggsave(paste0('M028_Bcell_monocle3.CytoTRACE.',Sys.Date(),'.pdf'),p1,width=5.5,height=5)

###################### Fig 4B MGUS01s B cnv score vlnplot ######################
M028_Bcell<-readRDS('M028.Bcell.filter.scs.rds')
Idents(M028_Bcell)<-factor(M028_Bcell$cluster,levels=c('CD34+','MKI67+','Mature','CD138+_C3','CD138+_C5','CD138+_C4','CD138+_C2','CD138+_C1'))

M028_Bcell.sub<-subset(M028_Bcell,cluster%in%c('CD138+_C3','CD138+_C5','CD138+_C4','CD138+_C2','CD138+_C1'))
plotx<-VlnPlot(M028_Bcell.sub,features='infercnv.score.sd.normal_CD138_ref',pt.size=0)

ggsave(paste0('M028.CD138_pos.cnv.score.vlnplot.',Sys.Date(),'.pdf'),plotx,height=6,width=6)

###################### Fig 4C MGUS01s B cell clonotype cluster level stack plot ######################
M028_Bcell<-readRDS('M028.Bcell.filter.scs.rds')
Idents(M028_Bcell)<-factor(M028_Bcell$cluster,levels=c('CD34+','MKI67+','Mature','CD138+_C3','CD138+_C5','CD138+_C4','CD138+_C2','CD138+_C1'))

bcr<-FetchData(M028_Bcell,vars=c('cluster','BCR_clonotype_raw'))

bcr$BCR_clonotype_group<-bcr$BCR_clonotype_raw
bcr$BCR_clonotype_group[bcr$BCR_clonotype_raw%in%paste0('clonotype',1,'_M028')]<-'Top 1'
bcr$BCR_clonotype_group[bcr$BCR_clonotype_raw%in%paste0('clonotype',2:5,'_M028')]<-'Top 2-5'
bcr$BCR_clonotype_group[bcr$BCR_clonotype_raw%in%paste0('clonotype',6:20,'_M028')]<-'Top 6-20'
bcr$BCR_clonotype_group[bcr$BCR_clonotype_raw%in%paste0('clonotype',21:100,'_M028')]<-'Top 21-100'
bcr$BCR_clonotype_group[bcr$BCR_clonotype_raw%in%paste0('clonotype',101:500,'_M028')]<-'Top 101-500'

data1=as.data.frame(table(bcr[,c('cluster','BCR_clonotype_group')]))
data2 = as.data.frame(table(bcr$cluster))
colnames(data2)<-c('cluster','count')
data3<-inner_join(data1,data2,by='cluster')
data3$percentage<- data3$Freq*100/data3$count;

data3$cluster<-factor(data3$cluster,levels=c('CD34+','MKI67+','Mature','CD138+_C3','CD138+_C5','CD138+_C4','CD138+_C2','CD138+_C1'))
data3$BCR_clonotype_group<-factor(data3$BCR_clonotype_group,levels=c('Top 1','Top 2-5','Top 6-20','Top 21-100','Top 101-500'))


color=c('#fee5d9','#fcae91','#fb6a4a','#de2d26','#a50f15')
names(color)=rev(c('Top 1','Top 2-5','Top 6-20','Top 21-100','Top 101-500'))

plotx<-ggplot(data3, aes(x=cluster,fill=BCR_clonotype_group,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values = color) +
  labs(x='',y = 'clonotype proportion') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 16),
        #legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))
ggsave('M028.Bcell.BCR.clonotype.stackbarplot.pdf',useDingbats=F,plotx,width = 7)

###################### Fig 4D MGUS01s CD138+ inferCNV visualization ######################
library('data.table')
fdat = fread('infercnv.observations.txt', sep = ' ')
dat_obs = data.frame(t(fdat[,-1]),check.names = F)
colnames(dat_obs)=fdat$V1
dat_obs = data.frame('barcode'=rownames(dat_obs),dat_obs,check.names = F)

infercnv_obj=readRDS('run.final.infercnv_obj')
gene_order=infercnv_obj@gene_order
gene_order$chr=as.character(gene_order$chr)

meta.data<-read.table('M028.CD138_pos.filter.meta.data.2021-08-12.txt',sep='\t',header=T)
meta.data$cluster<-'CD138+_C1'
meta.data$cluster[meta.data$seurat_clusters==1]<-'CD138+_C2'
meta.data$cluster[meta.data$seurat_clusters==2]<-'CD138+_C3'
meta.data$cluster[meta.data$seurat_clusters==3]<-'CD138+_C4'
meta.data$cluster[meta.data$seurat_clusters==4]<-'CD138+_C5'

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

# anno.row<-unique(meta.data[,c('cluster','treatment')])
# rownames(anno.row)<-anno.row$cluster
# anno.row<-anno.row[rownames(avg.cnv.df),'treatment',drop=FALSE]
# 
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


# left_ann = HeatmapAnnotation( df=anno.row,
#                               which='row',
#                               show_annotation_name = F,
#                               #annotation_name_side = NULL,
#                               col=mat_colors,
#                               #simple_anno_size=unit(1, "cm"),
#                               show_legend = F
# )

# lgd_list = list(
#   Legend(labels = c('pre','on','post'),
#          title = "treatment", type = "grid", pch = 16,
#          grid_height = unit(6, "mm"), grid_width = unit(6, "mm"), gap = unit(3, "mm"),
#          labels_gp = gpar(fontsize = 15),
#          title_gp = gpar(fontsize = 15, fontface = "bold"),
#          legend_gp = gpar(fill=c('grey','lightblue','#756bb1')),
#          title_position ='leftcenter',
#          nrow=1)
# )

avg.cnv.df<-avg.cnv.df[c('CD138+_C3','CD138+_C5','CD138+_C4','CD138+_C2','CD138+_C1'),]
pdf(paste0('M028.inferCNV.cluster.level.heatmap.',Sys.Date(),'.pdf'),height=2,width=15)
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
                  row_split = factor(rownames(avg.cnv.df),levels=c('CD138+_C1','CD138+_C2','CD138+_C3','CD138+_C4','CD138+_C5')),
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



pdf(paste0('inferCNV.cluster.level.heatmap.',Sys.Date(),'.pdf'),height=2,width=25)
Heatmap(avg.cnv.df[,rownames(anno.col)[anno.col=='chr9']], name = "inferCNV",
        col=rev(cnv.cols),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names=T,
        row_names_side = "left",
        show_column_names = T,
        column_names_gp = grid::gpar(fontsize = 6),
        #column_names_side = "bottom",
        show_heatmap_legend = T,
        row_title = NULL,
        #column_split = anno.col$chr,
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
dev.off()


aa<-avg.cnv.df['CD138+_C1',rownames(anno.col)[anno.col=='chr19']]
aa<-aa[which(aa>1.05)]

library('clusterProfiler')
library('msigdbr')

m_t2g_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% # hallmark gene sets
  select(gs_name, gene_symbol)

m_t2g<-m_t2g_hallmark

res<-enricher(names(aa),TERM2GENE=m_t2g,pvalueCutoff = 1)
write.table(names(aa),'aa.txt',sep='\t',row.names = F,col.names = F,quote = F)

res$GeneRatio2<-sapply(res$GeneRatio,function(x) {eval(parse(text=x)) })
res$FDR<- -log10(res$qvalue)

plotx<-ggplot(res, aes(x = cluster,y = ID)) +        ## global aes
  geom_point(aes(fill = FDR,size =GeneRatio2),color='black',shape=21)  +    ## geom_point for circle illusion
  #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
  #scale_fill_viridis()+
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b'),limits=c(0,max(pathway_er$FDR))) +
  scale_size(range = c(0,6), limits = c(0, 0.3), breaks = c(0,0.05,0.1,0.15,0.2,0.25,0.3)) +             ## to tune the size of circles
  #facet_grid(cols=vars(gene_function),scales='free',space='free')+
  labs(x='',y='')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), #panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text = element_text(angle = 0,size=12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95) )
ggsave(paste0('HPD.vs.translocation.pathway.enrichment.BubblePlot.',Sys.Date(),'.pdf'),plotx,height=5,width=5)

###################### Fig 4E MGUS01s B cell DEG heatmap ######################
M028_Bcell<-readRDS('M028.Bcell.filter.scs.rds')
Idents(M028_Bcell)<-factor(M028_Bcell$cluster,levels=c('CD34+','MKI67+','Mature','CD138+_C3','CD138+_C5','CD138+_C4','CD138+_C2','CD138+_C1'))
M028_Bcell$cluster<-Idents(M028_Bcell)

load('M028.Bcell.cluster.markers.RData')
cluster.markers$cluster<-factor(cluster.markers$cluster,levels=c('CD34+','MKI67+','Mature','CD138+_C3','CD138+_C5','CD138+_C4','CD138+_C2','CD138+_C1'))
cluster.markers<-cluster.markers%>%arrange(cluster)
top50 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.table(top50,file=paste0('M028','.Bcell.cluster.top50.markers.txt'),quote=F,sep="\t",row.names = F,col.names = T)

cluster.markers.rm.IG<-cluster.markers[!grepl('^IG.*V',cluster.markers$gene),]

top10 <- cluster.markers.rm.IG %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10<-unique(top10$gene)


# cell level
aa<-GetAssayData(M028_Bcell,slot='scale.data')
aa<-aa[top10,]
aa.anno<-FetchData(M028_Bcell,vars = 'cluster')
aa.anno<-aa.anno%>%arrange(cluster)
aa<-aa[,rownames(aa.anno)]


diagnosis_ann_v = HeatmapAnnotation( 'cluster' = aa.anno$cluster,
                                     col = list(cluster = c("CD34+" = "#f8766d", 
                                                            'MKI67+'='#cd9600',
                                                            'Mature'='#7cae00',
                                                            'CD138+_C3'='#00be67',
                                                            'CD138+_C5'='#00bfc4',
                                                            'CD138+_C4'='#00a9ff',
                                                            'CD138+_C2'='#c77cff',
                                                            'CD138+_C1'='#ff61cc')),
                                     height=unit(2,'cm'),
                                     which='column')

gene=unique(c('DNTT','VPREB1',
              'CD9','CD99','TYMS','STMN1',
              'TCL1A','CD24','CD37','CD79B','SDC1','MZB1',
              'JCHAIN','CCND2','KLF6','CDKN1A',
              'ASS1','FOS','JNU','DUSP1',
              'KLF2','ITGB7','MAF',
              'CKS1B','CTHRC1','LIME1','SDF2L1','LY6E','PTPRCAP',
              'SLAMF1','FAM30A','ITM2C','LAMP5','SULF2',
              'FGFR3','NSD2','CCND1','MAFB',
              'MYC','CDKN2C',
              'CD38','TNFRSF17','GPRC5D','FCRL5',
              'FRZB','CADM1','LAPTM5','WFDC2','TXNIP','IRF4','CD47',
              'HLA-A','HLA-B','HLA-C')) 

gene_ann_v = HeatmapAnnotation(foo = anno_mark(at = match(gene,rownames(aa)), 
                                               labels = gene, side='left', 
                                               labels_gp=gpar(col = "black", fontsize = 15)),
                               which='row')

pdf('M028.Bcell.seurat_clusters.TOP10.DEG.heatmap.pdf',height=8,width=16)
Heatmap(aa, name = "Z-score",
        col=circlize::colorRamp2(c(-2, 0, 2), c("blue", "black", "yellow")),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names=F,
        #row_names_side = "bottom",
        show_column_names = F,
        column_names_side = "bottom",
        column_split=aa.anno$cluster,
        column_gap = unit(1, "mm"),
        show_heatmap_legend = T,
        left_annotation = gene_ann_v,
        use_raster = T,
        #left_annotation = left_ann,
        top_annotation = diagnosis_ann_v) 
#draw(ht_list, annotation_legend_list = lgd_list, annotation_legend_side = "top")
dev.off()

###################### Fig 4F MGUS01s C0 vd C1 deg volcano ######################
M028_Bcell<-readRDS('M028.Bcell.filter.scs.rds')

cluster.markers<-FindMarkers(M028_Bcell,min.pct = 0.25,ident.1='CD138+_C0',ident.2='CD138+_C1',logfc.threshold = 0)
cluster.markers$symbol<-rownames(cluster.markers)

cluster.markers$note<-ifelse(cluster.markers$avg_log2FC>0,'High in C0','High in C1')
write.table(cluster.markers,file=paste0('M028','.C0.vs.C1.txt'),quote=F,sep="\t",row.names = F,col.names = T)

cluster.markers$log_p_val_adj= -log10(cluster.markers$p_val_adj)
result=cluster.markers[,c('avg_log2FC','log_p_val_adj')]
colnames(result)=c('log2FoldChange','log10fdr')
color<-rep('other',dim(result)[1])
color[result$log2FoldChange > 0.3 & result$log2FoldChange <= 0.5 & result$log10fdr > 4]<-'High in C1';
color[result$log2FoldChange > 0.5 & result$log10fdr > 2]<-'High in C1';
color[result$log2FoldChange >= -0.5 & result$log2FoldChange < -0.3 & result$log10fdr > 4]<-'High in C2';
color[result$log2FoldChange < -0.5 & result$log10fdr > 2]<-'High in C2';
result$color<-color;
#gn.selected <- (result$log2FoldChange > 0.5 & result$log10fdr > 2) | (result$log2FoldChange < -0.5 & result$log10fdr > 2) | (result$log2FoldChange < 0.5 & result$log2FoldChange >= 0.3 & result$log10fdr > 4) | (result$log2FoldChange < -0.3 & result$log2FoldChange >= -0.5 & result$log10fdr > 4)
gn.selected <- result$color!='other'
# gn.selected <- rownames(result)%in%c('IL7R','FOS','LTB','GZMK','JUNB','CXCL13',
#                                      'GZMB',
#                                      'ENTPD1',
#                                      'HLA-DRB1','HLA-DRA',
#                                      'HLA-DQA1','ITGAE','CTLA4','HLA-DPA1','TNFRSF9','LAG3','CD74')


plotx<-ggplot(result,aes(x=log2FoldChange, y=log10fdr)) +
  geom_point(aes(fill=color,size=log10fdr),color='black',shape=21) +
  scale_size(range = c(0,4), limits = c(0, 50), breaks = c(0,10,20,50)) +
  scale_fill_manual(values = c('red','blue','grey')) +
  ggtitle('C1 vs C2') +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  #ylim(0,150) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = subset(result, gn.selected),aes(label = rownames(subset(result, gn.selected))),size = 3,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  max.overlaps=20) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0('M028.CD138_C1.vs.CD138_C2.Volcano.plot.',Sys.Date(),'.pdf'),plotx,width = 6,height = 6.5)


###################### Fig 4G validation of upregulated genes in bulk ######################
fpkm_log2_combat<-read.csv('5-rm-samples-combat-adjusted-log2_tpm-2021-04-25.csv',row.names = 1,check.names = F)


Anno<-readRDS('Bulk.Anno.rds')
#Anno.use<-subset(Anno,study=='MS54')
Anno.use<-Anno

Anno_color=readRDS('Bulk.Anno_color.rds')

gene<-c('ASS1',
        'CST6',
        'PTP4A3',
        'PRPSAP2',
        'CCDC144A',
        'IL15',
        'IL15RA')

aa=fpkm_log2_combat[gene,]
aa=aa[,which(colnames(aa)%in%Anno.use$sampleID)]

library('reshape2')
aa<-melt(t(aa))
colnames(aa)=c('sampleID','gene','value')

library('dplyr')
aa=left_join(aa,Anno.use,by='sampleID')
aa$progression[aa$progression!='NonPD']<-'PD'

sel1=which(aa$CD138=='Pos' & aa$point=='baseline')
sel2=which(aa$CD138=='Pos' & aa$progression=='PD')
sel3=which(aa$CD138=='Neg' & aa$point=='baseline')
sel4=which(aa$CD138=='Neg' & aa$progression=='PD')

library('ggsignif')
pdf('M028.gene.in.bulk.boxplot.pdf',width=15,height=5)
p1<-ggplot(aa[sel1,], aes(x=factor(progression), y=value) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T, aes(fill=factor(progression)), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.3, height=0),aes(colour=factor(batch)), alpha=0.9) +
  geom_signif(test="wilcox.test", comparisons = list(levels(factor(aa$progression))),
              map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
              vjust=0.5,
              textsize=4,
              size=0.5,
              step_increase = 0.06)+
  facet_wrap(.~gene,ncol=length(gene)) +
  theme(strip.text.x = element_text(size=12, color="black", face="bold"),text = element_text(size=20)) +
  labs(x="", y="Batch corrected log2(TPM+1)", fill='CD138+ Baseline', color='CD138+ Baseline')
#ggsave('/Users/mdang1/Box/MS54/outputs/compare1.gene.exp.pdf',P,width=20,height=5)

p2<-ggplot(aa[sel2,], aes(x=factor(point), y=value) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T, aes(fill=factor(point)), outlier.shape=NA) +
  geom_jitter(position=position_jitter(width=0.3, height=0),aes(colour=factor(batch)), alpha=0.9) +
  geom_signif(test="wilcox.test", comparisons = list(levels(factor(aa$point))),
              map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
              vjust=0.5,
              textsize=4,
              size=0.5,
              step_increase = 0.06)+
  facet_wrap(.~gene,ncol=length(gene)) +
  theme(strip.text.x = element_text(size=12, color="black", face="bold"),text = element_text(size=20)) +
  labs(x="", y="Batch corrected log2(TPM+1)", fill='CD138+ Progression', color='CD138+ Progression')
#ggsave('/Users/mdang1/Box/MS54/outputs/compare2.gene.exp.pdf',P,width=20,height=5)

p3<-ggplot(aa[sel1,], aes(x=factor(Diagnosis_at_Presentation), y=value) ) + 
  geom_boxplot(alpha = 0.4, show.legend = T, aes(fill=factor(Diagnosis_at_Presentation)), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.3, height=0),aes(colour=factor(Diagnosis_at_Presentation))) +
  geom_signif(test="wilcox.test", comparisons = list(levels(factor(aa$Diagnosis_at_Presentation))),
              map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
              vjust=0.5,
              textsize=4,
              size=0.5,
              step_increase = 0.06)+
  facet_grid(rows=vars(CD138),cols=vars(gene),scale='free') +
  theme(strip.text.x = element_text(size=12, color="black", face="bold"),text = element_text(size=20)) +
  labs(x="", y="Expression") 

aa$group<-paste0(aa$Diagnosis_at_Presentation,'_',aa$progression,'_',aa$point)
aa$group[aa$group%in%c('MGUS_NonPD_baseline','MGUS_NonPD_end')]<-'MGUS'
aa$group[aa$group%in%c('SMM_NonPD_baseline','SMM_NonPD_end','SMM_PD_baseline')]<-'SMM'
aa$group[aa$group%in%c('SMM_PD_end')]<-'MM'
aa$group<-factor(aa$group,levels=c('MGUS','SMM','MM'))

p4<-ggplot(aa[aa$CD138=='Pos',], aes(x=factor(group), y=value) ) + 
  geom_boxplot(alpha = 0.4, show.legend = T, aes(fill=factor(group)), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.3, height=0)) + # ,aes(colour=factor(study))
  scale_fill_manual(values=c('MGUS'='#1c75bc','SMM'='#fbb040','MM'='#d01c8b')) +
  geom_signif(test="wilcox.test", comparisons = combn(levels(aa$group),2, simplify = F),
              map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
              vjust=0.5,
              textsize=4,
              size=0.5,
              step_increase = 0.06)+
  facet_grid(rows=vars(CD138),cols=vars(gene),scale='free') +
  theme(axis.text.x = element_text(size = 12,angle=90,vjust=0.2,hjust=0.95),
        strip.text.x = element_text(size=12, color="black", face="bold"),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1)) +
  labs(x="", y="Expression") 

p5<-ggplot(aa, aes(x=factor(Diagnosis_at_Presentation), y=value) ) + 
  geom_boxplot(alpha = 0.4, show.legend = T, aes(fill=factor(Diagnosis_at_Presentation)), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.3, height=0),aes(colour=factor(Diagnosis_at_Presentation))) +
  #geom_line(aes(group = id)) +
  geom_signif(test="wilcox.test", comparisons = list(levels(factor(aa$Diagnosis_at_Presentation))),
              map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
              vjust=0.5,
              textsize=4,
              size=0.5,
              step_increase = 0.06)+
  facet_grid(rows=vars(CD138),cols=vars(gene),scale='free') +
  theme(strip.text.x = element_text(size=12, color="black", face="bold"),text = element_text(size=20)) +
  labs(x="", y="Expression") 

print(p1);print(p2);print(p3);print(p4);
dev.off()

p1+p2+p3+p4+p5

###################### Fig 4H featureplot ######################
library('monocle3')

monocle_cds<-readRDS('M028_Bcell.monocle3.rds')

monocle_cds@int_colData@listData$reducedDims@listData$UMAP[,1]=-monocle_cds@int_colData@listData$reducedDims@listData$UMAP[,1]
monocle_cds <- learn_graph(monocle_cds,
                           close_loop = F,
                           learn_graph_control = list(minimal_branch_len=8,rann.k=NULL))

gene<-unique(c('CDR1', 'NSD2', 'ITGB7', 'CCND2', 'SPP1', 'LAMP5', 'FRZB', 'CCND1', 'ASS1',
               'KDM3A','KLF2','IRF4',
               'EGFL7','ITGB3',
               'PEBP1','GAS6','CIRBP','PERP','RRBP1','SEL1L',
               'CKS1B','CTHRC1','LIME1','SDF2L1','LY6E','DUSP1','CD99','PTPRCAP',
               'SLAMF1','FAM30A','ITM2C','LAMP5','SULF2',
               'FGFR3','NSD2','CCND1','MAF','MAFB','CCND2',
               'MYC','CDKN2C',
               'CD38','SDC1','TNFRSF17','GPRC5D','FCRL5',
               'ITGB7', 'SPP1', 'FRZB','CADM1','LAPTM5','WFDC2','TXNIP','KLF2','IRF4','ASS1','CD47',
               'HLA-A','HLA-B','HLA-C'));

pdf(paste0('M028_Bcell_plot_cells.',Sys.Date(),'.pdf'),height=5,width=5.5)
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
    #scale_fill_viridis(limits=c(0,3.1)) +
    #scale_color_gradientn(colors=c(viridis(3)[1],viridis(21)[c(9,10,11,12,20)],viridis(3)[3]),na.value = "#f0f0f0")+
    #scale_color_gradientn(colors=jet()[c(1:11,15:25)],na.value = "#e6e6e6")+
    theme(strip.background = element_blank(),strip.text.x = element_blank(),legend.title = element_blank()) +
    labs(title=i)
  # leg <- get_legend(plotx)
  plotx$data=plotx$data[order(plotx$data$value,na.last = F),]
  # plotx=AugmentPlot(plotx)
  #
  # aa=ggarrange(as.ggplot(plotx),as.ggplot(leg),ncol=2)
  print(plotx)
}
dev.off()

###################### Fig 4I MGUS01s B driver gene vlnplot ######################
M028_Bcell<-readRDS('M028.Bcell.filter.scs.rds')
Idents(M028_Bcell)<-factor(M028_Bcell$cluster,levels=c('CD34+','MKI67+','Mature','CD138+_C3','CD138+_C5','CD138+_C4','CD138+_C2','CD138+_C1'))

gene<-c('EGFL7','NSD2','LAMP5','IRF4','KLF2','ITGB7','MAF','CCND2','ASS1')

plotx<-VlnPlot(M028_Bcell,features=gene,pt.size=0,ncol=1)&
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),#axis.ticks=element_blank(),
        #strip.text.x = element_text(size=12, color="black", face="bold"),
        strip.background = element_blank(),strip.text = element_text(angle = 0),
        #strip.text.y = element_text(angle = 0),
        text = element_text(size=15),
        plot.title = element_blank(),#element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1))
ggsave(paste0('M028.Bcell.driver.gene.vlnplot.',Sys.Date(),'.pdf'),plotx,height=12,width=3)
###################### Fig 4J ######################
M028_Bcell<-readRDS('M028.Bcell.filter.scs.rds')
monocle_cds<-readRDS('M028_Bcell.monocle3.rds')
pseudotime=pseudotime(monocle_cds)
M028_Bcell$pseudotime<-pseudotime[colnames(M028_Bcell)];

gene<-c('CDR1', 'NSD2', 'ITGB7', 'CCND2', 'SPP1', 'LAMP5', 'FRZB', 'CCND1', 'ASS1',
        'KDM3A','KLF2','IRF4','MAF',
        'EGFL7','ITGB3');
gene<-intersect(gene,rownames(M028_Bcell));
res<-NULL;
for(i in 1:length(gene)){
  exp1=as.matrix(GetAssayData(M028_Bcell)[gene[i],]);
  plot.data <- data.frame(pse_time=M028_Bcell$pseudotime,exp=exp1,gene=gene[i]);
  plot.data<-plot.data[plot.data$pse_time != Inf,]
  res<-rbind(res,plot.data)
}

res<-subset(res,gene%in%c('EGFL7','NSD2','LAMP5','IRF4','KLF2','ITGB7','MAF','CCND2','ASS1'))

res$gene<-factor(res$gene,levels=c(#'CCND1','CDR1','FRZB','ITGB3','KDM3A',
  'EGFL7','NSD2','LAMP5','IRF4','KLF2','ITGB7','MAF','CCND2','ASS1'))
p2 <- ggplot(res,aes(x = pse_time, y = exp)) +
  #scale_color_manual(values = color.gene)+
  geom_smooth() +
  ylab('expression')+
  xlab('pseudotime')+
  facet_grid(rows=vars(gene),scales='free')+
  #theme_classic()+
  theme_bw()
ggsave('M028_Bcell.monocle3.pseudotime.line.marker.genes.pdf',p2,width =3.5,height=12)


###################### Fig 4K MGUS01s B cell hallmark pathway gsea ######################
library('msigdbr')
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% # hallmark gene sets
  select(gs_name, gene_symbol) %>% as.data.frame()

pathway<-hallmark
pathway<-split(pathway,pathway$gs_name)
pathway<-lapply(pathway,function(x){x$gene_symbol})


M028_Bcell<-readRDS('M028.Bcell.filter.scs.rds')

library('GSVA')
data <- GetAssayData(M028_Bcell)
gsva.pathway<-gsva(as.matrix(data),pathway,annotation,method=c('ssgsea'))

gsva.score<-data.frame('barcode'=colnames(gsva.pathway),t(gsva.pathway))

meta.data<-M028_Bcell@meta.data
meta.data<-data.frame('barcode'=rownames(meta.data),meta.data)

gsva.score<-left_join(meta.data[,c('barcode','cluster')],gsva.score,by=c('barcode'='barcode'))

gsva.score$cluster<-factor(gsva.score$cluster,levels=c('CD34+','MKI67+','Mature','CD138+_C5','CD138+_C3','CD138+_C4','CD138+_C2','CD138+_C1'))

gsva.score$barcode<-NULL

gsva.score<-melt(gsva.score)

P<-ggplot(gsva.score,aes(x=cluster,y=value,fill=cluster)) +
  geom_violin() + # adjust=1.5,bw=0.3,scale='area',color='white'
  stat_summary(fun.y=median,geom='point',
               position = position_dodge(width = 0.9)) +
  # scale_fill_manual(values = c('hot'='#d73027','warm'='#FFC53F','cold'='#4575b4')) +
  # scale_color_manual(values = c('hot'='#d73027','warm'='#FFC53F','cold'='#4575b4')) +
  # geom_signif(test="kruskal.test", comparisons = list(c('cold','warm','hot')),#combn(levels(gene.to.check.exp_value$symbol),2, simplify = F),
  #             map_signif_level = T,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
  #             vjust=0.5,
  #             textsize=4,
  #             size=0.5,
  #             step_increase = 0.06)+
  stat_compare_means(method = "kruskal.test")+
  facet_wrap(vars(variable),nrow=5,ncol=10,scales='free')+
  theme(strip.text.x = element_text(size=10, color="black", face="bold"),
        axis.text.x = element_blank(), #element_text(angle = 90,size = 15,vjust=0.5,hjust=0.95),
        text = element_text(size=20),
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x='', y="ssGSEA score")

ggsave('M028.B_cell.hallmark.pathway.gsea.violin.pdf',useDingbats=F,P,width=40,height=15)
