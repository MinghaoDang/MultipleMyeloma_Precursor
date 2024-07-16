library('harmony')
library('dplyr')
library('Seurat')
library('ggpubr')
library('RColorBrewer')
library('cowplot')
library('ggplot2')
library('pheatmap')
library('ggrepel')
library('DoubletFinder')
library('ggplotify')
library('viridis')
library('readxl')
library('ComplexHeatmap')
library('monocle3')
`%notin%`<-Negate(`%in%`)
options(warn=0)

setwd('/rsrch3/scratch/genomic_med/mdang1/Elisabet/SingleCell/seurat')

theme_umap<-function(){
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),axis.ticks=element_blank(),
        #strip.text.x = element_text(size=12, color="black", face="bold"),
        strip.background = element_blank(),strip.text = element_text(angle = 0),
        #strip.text.y = element_text(angle = 0),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA, size=1))
}

theme_boxplot<-function(){
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        #plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20))
}

theme_vlnplot<-function(){
  theme(strip.text.x = element_text(size=12, color="black", face="bold"),
        strip.text.y = element_text(size=12, color="black",angle = 0),
        axis.text.x = element_text(angle = 90,size = 12,vjust=0.5,hjust=0.95),
        text = element_text(size=20),
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA))
}

batch_feature_plot=function(gene.symbol,pbmc,out){
  
  gene.use<-intersect(gene.symbol,rownames(pbmc))
  if(length(gene.use)>0){
    pdf(out,height=8,width=9)
    for(i in gene.use){
      plotx<-FeaturePlot(object = pbmc,
                         pt.size = 0.5,
                         features = i,#max.cutoff = 3,
                         cols = c("grey", "red"),
                         reduction = "umap",
                         order=F,
                         combine = T
      )#+scale_color_gradientn(colors=jet(),na.value = "grey50");
      # leg <- get_legend(plotx)
      # plotx=AugmentPlot(plotx)
      # plotx=ggarrange(as.ggplot(plotx),as.ggplot(leg),ncol=2)
      print(plotx)
    }
    dev.off()  
  }
  
}

# cell.type.level1.color<-c('#1f78b4','#33a02c','#a6cee3','#ff7f00','#b15928','#e31a1c','#6a3d9a','#ffed6f','#525252')
# names(cell.type.level1.color)<-c('Mature B cells','T cells','NK','Myeloid','Endothelial','Fibroblast','Platelet','Immature/Progeniors','CD138_pos')
# saveRDS(cell.type.level1.color,'cell.type.level1.color.rds')
cell.type.level1.color<-readRDS('cell.type.level1.color.rds')
cell.type.level2.color<-readRDS('cell.type.level2.color.rds')