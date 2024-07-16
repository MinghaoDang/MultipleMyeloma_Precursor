#--------------------------------------------------------------
# filename : Figure6.R
# Date : 2023-09-01
# contributor : Minghao Dang, PhD
# function: 
# R version: R/4.2.0
#--------------------------------------------------------------

source('global_settings.R')

###################### Fig 6A TME UMAP ######################
pbmc.f.TME<-readRDS('MM.TME.filter.scs.rds')

cell.type.level2.color<-readRDS('cell.type.level2.color.rds')
names(cell.type.level2.color)[10]<-'cDC2'
cell.type.level2.color<-c(cell.type.level2.color,'pDC'='#dc9ae9')
saveRDS(cell.type.level2.color,'cell.type.level2.color.rds')

p1<-DimPlot(object = pbmc.f.TME, reduction = "umap",pt.size = 0.5,group.by = 'orig.ident',label = F,label.size = 4)+theme(legend.position = "none")


p2<-DimPlot(object = pbmc.f.TME, reduction = "umap",pt.size = 0.5,group.by = 'cell.type.level2',label = F,label.size = 4)+ scale_color_manual(values=cell.type.level2.color)
p3<-DimPlot(object = pbmc.f.TME, reduction = "umap",pt.size = 0.5,group.by = 'cell.type.level2',label = F,label.size = 4)+ scale_color_manual(values=cell.type.level2.color)+
  theme(legend.position = "none")
p4<-DimPlot(object = pbmc.f.TME, reduction = "umap",pt.size = 0.5,group.by = 'cell.type.level2',label = F,label.size = 4,raster = F)+ scale_color_manual(values=cell.type.level2.color)+
  theme(legend.position = "none")

p5<-DimPlot(object = pbmc.f.TME, reduction = "umap",pt.size = 0.5,group.by = 'diagnosis',label = F,label.size = 4)+ scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'))
p6<-DimPlot(object = pbmc.f.TME, reduction = "umap",pt.size = 0.5,group.by = 'diagnosis',label = F,label.size = 4)+ scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'))+
  theme(legend.position = "none")

p7<-DimPlot(object = pbmc.f.TME, reduction = "umap",pt.size = 0.5,group.by = 'FISH_2',label = F,label.size = 4)+ 
  scale_color_manual(values=c('Normal'='#00a651','t(4;14)'='#8dd3c7','t(11;14)'='#bebada','t(14;16)'='#80b1d3','t(14;20)'='#fccde5','Hyperdiploidy'='#ff7f00','Other'='#7f3f98','Not available'='#878787'))+
  theme(legend.position = "none")

do.call(ggarrange,c(list(p1,p2,p3,p4,p5,p6,p7),ncol = 1,nrow = 1)) -> combined.gp
pdf(paste0('TME.umap.',Sys.Date(),'.pdf'),height=8,width = 8.5)
print(combined.gp)
dev.off()


# p<-DimPlot(object = pbmc.f.TME, reduction = "umap",pt.size = 1,
#            group.by = 'cell.type',label = F,label.size = 4, split.by = 'diagnosis2') + 
#   scale_color_manual(values=cell.type.color)
# 
# ggsave('TME.umap.split.by.diagnosis.pdf',p,height=8,width = 35)

###################### Fig 6B TME samples cell composition stackbarplot ######################
pbmc.f.TME<-readRDS('MM.TME.filter.scs.rds')

pbmc.f.TME.use<-subset(pbmc.f.TME,orig.ident%in%names(which(table(pbmc.f.TME$orig.ident)>100)))

cell.type.level2.color<-readRDS('cell.type.level2.color.rds')

data=pbmc.f.TME.use[[c('New.ID','cell.type.level2')]]

data1=as.data.frame(table(data))
data2 = data1 %>% group_by(New.ID) %>% summarise(count=sum(Freq))
data3<-inner_join(data1,data2,by='New.ID')
data3$percentage<- data3$Freq*100/data3$count;

data3$cell.type.level2<-factor(data3$cell.type.level2,levels=c('Mature B cells',
                                                               'CD8T','gdT','MAIT','CD4T','NK',
                                                               'CD14+ Monocytes','CD16+ Monocytes','TAM','cDC1','cDC2','pDC',
                                                               'Platelet',
                                                               'Immature/Progeniors',
                                                               'Endothelial','Fibroblast'))

diagnosis<-unique(pbmc.f.TME.use[[c('New.ID','diagnosis')]])

data3<-left_join(data3,diagnosis,by=c('New.ID'='New.ID'))
data3$diagnosis<-factor(data3$diagnosis,levels=c("normal","MGUS","SMM","NDMM","RRMM"))

Bcell<-subset(data3,cell.type.level2=='Mature B cells')
Bcell<-Bcell%>%group_by(diagnosis)%>%arrange(desc(percentage),.by_group=T)

data3$New.ID<-factor(data3$New.ID,levels=Bcell$New.ID)


plotx<-ggplot(data3, aes(x=New.ID,fill=cell.type.level2,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values = cell.type.level2.color) +
  #scale_x_continuous(breaks=c(0,25,50,75,100)) +
  labs(x='',y = '% (of TME cells)') +
  facet_grid(.~diagnosis,scales = 'free',space = 'free') +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 90,size = 15,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

ggsave('TME.cell.type.level2.composition.stackbarplot.pdf',useDingbats=F,plotx,width = 18,height=8)

###################### Fig 6C TME unsorted samples cell type level2 proportion in TME boxplot ######################
meta<-read.table('MM.filter.meta.data.2021-12-18.txt',sep='\t',header = T)

meta.unsorted<-subset(meta,CD138_sorted=='No' & cell.type.level2!='CD138_pos')

data=meta.unsorted[,c('orig.ident','cell.type.level2')]
data<-subset(data,orig.ident%in%names(which(table(data$orig.ident)>=100)))

data1=as.data.frame(table(data))
data2 = data1 %>% group_by(orig.ident) %>% summarise(count=sum(Freq))
data3<-inner_join(data1,data2,by='orig.ident')
data3$percentage<- data3$Freq*100/data3$count;

data4<-reshape2::dcast(data3,orig.ident~cell.type.level2,value.var='percentage')
rownames(data4)<-data4$orig.ident

data5<-t(scale(data4[,-1]))

data5[data5>2]<-2
data5[data5< -2]<- -2

clinic.info<-readRDS('clinic.info.2.rds')
clinic.info.use<-subset(clinic.info,sample.ID%in%colnames(data5))
rownames(clinic.info.use)<-clinic.info.use$sample.ID
anno.color<-readRDS('CD138_pos.sample_cluster.anno_color.rds')

diagnosis<-unique(meta.unsorted[,c('orig.ident','diagnosis2')])

data6<-left_join(subset(data3,cell.type.level2%in%c('Mature B cells','CD4T','CD8T','MAIT','gdT','NK','CD14+ Monocytes','CD16+ Monocytes','DC','TAM','Immature/Progeniors','Platelet','Endothelial','Fibroblast')),diagnosis,by=c('orig.ident'='orig.ident'))
data6$diagnosis<-data6$diagnosis2
data6$diagnosis<-factor(data6$diagnosis,levels=c("normal","MGUS","SMM","MM"))
data6$cell.type.level2<-factor(data6$cell.type.level2,
                               levels=c('Mature B cells','CD4T','CD8T','MAIT','gdT','NK','CD14+ Monocytes','CD16+ Monocytes','DC','TAM','Immature/Progeniors','Platelet','Endothelial','Fibroblast'))

library('gridExtra')
library('grid')
library('purrr')
library('ggsignif')

signif_plot = function(data,x,y,split_by=NULL,signif.cutoff=0.05, height.factor=0.23,
                       ncol, fill_by=NULL, fill_color='black',
                       color_by=NULL,color='black', ylab) {
  
  # Get full range of y-values
  y_rng = range(data[,y])
  
  # Generate a list of three plots, one for each Species (these are the facets)
  plot_list = lapply(split(data, data[,split_by]), function(d) {
    
    # Get pairs of x-values for current facet
    pairs = combn(levels(data[,x]),2, simplify = F)
    
    # Run wilcox test on every pair
    w.tst =  pairs %>% 
      map_df(function(lv) { 
        p.value = wilcox.test(d$percentage[d[,x]==lv[1]], d$percentage[d[,x]==lv[2]])$p.value
        data.frame(levs=paste(lv, collapse=" "), p.value)
      })
    
    # Record number of significant p.values. We'll use this later to adjust the top of the
    # y-range of the plots
    num_signif = sum(w.tst$p.value <= signif.cutoff)
    
    # Plot significance levels only for combinations with p <= signif.cutoff
    p = ggplot(d, aes_string(x=x, y=y)) +
      geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, aes_string(fill=fill_by), outlier.shape=NA) + 
      geom_jitter(position=position_jitter(width=0.1, height=0), size=3, aes_string(color=color_by)) +
      facet_grid(as.formula(paste("~", split_by)), scales="free", space="free_x") +
      scale_y_continuous(expand = expansion(mult = c(0.15)))+
      geom_signif(test="wilcox.test", comparisons = pairs[which(w.tst$p.value <= signif.cutoff)],
                  map_signif_level = F,            
                  vjust=0.1,
                  textsize=5,
                  size=0.75,
                  step_increase = 0.15) +
      scale_fill_manual(values=fill_color)+
      scale_color_manual(values=color)+
      theme(strip.text.x = element_text(size=14, color="black", face="bold"),
            axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_blank(),
            axis.title=element_blank(),legend.position="none",
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            text = element_text(size=20))
    
    # Return the plot and the number of significant p-values
    return(list(num_signif, p))
  })
  
  # Get the highest number of significant p-values across all three "facets"
  max_signif = max(sapply(plot_list, function(x) x[[1]]))
  
  # Lay out the three plots as facets (one for each Species), but adjust so that y-range is same
  # for each facet. Top of y-range is adjusted using max_signif.
  grid.arrange(grobs=lapply(plot_list, function(x) x[[2]]), 
               ncol=ncol, left="% (of Lymphoid cells)")
}

signif_plot(data6,
            x='diagnosis',
            y='percentage',
            split='cell.type.level2',
            signif.cutoff=0.1,
            ncol=7,
            fill_by='diagnosis',
            fill_color=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','MM'='#d01c8b'),
            ylab='% (of TME cells)')

# wilcox.test(data6$percentage[data6$cell.type.level2=='CD4T' & data6$diagnosis2%in%c('MGUS','SMM')],
#             data6$percentage[data6$cell.type.level2=='CD4T' & data6$diagnosis2=='MM'])

aa<-clinic.info.use[,c('sample.ID','t(4;14)','t(11;14)', 't(14;16)', 't(14;20)', 'Hyperdiploidy','CKS1B/chr1 abnormalities','FISH_1')]
data7<-left_join(data6,aa,by=c('orig.ident'='sample.ID'))
data7$`t(11;14)`<-factor(data7$`t(11;14)`,levels=c('0','1'))
data7$Hyperdiploidy<-factor(data7$Hyperdiploidy,levels=c('0','1'))
data7$`CKS1B/chr1 abnormalities`<-factor(data7$`CKS1B/chr1 abnormalities`,levels=c('0','1'))


dd<-data7[data7$FISH_1%in%c('Hyperdiploidy','Translocation') & data7$diagnosis%in%c('SMM'),]
dd$FISH_1<-factor(dd$FISH_1,levels=c('Translocation', 'Hyperdiploidy'))

signif_plot(dd,
            x='FISH_1',
            y='percentage',
            split='cell.type.level2',
            signif.cutoff=0.2,
            ncol=7,
            color_by = 'diagnosis',
            color=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','MM'='#d01c8b'),
            ylab='% (of TME cells)')


P2<-ggplot(subset(data7,diagnosis%in%c('MGUS','SMM','MM')), aes(x=`t(11;14)`, y=percentage) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=3, aes(color=diagnosis)) +
  #geom_text_repel(data=data7,aes(label=orig.ident),size=2.5) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = combn(levels(data7$`t(11;14)`),2, simplify = F),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
                        vjust=0.1,
                        textsize=5,
                        size=0.75,
                        step_increase = 0.15)+
  scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','MM'='#d01c8b'))+
  facet_wrap(.~cell.type.level2,ncol=7,scales='free') +
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="% (of TME cells)",title='By t(11;14)') 
#ggsave('unsorted.sample.cell.type.level1.composition.t_11_14.boxplot.pdf',useDingbats=F,P,width = 24,height=6)

P3<-ggplot(subset(data7,diagnosis%in%c('MGUS','SMM','MM')), aes(x=Hyperdiploidy, y=percentage) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=3, aes(color=diagnosis)) +
  #geom_text_repel(data=data7,aes(label=orig.ident),size=2.5) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = combn(levels(data7$Hyperdiploidy),2, simplify = F),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
                        vjust=0.1,
                        textsize=5,
                        size=0.75,
                        step_increase = 0.15)+
  scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','MM'='#d01c8b'))+
  facet_wrap(.~cell.type.level2,ncol=7,scales='free') +
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="% (of TME cells)",title='By Hyperdiploidy') 
#ggsave('unsorted.sample.cell.type.level1.composition.Hyperdiploidy.boxplot.pdf',useDingbats=F,P,width = 24,height=6)


P4<-ggplot(subset(data7,diagnosis%in%c('MGUS','SMM','MM')), aes(x=`CKS1B/chr1 abnormalities`, y=percentage) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=3, aes(color=diagnosis)) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = combn(levels(data7$`CKS1B/chr1 abnormalities`),2, simplify = F),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
                        vjust=0.1,
                        textsize=5,
                        size=0.75,
                        step_increase = 0.15)+
  scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','MM'='#d01c8b'))+
  facet_wrap(.~cell.type.level2,ncol=7,scales='free') +
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="% (of TME cells)",title='By CKS1B/chr1 abnormalitie') 
#ggsave('unsorted.sample.cell.type.level1.composition.CKS1B_chr1_abnormalitie.boxplot.pdf',useDingbats=F,P,width = 24,height=6)

P5<-ggplot(data7[data7$FISH_1%in%c('Hyperdiploidy','Translocation') & data7$diagnosis%in%c('MGUS','SMM','MM'),], aes(x=FISH_1, y=percentage) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=3, aes(color=diagnosis)) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = list(c('Hyperdiploidy','Translocation')),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
                        vjust=0.1,
                        textsize=5,
                        size=0.75,
                        step_increase = 0.15)+
  scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','MM'='#d01c8b'))+
  facet_wrap(.~cell.type.level2,ncol=7,scales='free') +
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="% (of TME cells)",title='Hyperdiploidy vs Translocation (MGUS+SMM+MM)') 
#ggsave('unsorted.sample.cell.type.level2.composition.in.TME.primary_events.boxplot.pdf',useDingbats=F,P,width = 24,height=6)

P6<-ggplot(data7[data7$FISH_1%in%c('Hyperdiploidy','Translocation') & data7$diagnosis%in%c('MGUS','SMM'),], aes(x=FISH_1, y=percentage) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=3, aes(color=diagnosis)) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = list(c('Hyperdiploidy','Translocation')),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
                        vjust=0.1,
                        textsize=5,
                        size=0.75,
                        step_increase = 0.15)+
  scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','MM'='#d01c8b'))+
  facet_wrap(.~cell.type.level2,ncol=7,scales='free') +
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="% (of TME cells)",title='Hyperdiploidy vs Translocation (MGUS+SMM)') 

data7$FISH_1<-factor(data7$FISH_1,levels=c('Translocation','Hyperdiploidy','Other','Not available','Normal'))

P7<-ggplot(data7[data7$FISH_1%in%c('Hyperdiploidy','Translocation') & data7$diagnosis%in%c('SMM'),], aes(x=FISH_1, y=percentage) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=3, aes(color=diagnosis)) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = list(c('Hyperdiploidy','Translocation')),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
                        vjust=0.1,
                        textsize=5,
                        size=0.75,
                        step_increase = 0.15)+
  scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','MM'='#d01c8b'))+
  facet_wrap(.~cell.type.level2,ncol=7,scales='free') +
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="% (of TME cells)",title='Hyperdiploidy vs Translocation (SMM)') 

P8<-ggplot(data7[data7$FISH_1%in%c('Hyperdiploidy','Translocation') & data7$diagnosis%in%c('MGUS'),], aes(x=FISH_1, y=percentage) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), size=3, aes(color=diagnosis)) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = list(c('Hyperdiploidy','Translocation')),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
                        vjust=0.1,
                        textsize=5,
                        size=0.75,
                        step_increase = 0.15)+
  scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','MM'='#d01c8b'))+
  facet_wrap(.~cell.type.level2,ncol=7,scales='free') +
  theme(strip.text.x = element_text(size=14, color="black", face="bold"),
        axis.text.x = element_text(angle = 90,size = 14,vjust=0.2,hjust=0.95,colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size=20)) +
  labs(x="", y="% (of TME cells)",title='Hyperdiploidy vs Translocation (MGUS)') 

do.call(ggarrange,c(list(P1,P2,P3,P4,P5,P6,P7,P8),ncol = 1,nrow = 1)) -> combined.gp
pdf('unsorted.sample.cell.type.level2.composition.in.TME.boxplot.pdf',height=12,width = 24)
print(combined.gp)
dev.off()


pdf(paste0('Lymphoid.cell.status.boxplot.',Sys.Date(),'.pdf'),height=12,width=24)
signif_plot(0.05)
dev.off()

###################### Fig 6D CD8T hyperdiploid vs translocation volcano plot ######################
pbmc.f.CD8T<-readRDS('MM.CD8T_gdT_MAIT.filter.scs.rds')

CD8T<-subset(pbmc.f.CD8T,cell.type=='CD8T' & FISH_1%in%c('Translocation','Hyperdiploidy'))
Idents(CD8T)<-CD8T$orig.ident
CD8T<-subset(CD8T,downsample=500)
Idents(CD8T)<-factor(CD8T$FISH_1,levels = c('Hyperdiploidy','Translocation'))

cluster.markers<-FindMarkers(CD8T,min.pct = 0.25,ident.1='Hyperdiploidy',ident.2='Translocation',logfc.threshold = 0)
cluster.markers$symbol<-rownames(cluster.markers)


cluster.markers$log10_p_val_adj= -log10(cluster.markers$p_val_adj)
result=cluster.markers[,c('avg_log2FC','log10_p_val_adj')]
colnames(result)=c('log2FoldChange','log10fdr')
color<-rep('other',dim(result)[1])
color[result$log2FoldChange > 0.3 & result$log2FoldChange <= 0.5 & result$log10fdr > 4]<-'High in Hyperdiploidy';
color[result$log2FoldChange > 0.5 & result$log10fdr > 2]<-'High in Hyperdiploidy';
color[result$log2FoldChange >= -0.5 & result$log2FoldChange < -0.3 & result$log10fdr > 4]<-'High in Translocation';
color[result$log2FoldChange < -0.5 & result$log10fdr > 2]<-'High in Translocation';
result$color<-color;
#gn.selected <- (result$log2FoldChange > 0.5 & result$log10fdr > 2) | (result$log2FoldChange < -0.5 & result$log10fdr > 2) | (result$log2FoldChange < 0.5 & result$log2FoldChange >= 0.3 & result$log10fdr > 4) | (result$log2FoldChange < -0.3 & result$log2FoldChange >= -0.5 & result$log10fdr > 4)
#gn.selected <- result$color!='other'
gn.selected <- rownames(result)%in%c('MTRNR2L8','LTB','RPS4Y1','CD27','CCR7','NOSIP','SELL','TCF7','LEF1','NELL2','SMDT1','TCF7','LDHB',
                                     'NKG7','GZMH','CCL5','GZMB','PRF1','FGFBP2','KLRD1','LGALS1','CCL5','GZMA','AL138963.4','MTRNR2L12','PLEK','GNLY','LGALS1',
                                     'GZMA','FCGR3A','CTSW','CST7')

result$log2FoldChange[result$log2FoldChange< -1]<- -1
plotx<-ggplot(result,aes(x=log2FoldChange, y=log10fdr)) +
  geom_point(aes(fill=color,size=log10fdr),color='black',shape=21) +
  scale_size(range = c(0,4), limits = c(0, 150), breaks = c(0,10,100,150)) +
  scale_fill_manual(values = c('red','blue','grey')) +
  ggtitle(paste0('Hyperdiploidy',' vs Translocation in ','CD8T')) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  #ylim(0,40) +
  xlim(-1,1)+
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = subset(result, gn.selected),aes(label = rownames(subset(result, gn.selected))),size = 3,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0('CD8T_hyperdiploid_vs_translocation.Volcano.plot.',Sys.Date(),'.pdf'),plotx,width = 6,height = 6.5)

write.table(cluster.markers,'CD8T_hyperdiploid_vs_translocation.deg.txt',sep='\t',quote=F,col.names = T,row.names = F)

###################### Fig 6E CD8 UMAP ######################
pbmc.f.CD8T<-readRDS('MM.CD8T.filter.scs.rds')

p1<-DimPlot(object = pbmc.f.CD8T, reduction = "umap",pt.size = 0.5,group.by = 'orig.ident',label = F,label.size = 4)+ guides(color=guide_legend(ncol=2,override.aes = list(size=2)))
p1<-AugmentPlot(p1)
p2<-DimPlot(object = pbmc.f.CD8T, reduction = "umap",pt.size = 0.5,group.by = 'orig.ident',label = F,label.size = 4)+ guides(color=guide_legend(ncol=1,override.aes = list(size=2)))

p3<-DimPlot(object = pbmc.f.CD8T, reduction = "umap",pt.size = 1,group.by = 'cell.status',label = F,label.size = 10)+ guides(color=guide_legend(ncol=2,override.aes = list(size=2)))
p3<-AugmentPlot(p3)
p4<-DimPlot(object = pbmc.f.CD8T, reduction = "umap",pt.size = 0.5,group.by = 'cell.status',label = T,label.size = 4)+ guides(color=guide_legend(ncol=1,override.aes = list(size=2)))

p5<-DimPlot(object = pbmc.f.CD8T, reduction = "umap",pt.size = 0.5,group.by = 'diagnosis',label = F,label.size = 4)+ guides(color=guide_legend(ncol=2,override.aes = list(size=2)))
p5<-AugmentPlot(p5)
p6<-DimPlot(object = pbmc.f.CD8T, reduction = "umap",pt.size = 0.5,group.by = 'diagnosis',label = F,label.size = 4)+ guides(color=guide_legend(ncol=1,override.aes = list(size=2)))

p7<-DimPlot(object = pbmc.f.CD8T, reduction = "umap",pt.size = 0.5,group.by = 'time_point',label = F,label.size = 4)+ guides(color=guide_legend(ncol=2,override.aes = list(size=2)))
p7<-AugmentPlot(p7)
p8<-DimPlot(object = pbmc.f.CD8T, reduction = "umap",pt.size = 0.5,group.by = 'time_point',label = F,label.size = 4)+ guides(color=guide_legend(ncol=1,override.aes = list(size=2)))

do.call(ggarrange,c(list(p1,p2,p3,p4,p5,p6,p7,p8),ncol = 1,nrow = 1)) -> combined.gp
pdf(paste0('CD8T.umap.',Sys.Date(),'.pdf'),height=8,width = 8.5)
print(combined.gp)
dev.off()

###################### Fig 6F CD8 functional score ######################
pbmc.f.CD8T<-readRDS('MM.CD8T.filter.scs.rds')

gp.list<-NULL
for(i in c('TRM','predysfunctional','dysfunctional','naive.like','cytotoxic')){
  
  p<-FeaturePlot(object = pbmc.f.CD8T, reduction = "umap",features=i,pt.size = 1.5,label = F,label.size = 4,max.cutoff = 1.5,raster = T,order=F) +
    scale_color_gradientn(colors=pals::jet()[c(1:11,15:25)],na.value = "#e6e6e6")
  gp.list[[length(gp.list)+1]]<-p
}

do.call(ggpubr::ggarrange,c(gp.list,ncol = 1,nrow = 1)) -> combined.gp
pdf(paste0('CD8T.functional.score.umap.',Sys.Date(),'.pdf'),height=8,width = 9)
print(combined.gp)
dev.off()
###################### Fig 6G CD8 naive like vs cytotoxic density plot ######################
pbmc.f.CD8T<-readRDS('MM.CD8T.filter.scs.rds')

library('pals')
gdat<-FetchData(pbmc.f.CD8T,vars=c('naive.like','cytotoxic','diagnosis2','FISH_1'))
gdat<-as.data.frame(gdat);
gdat$diagnosis2<-factor(gdat$diagnosis2,levels=c('normal','MGUS','SMM','MM'))
gdat$FISH_1<-factor(gdat$FISH_1,levels=c('Normal','Hyperdiploidy','Translocation','Other','Not available'))

g.overlay = ggplot(data = gdat,
                   aes(x = naive.like, y = cytotoxic)) #+ xlim(-20,20) + ylim(-20,20)
g.overlay = g.overlay + stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)
#g.overlay = g.overlay + geom_point(color = 'white',size = .0001)
g.overlay = g.overlay + facet_wrap(~diagnosis2,ncol = 4)
g.overlay = g.overlay + scale_fill_gradientn(colours = turbo()[c(5:23)],limits=c(0,10),oob = scales::squish)
g.overlay = g.overlay + theme_bw()
ggsave(paste0('CD8T.diagnosis.density.plot.',Sys.Date(),'.pdf'),height = 4,width=15,g.overlay)

g.overlay = ggplot(data = gdat,
                   aes(x = naive.like, y = cytotoxic)) #+ xlim(-20,20) + ylim(-20,20)
g.overlay = g.overlay + stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)
#g.overlay = g.overlay + geom_point(color = 'white',size = .0001)
g.overlay = g.overlay + facet_wrap(~FISH_1,ncol = 5)
g.overlay = g.overlay + scale_fill_gradientn(colours = turbo()[c(5:23)],limits=c(0,10),oob = scales::squish)
g.overlay = g.overlay + theme_bw()
ggsave(paste0('CD8T.FISH.density.plot.',Sys.Date(),'.pdf'),height = 4,width=19,g.overlay)
