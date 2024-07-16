#--------------------------------------------------------------
# filename : Figure1.R
# Date : 2023-09-01
# contributor : Minghao Dang, PhD
# function: 
# R version: R/4.2.0
#--------------------------------------------------------------

source('global_settings.R')

###################### Fig 1A patient pie chart ######################

data1<-data.frame('diagnosis'=c('nBM','MGUS','SMM','NDMM','RRMM'),
                  'No'=c(4,21,32,7,1))

data1$diagnosis<-factor(data1$diagnosis,levels=c('nBM','MGUS','SMM','NDMM','RRMM'))

plotx<-ggplot(data1, aes(x = '', y = No, fill = diagnosis)) +
  geom_bar(width = 0.01, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(label = diagnosis), color = "black",position = position_stack(vjust = 0.5))+
  scale_fill_manual(values=c('nBM'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7')) +
  #facet_wrap(~k4,ncol = 4)+
  theme_void()

ggsave('patient.piechart.pdf',useDingbats=F,plotx,width = 6,height=6)

###################### Fig 1B draw clinical info figure group by diagnosis ######################
library('readxl')
library('dplyr')

data<-read_excel('../SCS\ data\ 9.20.21\ DE-IDENTIFIED\ ZB\ finalized\ with\ new\ ID.xlsx',sheet=1,skip=1)

fish<-data[,c('New ID','FISH BM aspirate positive abnormalities')]

aa<-strsplit(fish$`FISH BM aspirate positive abnormalities`,split=',')
names(aa)<-fish$`New ID`

bb<-as.data.frame(plyr::ldply(aa, rbind))

cc<-reshape2::melt(bb,id.vars = '.id')
cc<-cc[!is.na(cc$value),c('.id','value')]
#cc<-cc[cc$value!='NA',]
colnames(cc)<-c('New.ID','FISH')

dd<-reshape2::dcast(cc,New.ID~FISH)
dd[is.na(dd)]<-0
rownames(dd)<-dd$New.ID
dd<-dd[,-1]
dd[dd!=0]<-1

ee<-apply(dd,2,as.numeric)
rownames(ee)<-rownames(dd)

clinic.info<-readRDS('clinic.info.2.rds')
rownames(clinic.info)<-clinic.info$New.ID
clinic.color<-readRDS('clinic.color.rds')


gaps_col<-cumsum(table(clinic.info$Diagnosis))

ff<-t(ee[clinic.info$New.ID,])

ff<-ff[c('Normal',
         't(4;14)','t(11;14)','t(14;16)','t(14;20)',
         'trisomy/tetrasomy 4','trisomy/tetrasomy 9','trisomy/tetrasomy 11','trisomy/tetrasomy 13','trisomy/tetrasomy 14','tetrasomy 16','trisomy/tetrasomy 17',
         'trisomy/tetrasomy 8/MYC rearrangement',
         'trisomy 1','gain(1q)','del(1p)','monosomy 4','monosomy 13','monosomy 14','monosomy 16','monosomy 17','NA'),]

p<-pheatmap::pheatmap(ff[,clinic.info$New.ID],
                      cluster_cols = F,
                      cluster_rows = F,
                      show_colnames = T,
                      show_rownames = T,
                      border_color = 'grey',
                      color = c('white','black'),
                      gaps_col=gaps_col,
                      annotation_col = clinic.info[,c('Diagnosis','CD138 sorted')],
                      annotation_colors = clinic.color)

ggsave('FISH.pdf',p,width=20,height=10)


meta<-read.table('MM.filter.meta.data.2022-12-15.txt',sep='\t',header = T)

data=meta[,c('New.ID','cell.type')]
data$cell.type[data$cell.type%in%c('Malignant_PC','Normal_PC')]<-'CD138_pos'
data$cell.type[data$cell.type!='CD138_pos']<-'TME'

data1=as.data.frame(table(data))
data1$New.ID<-factor(data1$New.ID,levels=clinic.info$New.ID)
data1$cell.type<-factor(data1$cell.type,levels=c('TME','CD138_pos'))

data1$diagnosis<-clinic.info[data1$New.ID,'Diagnosis']

plotx<-ggplot(data1, aes(x=New.ID,fill=cell.type,y=Freq)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values = c('CD138_pos'='#525252','TME'='#5bbcc2')) +
  #scale_y_continuous(breaks=c(0,25,50,75,100)) +
  #labs(x='',y = '% cells') +
  facet_grid(.~diagnosis,space='free',scales = 'free') +
  geom_hline(yintercept=50, linetype="dashed", color = "red") +
  #coord_trans(y="log10") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(angle = 90,size = 15,vjust=0.2,hjust=0.95),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5,size=20))

ggsave('cell.count.per.sample.pdf',plotx,width=20,height=10)

###################### Fig 1C overall UMAP ######################
# pbmc.f<-readRDS('MM.filter.scs.rds')
# pbmc.f<-DietSeurat(pbmc.f,counts=F,dimreducs=c('umap'))
pbmc.f<-readRDS('MM.filter.slim.scs.rds')

p1<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'orig.ident',label = F,label.size = 4)+theme(legend.position = "none")
p2<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'cell.type.level1',label = F,label.size = 4)+ 
  scale_color_manual(values=cell.type.level1.color)+
  theme(legend.position = "none")

p3<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'diagnosis',label = F,label.size = 4) +
  scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'))+
  theme(legend.position = "none")

p4<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'FISH_2',label = F,label.size = 4)+ 
  scale_color_manual(values=c('Normal'='#00a651','t(4;14)'='#8dd3c7','t(11;14)'='#bebada','t(14;16)'='#80b1d3','t(14;20)'='#fccde5','Hyperdiploidy'='#ff7f00','Other'='#7f3f98','NA'='#878787'))+
  theme(legend.position = "none")

p5<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'FISH_3',label = F,label.size = 4)+ 
  scale_color_manual(values=c('Normal'='#00a651','chr1_1q_amp'='#8dd3c7','t(11;14)'='#bebada','t(14;16)'='#80b1d3','Hyperdiploidy'='#ff7f00','Hypodiploidy'='#7f3f98','Not available'='#878787'))+
  theme(legend.position = "none")

p6<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'cell.type.level2',label = F,label.size = 4)+ 
  scale_color_manual(values=cell.type.level2.color)+
  theme(legend.position = "none")

p7<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'cell.type.level2',label = F,label.size = 4,raster = F)+ 
  scale_color_manual(values=cell.type.level2.color)+
  theme(legend.position = "none")

p8<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'Translocation',label = F,label.size = 4,order = T)+ 
  scale_color_manual(values=c('t(4;14)'='#8dd3c7','t(11;14)'='#bebada','t(14;16)'='#80b1d3','t(14;20)'='#fccde5','NA'='#878787'))+
  theme(legend.position = "none")

p9<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'chr1_1q_amp',label = F,label.size = 4,order = F)+ 
  scale_color_manual(values=c('chr1_1q_amp'='#ff7f00','NA'='#878787'))+
  theme(legend.position = "none")

p10<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'Hyperdiploidy',label = F,label.size = 4,order = F)+ 
  scale_color_manual(values=c('Hyperdiploidy'='#ff7f00','NA'='#878787'))+
  theme(legend.position = "none")

p11<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'Hypodiploidy',label = F,label.size = 4,order = F)+ 
  scale_color_manual(values=c('Hypodiploidy'='#ff7f00','NA'='#878787'))+
  theme(legend.position = "none")


library(randomcoloR)
set.seed(5555)
ucols=c(distinctColorPalette(1360),distinctColorPalette(1360),distinctColorPalette(1361))
p12<-DimPlot(object = pbmc.f, reduction = "umap",pt.size = 0.5,group.by = 'BCR_clonotype',label = F,label.size = 4,raster = T, order=F)+ 
  theme(legend.position = "none")


do.call(ggarrange,c(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12),ncol = 1,nrow = 1)) -> combined.gp
pdf(paste0('MM.overall.umap.',Sys.Date(),'.pdf'),height=8,width = 8.5)
print(combined.gp)
dev.off()


gene.check=unique(c('SDC1','CD19','CD79A','CD79B','MS4A1','CD3D','CD3E','CD4','CD8A','CD14','FCGR3A','CD68','CD163','HBB','HBA1','HBA2','CD34','IGLL1','FCGR1A','CD1C','MPO','LILRA4','PPBP','GP1BB','PF4','MKI67'))
batch_feature_plot(gene.check,pbmc.f,paste0('featureplot.',Sys.Date(),'.pdf'))

###################### Fig 1D CD138+% boxplot by diagnosis ######################
pbmc.f<-readRDS('MM.filter.scs.rds')

data=pbmc.f[[c('New.ID','cell.type.level2','CD138_sorted')]]
data<-subset(data,CD138_sorted=='No')
data$CD138_sorted<-NULL

data$cell.type.level2[data$cell.type.level2%in%c('Normal','Tumor')]<-'CD138_pos'

data1=as.data.frame(table(data))
data2 = data1 %>% group_by(New.ID) %>% summarise(count=sum(Freq))
data3<-inner_join(data1,data2,by='New.ID')
data3$percentage<- data3$Freq*100/data3$count;

diagnosis<-unique(pbmc.f[[c('New.ID','diagnosis')]])

data3<-left_join(data3,diagnosis,by=c('New.ID'='New.ID'))
data3$diagnosis<-factor(data3$diagnosis,levels=c("normal","MGUS","SMM","NDMM","RRMM"))

data3_unsorted<-data3

cell.count<-unique(data3_unsorted[,c('New.ID','count')])
data3_unsorted<-subset(data3_unsorted,New.ID%in%cell.count$New.ID[cell.count$count>50])

data4<-data3_unsorted[data3_unsorted$cell.type.level2=='CD138_pos',]
p1<-ggplot(data4, aes(x=factor(diagnosis), y=percentage) ) + 
  geom_boxplot(alpha = 1, show.legend = T, aes(fill=factor(diagnosis)), outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), alpha=0.9,shape=16) +
  scale_fill_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b', 'RRMM'='#8856a7'))+
  geom_signif(test="wilcox.test",
              comparisons = list(c('normal','MGUS'),c('MGUS','SMM'),c('SMM','NDMM'),c('NDMM','RRMM')),
              map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
              vjust=0.5,
              textsize=4,
              size=0.5,
              step_increase = 0.06)+
  theme_boxplot()
ggsave('CD138_pos.pct.by.diagnosis.pdf',useDingbats=F,p1,width = 5,height=5)

###################### Fig 1E & Fig S3A CD138+ UMAP ######################
pbmc.f.sdc1<-readRDS('MM.CD138_pos.filter.scs.rds')

p1<-DimPlot(object = pbmc.f.sdc1, reduction = "umap",pt.size = 1,group.by = 'orig.ident',label = F,label.size = 4,raster = T)+ 
  theme(legend.position = "none")

p2<-DimPlot(object = pbmc.f.sdc1, reduction = "umap",pt.size = 1,group.by = 'diagnosis',label = F,label.size = 4,raster = T)+ 
  scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b'))+
  theme(legend.position = "none")

pbmc.f.sdc1$type<-factor(pbmc.f.sdc1$type,levels=c('Normal','Malignant'))
p3<-DimPlot(object = pbmc.f.sdc1, reduction = "umap",pt.size = 1,group.by = 'type',label = F,label.size = 4,raster = T)+ 
  scale_color_manual(values=c('Normal'='#4dac26','Malignant'='#525252'))+
  theme(legend.position = "none")

p4<-DimPlot(object = pbmc.f.sdc1, reduction = "umap",pt.size = 1,group.by = 'FISH_1',label = F,label.size = 4,raster = T)+ 
  scale_color_manual(values=c('Normal'='#00a651','Translocation'='#fbb040','Hyperdiploidy'='#be1e2d','Other'='#7f3f98','NA'='#878787'))+
  theme(legend.position = "none")

p5<-DimPlot(object = pbmc.f.sdc1, reduction = "umap",pt.size = 1,group.by = 'FISH_2',label = F,label.size = 4,raster = T)+ 
  scale_color_manual(values=c('Normal'='#00a651','t(4;14)'='#8dd3c7','t(11;14)'='#bebada','t(14;16)'='#80b1d3','t(14;20)'='#fccde5','Hyperdiploidy'='#ff7f00','Other'='#7f3f98','NA'='#878787'))+
  theme(legend.position = "none")


library(randomcoloR)
set.seed(5555)
ucols=c(distinctColorPalette(1360),distinctColorPalette(1360),distinctColorPalette(1361))
p6<-DimPlot(object = pbmc.f.sdc1, reduction = "umap",pt.size = 1,group.by = 'BCR_clonotype_raw',label = F,label.size = 4,raster = T,cols=ucols, order=T)+ 
  theme(legend.position = "none")

p7<-FeaturePlot(object = pbmc.f.sdc1,
                pt.size = 1,
                features = 'infercnv.score.sd.normal_CD138_ref',#max.cutoff = 3,
                cols = c("white", "#2CA25F"),
                max.cutoff=0.075,
                reduction = "umap",
                order=F,
                combine = T,
                raster = T
) + scale_color_gradientn(colors=pals::jet()[c(1:11,15:25)],na.value = "#00007F") + theme(legend.position = "none")

do.call(ggarrange,c(list(p1,p2,p3,p4,p5,p6,p7),ncol = 1,nrow = 1)) -> combined.gp
pdf(paste0('CD138_pos.umap.',Sys.Date(),'.pdf'),height=8,width = 8.5)
print(combined.gp)
dev.off()

###################### Fig 1F CD138+ malignant sample level driver genes violin plot group by cytogenetics ######################
pbmc.f.sdc1<-readRDS('MM.CD138_pos.filter.scs.rds')
pbmc.f.sdc1.sub<-subset(pbmc.f.sdc1,type=='Malignant' | diagnosis=='normal')
data<-FetchData(pbmc.f.sdc1.sub,vars=c('New.ID','diagnosis','FISH_2','CCND1','CCND2','MAF','MAFB','CCND3','FGFR3','NSD2','MYC'))

aa<-reshape2::melt(data)
aa$FISH<-factor(aa$FISH_2,levels=c('Normal','t(4;14)','t(11;14)','t(14;16)','t(14;20)','Hyperdiploidy','Other','Not available'))
aa$diagnosis<-factor(aa$diagnosis,levels=c('normal','MGUS','SMM','NDMM'))
aa$variable<-factor(aa$variable,levels=c('CCND1','CCND2','MAF','MAFB','CCND3','FGFR3','NSD2','MYC'))
aa<-aa%>%arrange(FISH,diagnosis,variable)
aa$New.ID<-factor(aa$New.ID,levels=c('nBM01s','MGUS03s','SMM06s',
                                     'SMM29','NDMM05s',
                                     'MGUS08','MGUS05s','SMM33','SMM03s','SMM19','SMM12','SMM01s','SMM17','SMM23','SMM02s','SMM28','NDMM06','NDMM01s',
                                     'MGUS01s','MGUS09','MGUS13',
                                     'MGUS10',
                                     'MGUS07','MGUS02s','SMM11','SMM27','SMM08','SMM05s','SMM07s','SMM25','SMM20','SMM16','SMM22','SMM21','SMM04s','SMM18','SMM15','SMM24','NDMM02s','NDMM04s','NDMM07',
                                     'MGUS06s','SMM26','NDMM03s',
                                     'SMM13'))

P<-ggplot(aa,aes(x=New.ID,y=value,fill=diagnosis)) +
  geom_violin(scale='width') + # adjust=1.5,bw=0.3,scale='area',color='white'
  stat_summary(fun=median,geom='point',
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b'))+
  facet_grid(rows=vars(variable),cols = vars(FISH),scales = 'free',space = 'free_x') +
  theme_vlnplot() +
  labs(x='', y="Expression")

ggsave('malignant.translocation.gene.expression.pdf',useDingbats=F,P,width=20,height=8)
