#--------------------------------------------------------------
# filename : Figure3.R
# Date : 2023-09-01
# contributor : Minghao Dang, PhD
# function: 
# R version: R/4.2.0
#--------------------------------------------------------------

source('global_settings.R')

###################### Fig 3A CD138+ sample cluster level driver genes heatmap order by primary driver events ######################
pbmc.f.sdc1<-readRDS('MM.CD138_pos.filter.scs.rds')

gene<-c('CKS1B','CTHRC1','LIME1','SDF2L1','LY6E','DUSP1','CD99','PTPRCAP',
        'SLAMF1','FAM30A','ITM2C','LAMP5','SULF2',
        'FGFR3','NSD2','CCND1','MAF','MAFB','CCND2',
        'MYC','CDKN2C',
        'CD38','SDC1','TNFRSF17','GPRC5D','FCRL5',
        'ITGB7', 'SPP1', 'FRZB','CADM1','LAPTM5','WFDC2','TXNIP','KLF2','IRF4','ASS1','CD47',
        'HLA-A','HLA-B','HLA-C')
p<-DotPlot(pbmc.f.sdc1,features = gene,group.by = 'New.Cluster')
data<-p$data[,c('id','features.plot','avg.exp.scaled')]
colnames(data)=c('cluster','gene','avg.exp.scaled')

cluster_orig<-unique(pbmc.f.sdc1[[c('New.ID','New.Cluster')]])
data1<-left_join(data,cluster_orig,by=c('cluster'='New.Cluster'))

data1$New.ID<-factor(data1$New.ID,levels=c('nBM01s','MGUS04s','MGUS03s','SMM06s','SMM29','NDMM05s',
                                           'MGUS08','MGUS05s','SMM33','SMM03s','SMM19','SMM12',
                                           'SMM01s','SMM17','SMM23','SMM02s','SMM28','NDMM06',
                                           'NDMM01s','MGUS01s','MGUS09','MGUS13','MGUS10',
                                           'MGUS07','MGUS02s','SMM11','SMM27','SMM08','SMM05s',
                                           'SMM07s','SMM25','SMM20','SMM16','SMM22','SMM21',
                                           'SMM04s','SMM18','SMM15','SMM24','NDMM02s','NDMM04s',
                                           'NDMM07','MGUS06s','SMM26','SMM30','NDMM03s','MGUS11','SMM13'))

data1<-data1[with(data1,order(New.ID,cluster)),]

data2=reshape2::dcast(data1,gene~cluster,value.var="avg.exp.scaled")
rownames(data2)=data2$gene

data3<- as.matrix(data2[,-1])
data3<-data3[,unique(data1$cluster)]
data3[which(data3>1.5)]=1.5
data3[which(data3< -1.5)]=-1.5

anno.row<-readRDS('CD138_pos.sample_cluster.anno.rds')
#anno.row$top1.clonotype.proportion<-data4[rownames(anno.row),'clonotype.proportion']

anno.color<-readRDS('CD138_pos.sample_cluster.anno_color.rds')


orig.ident.color=rep(c("lightgrey", "darkgrey"), 24)
names(orig.ident.color)<-levels(data1$New.ID)
anno.color[['New.ID']]=orig.ident.color 

data4_order<-readRDS('CD138_pos.sample_cluster.order.rds')

data4<-data3[,data4_order$New.Cluster]

gap<-readRDS('CD138_pos.sample_cluster.order.gap.rds')


aa<-unique(subset(pbmc.f.sdc1,type=='Malignant')[[c('cluster','FISH_2')]])
aa$FISH_2<-factor(aa$FISH_2,levels=c('Normal','t(4;14)','t(11;14)','t(14;16)','t(14;20)','Hyperdiploidy','Other','Not available'))
gap<-cumsum(c(ncol(data4_1),table(aa$FISH_2)))

plotx<-pheatmap::pheatmap(data4, fontsize_col = 10,
                          color=colorRampPalette(c('blue','black','yellow'))(100),
                          border_color = NA,
                          treeheight_row =15,
                          treeheight_col =30,
                          cluster_rows = F,
                          show_rownames =T,
                          show_colnames = T,
                          display_numbers=F,
                          #legend_breaks = c(0,10,20,30,40,50,60),
                          #legend_labels = c('0%','10%','20%','30%','40%','50%','60%'),
                          cluster_cols = F,
                          gaps_col = gap,
                          annotation_col = anno.row[,c('New.ID','Diagnosis','cnv.score','top1.clonotype.proportion','cluster2','FISH_2')],
                          annotation_colors = append(anno.color,list('cluster2'=c('C0'='#543005',
                                                                                  'C1'='#8c510a',
                                                                                  'C2'='#bf812d',
                                                                                  'C3'='#dfc27d',
                                                                                  'C4'='#f6e8c3',
                                                                                  'C5'='#c7eae5',
                                                                                  'C6'='#80cdc1',
                                                                                  'C7'='#35978f',
                                                                                  'C8'='#01665e',
                                                                                  'C9'='#003c30')))
);
pdf('CD138_pos.sample.cluster.level.driver.gene.heatmap.order.by.normal_malignant.primary.driver.event.pdf',width=20,height=10)
print(plotx)
dev.off()

###################### Fig 3B correlation of genes with CytoTRACE by cytogenetics (heatmap) ######################

pbmc.f.sdc1<-readRDS('MM.CD138_pos.filter.scs.rds')

set.seed(12345)
pbmc.f.sdc1.sub<-subset(pbmc.f.sdc1,downsample=1000)


gene<-c('SLAMF1','ITM2C','LAMP5','SULF2',
        'CKS1B','CTHRC1','LIME1','SDF2L1','LY6E','DUSP1','CD99',#'PTPRCAP',
        'NSD2','CCND1','MAF','MAFB','CCND2',#'FGFR3',
        'MYC','CDKN2C',
        'CD38','SDC1','TNFRSF17','GPRC5D','FCRL5',
        'ITGB7', 'SPP1', 'FRZB','CADM1','LAPTM5','WFDC2','TXNIP','KLF2','IRF4','ASS1','CD47',
        'HLA-A','HLA-B','HLA-C')

gene.exp<-t(GetAssayData(pbmc.f.sdc1.sub,slot = 'data'))
aa<-read.table('./CytoTRACE/CytoTRACE_plot_table.txt',sep='\t',header = T,row.names = 1)

aa.sub<-aa[rownames(aa)%in%colnames(pbmc.f.sdc1.sub),]

aa.use<-subset(aa.sub,Phenotype=='Normal')
gene.exp.use<-gene.exp[rownames(aa.use),]
bb<-cor(as.matrix(gene.exp.use),aa.use[,c('CytoTRACE')],method = 'spearman')
colnames(bb)<-'Normal_cor_with_CytoTRACE'


aa.use<-subset(aa.sub,Phenotype=='Translocation')
gene.exp.use<-gene.exp[rownames(aa.use),]
bb_2<-cor(as.matrix(gene.exp.use),aa.use[,c('CytoTRACE')],method = 'spearman')
colnames(bb_2)<-'Translocation_cor_with_CytoTRACE'


aa.use<-subset(aa.sub,Phenotype=='Hyperdiploidy')
gene.exp.use<-gene.exp[rownames(aa.use),]
#gene.exp.use<-gene.exp.use[,apply(gene.exp.use,2,function(x){(sum(x>0)/length(x))>0.1})]
bb_3<-cor(as.matrix(gene.exp.use),aa.use[,c('CytoTRACE')],method = 'spearman')
colnames(bb_3)<-'HPD_cor_with_CytoTRACE'


cor_mtx<-merge(bb,bb_2,by='row.names',all=T)
rownames(cor_mtx)<-cor_mtx$Row.names
cor_mtx$Row.names<-NULL

cor_mtx<-merge(cor_mtx,bb_3,by='row.names',all=T)
rownames(cor_mtx)<-cor_mtx$Row.names
cor_mtx$Row.names<-NULL

cor_mtx.output<-data.frame('gene'=rownames(cor_mtx),cor_mtx)
write.table(cor_mtx.output,'gene.expr.spearman_cor.with.CytoTRACE.by.cytogenetics.downsampling.txt',row.names = F,col.names = T,sep='\t')
# cor_mtx<-cor_mtx[apply(cor_mtx,1,function(x){sum(is.na(x))<3}),]
# 
# cor_mtx[cor_mtx>0.5]<-0.5
# cor_mtx[cor_mtx< -0.5]<- -0.5
# cor_mtx[is.na(cor_mtx)]<-0
# 
# cor_mtx<-cor_mtx[apply(cor_mtx,1,function(x){sum(abs(x)>0.2)>0}),]
# cor_mtx[abs(cor_mtx)<0.2]<-0

cor_mtx.use<-cor_mtx[gene,]
cor_mtx.use[cor_mtx.use>0.5]<-0.5

col_fun<-circlize::colorRamp2(c(-0.5, 0, max(cor_mtx.use,na.rm=T)), c("#2166ac", "white", "#b2182b"))
lgd = Legend(col_fun = col_fun, title = "spearman cor with CytoTRACE", at = round(c(-0.5, 0, max(cor_mtx.use,na.rm=T)),1))

pdf('gene.rank.by.cor.with.CytoTRACE.downsampling.pdf',height=10,width=5)
ht_list<-Heatmap(cor_mtx.use, #name = "spearman cor with CytoTRACE",
                 col=col_fun,
                 cluster_rows = T,
                 cluster_columns = F,
                 show_row_names=T,
                 #row_names_side = "bottom",
                 show_column_names = T,
                 column_names_side = "bottom",
                 #column_split=aa.anno$cluster,
                 column_gap = unit(1, "mm"),
                 #right_annotation = gene_ann_v,
                 show_heatmap_legend = F) 
draw(ht_list, annotation_legend_list = lgd, annotation_legend_side = "right")
dev.off()

###################### Fig 3C correlation of genes with CytoTRACE by cytogenetics (smoothed curve)  ######################
pbmc.f.sdc1<-readRDS('MM.CD138_pos.filter.scs.rds')

set.seed(12345)
pbmc.f.sdc1.sub<-subset(pbmc.f.sdc1,downsample=1000)
#pbmc.f.sdc1.sub<-subset(pbmc.f.sdc1.sub,Phase=='G1')

gene<-c('SLAMF1','ITM2C','LAMP5','SULF2',
        'CKS1B','CTHRC1','LIME1','SDF2L1','LY6E','DUSP1','CD99',#'PTPRCAP',
        'NSD2','CCND1','MAF','MAFB','CCND2',#'FGFR3',
        'MYC','CDKN2C',
        'CD38','SDC1','TNFRSF17','GPRC5D','FCRL5',
        'ITGB7', 'SPP1', 'FRZB','CADM1','LAPTM5','WFDC2','TXNIP','KLF2','IRF4','ASS1','CD47',
        'HLA-A','HLA-B','HLA-C','ICAM1','ICAM3','ICAM4','CD74')

gene.exp<-FetchData(pbmc.f.sdc1.sub,vars = gene)
diagnosis<-FetchData(pbmc.f.sdc1.sub,vars = 'diagnosis')
aa<-read.table('./CytoTRACE/CytoTRACE_plot_table.txt',sep='\t',header = T,row.names = 1)
aa.sub<-aa[rownames(aa)%in%colnames(pbmc.f.sdc1.sub),]


bb<-cbind(aa.sub[,c('Phenotype','CytoTRACE')],diagnosis[rownames(aa.sub),],gene.exp[rownames(aa.sub),])
colnames(bb)[3]<-'diagnosis'
bb<-subset(bb,Phenotype%in%c('Hyperdiploidy','Normal','Translocation'))

cc<-reshape2::melt(bb,id.vars=c('Phenotype','CytoTRACE','diagnosis'))

scaleFUN <- function(x) sprintf("%.1f", x)
p1 <- ggplot(cc,aes(x = CytoTRACE, y = value, color=Phenotype)) +
  scale_color_manual(values=c('Normal'='#00a651','Translocation'='#8dd3c7',
                              'Hyperdiploidy'='#ff7f00'))+
  geom_smooth() +
  ylab('expression')+
  xlab('CytoTRACE')+
  facet_wrap(.~variable,ncol=10,scales='free')+
  scale_x_continuous(labels=scaleFUN,breaks = c(0,0.5,1))+
  #theme_classic()+
  theme(#axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.title = element_blank(),
    axis.line = element_blank(),#axis.ticks=element_blank(),
    #strip.text.x = element_text(size=12, color="black", face="bold"),
    strip.background = element_blank(),strip.text = element_text(angle = 0),
    #strip.text.y = element_text(angle = 0),
    text = element_text(size=12),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black",fill=NA, size=0.5))

ggsave(paste0('gene.exp.vs.CytoTRACE.downsampling.',Sys.Date(),'.pdf'),p1,height=8,width=15)
