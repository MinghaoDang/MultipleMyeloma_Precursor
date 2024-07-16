#--------------------------------------------------------------
# filename : Figure7.R
# Date : 2023-09-01
# contributor : Minghao Dang, PhD
# function: 
# R version: R/4.2.0
#--------------------------------------------------------------

source('global_settings.R')


###################### Fig 7A CellPhoneDB count heatmap ######################

hyperdiploid_count<-read.table('./CellPhoneDB/out/hyperdiploid/count_network.txt',header = T,sep = '\t')
hyperdiploid_count$SOURCE<-factor(hyperdiploid_count$SOURCE,levels=c('CD4T','CD8T','gdT','MAIT','NK','Mature B cells','CD14+ Monocytes','CD16+ Monocytes','DC','TAM','Platelet','Malignant_PC'))
hyperdiploid_count$TARGET<-factor(hyperdiploid_count$TARGET,levels=c('CD4T','CD8T','gdT','MAIT','NK','Mature B cells','CD14+ Monocytes','CD16+ Monocytes','DC','TAM','Platelet','Malignant_PC'))

translocation_count<-read.table('./CellPhoneDB/out/translocation/count_network.txt',header = T,sep = '\t')
translocation_count$SOURCE<-factor(translocation_count$SOURCE,levels=c('CD4T','CD8T','gdT','MAIT','NK','Mature B cells','CD14+ Monocytes','CD16+ Monocytes','DC','TAM','Platelet','Malignant_PC'))
translocation_count$TARGET<-factor(translocation_count$TARGET,levels=c('CD4T','CD8T','gdT','MAIT','NK','Mature B cells','CD14+ Monocytes','CD16+ Monocytes','DC','TAM','Platelet','Malignant_PC'))

T_H<-left_join(hyperdiploid_count,translocation_count,by=c('SOURCE'='SOURCE','TARGET'='TARGET'))
T_H$`H-T`<-T_H$count.x-T_H$count.y
plox<-ggplot(T_H, aes(SOURCE, TARGET, fill= `H-T`)) + 
  geom_tile() +
  scale_fill_gradient2(low = scales::muted("blue",l=30,c=100),
                       mid = "white",
                       high = scales::muted("red",l=20,c=100))+
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        #strip.text.x = element_text(size=12, color="black", face="bold"),
        strip.background = element_blank(),strip.text = element_text(angle = 0),
        #strip.text.y = element_text(angle = 0),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())
ggsave('hyperdiploid_translocation.diff.CellPhoneDB.count.pdf',plox,width=6.5,height=5.5)
###################### Fig 7B CellPhoneDB mean exp, p value ######################
hyperdiploid_mean<-read.table('./CellPhoneDB/out/hyperdiploid/means.txt',header = T,sep = '\t', check.names = F)
hyperdiploid_mean<-reshape2::melt(hyperdiploid_mean[,c(2,12:ncol(hyperdiploid_mean))])
colnames(hyperdiploid_mean)<-c('interacting_pair','Source|Target','mean')
hyperdiploid_p<-read.table('./CellPhoneDB_v4/out/hyperdiploid/pvalues.txt',header = T,sep = '\t', check.names = F)
hyperdiploid_p<-reshape2::melt(hyperdiploid_p[,c(2,12:ncol(hyperdiploid_p))])
colnames(hyperdiploid_p)<-c('interacting_pair','Source|Target','p.value')
hyperdiploid_mean_p<-left_join(hyperdiploid_mean,hyperdiploid_p,by=c('interacting_pair'='interacting_pair','Source|Target'='Source|Target'))
hyperdiploid_mean_p$diagnosis<-'hyperdiploid'

translocation_mean<-read.table('./CellPhoneDB/out/translocation/means.txt',header = T,sep = '\t', check.names = F)
translocation_mean<-reshape2::melt(translocation_mean[,c(2,12:ncol(translocation_mean))])
colnames(translocation_mean)<-c('interacting_pair','Source|Target','mean')
translocation_p<-read.table('./CellPhoneDB/out/translocation/pvalues.txt',header = T,sep = '\t', check.names = F)
translocation_p<-reshape2::melt(translocation_p[,c(2,12:ncol(translocation_p))])
colnames(translocation_p)<-c('interacting_pair','Source|Target','p.value')
translocation_mean_p<-left_join(translocation_mean,translocation_p,by=c('interacting_pair'='interacting_pair','Source|Target'='Source|Target'))
translocation_mean_p$diagnosis<-'translocation'

data<-rbind(hyperdiploid_mean_p,translocation_mean_p)
data$diagnosis<-factor(data$diagnosis,levels=c('hyperdiploid','translocation'))
data$log10pvalue<- -log10(data$p.value)
data$log10pvalue[data$log10pvalue>3]<-3
data$log2mean<- log2(data$mean+1)


# tumor -> TME
data1<-subset(data,grepl('Malignant_PC\\|',data$`Source|Target`))
data1<-subset(data1,!grepl('\\|Malignant_PC',data1$`Source|Target`))

sel1<-data1 %>% group_by(interacting_pair) %>% 
  top_n(n=1, wt = log10pvalue) %>% 
  filter(log10pvalue> -log10(0.05)) %>% 
  pull(interacting_pair) %>% unique() %>%
  as.character()
data4<-subset(data1,interacting_pair%in%sel1)

sel2<-data4 %>% group_by(`Source|Target`) %>% 
  top_n(n=1, wt = log10pvalue) %>% 
  filter(log10pvalue> -log10(0.01)) %>% 
  pull(`Source|Target`) %>% unique() %>%
  as.character()
data4<-subset(data4,`Source|Target`%in%sel2)

sel3<-data4[,c('interacting_pair','Source|Target','diagnosis','log10pvalue')] %>% group_by(interacting_pair,`Source|Target`) %>% 
  mutate(p_diff=log10pvalue[diagnosis=='hyperdiploid']-log10pvalue[diagnosis=='translocation']) %>% 
  group_by(interacting_pair)%>%summarise(n=sum(p_diff))%>%
  filter(abs(n)>1) %>% 
  arrange(desc(n)) %>% 
  pull(interacting_pair) %>% unique() %>%
  as.character()
saveRDS(sel3,'tumor_TME.L-R.pairs.rds')
data4<-subset(data4,interacting_pair%in%sel3)

data4$`Source|Target`<-factor(data4$`Source|Target`,levels=c('Malignant_PC|CD4T','Malignant_PC|CD8T','Malignant_PC|gdT','Malignant_PC|MAIT','Malignant_PC|NK','Malignant_PC|Mature B cells','Malignant_PC|CD14+ Monocytes','Malignant_PC|CD16+ Monocytes','Malignant_PC|DC','Malignant_PC|TAM','Malignant_PC|Platelet','Malignant_PC|Malignant_PC'))
data4$interacting_pair<-factor(data4$interacting_pair,levels=rev(sel3))

plotx<-ggplot(data4, aes(x = diagnosis,y = interacting_pair)) +        ## global aes
  geom_point(aes(fill = log2mean,size =log10pvalue),color='black',shape=21)  +    ## geom_point for circle illusion
  #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
  scale_fill_viridis()+
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  scale_size(range = c(0,6), limits = c(0, 3), breaks = c(0,1,2,3)) +             ## to tune the size of circles
  facet_grid(cols=vars(`Source|Target`),scales='free',space='free')+
  labs(x='',y='')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), #panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text = element_text(angle = 90,vjust=0.5,hjust=0.05),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95),
        legend.position="right")
ggsave(paste0('CellPhoneDB.tumor_TME.interaction.',Sys.Date(),'.pdf'),plotx,height=11,width=7.5,limitsize = FALSE)

# TME -> tumor
data1<-subset(data,grepl('\\|Malignant_PC',data$`Source|Target`))
data1<-subset(data1,!grepl('Malignant_PC\\|',data1$`Source|Target`))
sel1<-data1 %>% group_by(interacting_pair) %>% 
  top_n(n=1, wt = log10pvalue) %>% 
  filter(log10pvalue> -log10(0.01)) %>% 
  pull(interacting_pair) %>% unique() %>%
  as.character()
data4<-subset(data1,interacting_pair%in%sel1)

sel2<-data4 %>% group_by(`Source|Target`) %>% 
  top_n(n=1, wt = log10pvalue) %>% 
  filter(log10pvalue> -log10(0.01)) %>% 
  pull(`Source|Target`) %>% unique() %>%
  as.character()
data4<-subset(data4,`Source|Target`%in%sel2)

sel3<-data4[,c('interacting_pair','Source|Target','diagnosis','log10pvalue')] %>% group_by(interacting_pair,`Source|Target`) %>% 
  mutate(p_diff=log10pvalue[diagnosis=='hyperdiploid']-log10pvalue[diagnosis=='translocation']) %>% 
  group_by(interacting_pair)%>%summarise(n=sum(p_diff))%>%
  filter(abs(n)>0.5) %>% 
  arrange(desc(n)) %>% 
  pull(interacting_pair) %>% unique() %>%
  as.character()
saveRDS(sel3,'TME_tumor.L-R.pairs.rds')
data4<-subset(data4,interacting_pair%in%sel3)

data4$`Source|Target`<-factor(data4$`Source|Target`,levels=c('CD4T|Malignant_PC','CD8T|Malignant_PC','gdT|Malignant_PC','MAIT|Malignant_PC','NK|Malignant_PC','Mature B cells|Malignant_PC','CD14+ Monocytes|Malignant_PC','CD16+ Monocytes|Malignant_PC','DC|Malignant_PC','TAM|Malignant_PC','Platelet|Malignant_PC','Malignant_PC|Malignant_PC'))
data4$interacting_pair<-factor(data4$interacting_pair,levels=rev(sel3))
plotx<-ggplot(data4, aes(x = diagnosis,y = interacting_pair)) +        ## global aes
  geom_point(aes(fill = log2mean,size =log10pvalue),color='black',shape=21)  +    ## geom_point for circle illusion
  #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
  scale_fill_viridis()+
  scale_fill_gradientn(colours=c('#f7fcfd','#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')) +
  scale_size(range = c(0,6), limits = c(0, 3), breaks = c(0,1,2,3)) +             ## to tune the size of circles
  facet_grid(cols=vars(`Source|Target`),scales='free',space='free')+
  labs(x='',y='')+
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2), #panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text = element_text(angle = 90,vjust=0.5,hjust=0.05),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95),
        legend.position="right")
ggsave(paste0('CellPhoneDB.TME_tumor.interaction.',Sys.Date(),'.pdf'),plotx,height=11,width=8,limitsize = FALSE)

###################### Fig 7C expression of ICAM1 ######################
pbmc.f.sdc1.malignant<-readRDS('MM.CD138_pos.malignant.filter.scs.rds')

gene<-c('ICAM1')

p<-DotPlot(pbmc.f.sdc1.malignant,features = gene,group.by = 'orig.ident')
data<-p$data[,c('id','features.plot','avg.exp','pct.exp')]
colnames(data)=c('orig.ident','gene','avg.exp','pct.exp')

anno.row<-readRDS('CD138_pos.sample_cluster.anno.rds')
anno.color<-readRDS('CD138_pos.sample_cluster.anno_color.rds')

data<-left_join(data,unique(anno.row[,c('orig.ident','Diagnosis','FISH_1')]),by=c('orig.ident'='orig.ident'))

data<-subset(data,Diagnosis!='normal bone marrow')

data$avg.exp<-log2(data$avg.exp+1)
scaleFUN <- function(x) sprintf("%.2f", x)


P<-ggplot(subset(data,FISH_1%in%c('Hyperdiploidy','Translocation')), aes(x=FISH_1, y=avg.exp) ) + 
  geom_boxplot(alpha = 0.6, show.legend = T,  outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0),aes(color=Diagnosis)) +
  # geom_text_repel(data = subset(data,orig.ident%in%c('NDMM-07','SMM-06')),aes(label = orig.ident),size = 3,
  #                 box.padding = unit(0.35, "lines"),
  #                 point.padding = unit(0.3, "lines"))+
  ggsignif::geom_signif(test="wilcox.test", comparisons = list(c('Hyperdiploidy','Translocation')),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)
                        vjust=0,
                        textsize=4,
                        size=0.6,
                        step_increase = 0.1) +
  facet_wrap(vars(gene),ncol=5,scales='free') +
  scale_color_manual(values=anno.color$Diagnosis)+
  theme_boxplot() +
  labs(x="", y="Average Expression") +
  scale_y_continuous(labels=scaleFUN,expand = expansion(mult = c(0.15)))
ggsave('ICAM1.malignant.gene.expr.H.vs.T.compare.boxplot.pdf',P,width=10,height=5)
