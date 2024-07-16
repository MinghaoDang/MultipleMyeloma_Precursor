#--------------------------------------------------------------
# filename : Figure1.R
# Date : 2023-09-01
# contributor : Minghao Dang, PhD
# function: 
# R version: R/4.2.0
#--------------------------------------------------------------

###################### Fig 2A malignant cell sample pseudobulk phylogenetic tree ######################
sdc1.malignant<-readRDS('MM.CD138_pos.malignant.filter.scs.rds')
Idents(sdc1.malignant)<-factor(sdc1.malignant$New.ID)

aa=AverageExpression(sdc1.malignant)
aa=aa$RNA

vg <- FindVariableFeatures(aa)
vg.rm.IgV<-vg[!grepl('^IG.*V',rownames(vg)),]

vg.rm.IgV<-rownames(vg.rm.IgV[order(vg.rm.IgV$vst.variance.standardized,decreasing = T),])[1:1000]

aa1<-aa[vg.rm.IgV,]
aa1<-log2(aa1+1)

bb=cor(aa1,method ='spearman')
cc=hclust(dist(bb))

library('ggtree')
library('ape')
clinic.info<-readRDS('clinic.info.2.rds')
group=clinic.info[,c('New.ID','FISH_1','Diagnosis')]
group<-group[rownames(bb),]

tree=as.phylo(cc)
node_order=fortify(tree)
node_order = subset(node_order, isTip)
node_order<-node_order$label[order(node_order$y, decreasing=F)]

pp <- ggtree(tree,size=1,branch.length='none') +
  geom_tiplab(size = 0) +
  xlim(-5, NA) +
  layout_dendrogram()

pp <- pp %<+% group + 
  geom_label(aes(label=Diagnosis,fill=Diagnosis),size=0) +
  scale_fill_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b'))

metadata<-sdc1.malignant@meta.data

d1<-unique(metadata[,c('New.ID','Translocation','Hyperdiploidy','chr1_1q_amp','Hypodiploidy')])
d1$Hyperdiploidy[!is.na(d1$Translocation)]<-NA
rownames(d1)<-d1$New.ID
d2<-t(d1[node_order,-1])

d3<-reshape2::melt(t(d2))
d3$Var1<-factor(d3$Var1,levels=node_order)
d3$Var2<-factor(d3$Var2,levels=c('Hypodiploidy','chr1_1q_amp','Hyperdiploidy','Translocation'))

plotx<-ggplot(data = d3, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color='black') +
  scale_fill_manual(values=c('t(4;14)'='#8dd3c7',
                             't(11;14)'='#bebada',
                             't(14;16)'='#80b1d3',
                             't(14;20)'='#fccde5',
                             'Hyperdiploidy'='#ff7f00',
                             'chr1_1q_amp'='#878787',
                             'Hypodiploidy'='#878787'), na.value = "white")+
  theme(axis.text.x = element_text(angle = 90,color='black',vjust=0.5,hjust=1),
        axis.text.y = element_text(color='black'),
        axis.title = element_blank(),
        axis.line = element_blank(),axis.ticks=element_blank(),
        strip.background = element_blank(),strip.text = element_text(angle = 0),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

pdf('malignant.sample.pseudobulk.cytogenetics.heatmap.pdf',width=15,height=5)
plot_grid(pp,plotx,align = 'v',ncol=1)
dev.off()


###################### Fig 2B Bhat distance ######################
library('fpc')
sdc1.malignant<-readRDS('MM.CD138_pos.malignant.filter.scs.rds')
Idents(sdc1.malignant)<-factor(sdc1.malignant$New.ID)

Embeddings<-Embeddings(sdc1.malignant)[,1:50]
sample.id<-levels(sdc1.malignant)

bd.mtx=matrix(0,nrow=length(sample.id),ncol=length(sample.id))

for(i in 1:length(sample.id)){
  for(j in 1:length(sample.id)){
    
    if(i==j){bd.mtx[i,j]=0}else{
      
      cells1=colnames(sdc1.malignant)[sdc1.malignant$New.ID==sample.id[i]]
      dat1<-Embeddings[cells1,]
      cells2=colnames(sdc1.malignant)[sdc1.malignant$New.ID==sample.id[j]]
      dat2<-Embeddings[cells2,]
      mu1<-apply(dat1,2,mean);
      Sigma1<-cov(dat1);
      mu2<-apply(dat2,2,mean);
      Sigma2<-cov(dat2);
      bd<-bhattacharyya.dist(mu1, mu2, Sigma1, Sigma2)
      
      bd.mtx[i,j]=bd
    }
  }
  
}

colnames(bd.mtx)=sample.id
rownames(bd.mtx)=sample.id

saveRDS(bd.mtx,'bd.mtx.rds')


bd.mtx[upper.tri(bd.mtx,diag = T)]<-NA

bd.mtx.pair<-reshape2::melt(bd.mtx)
colnames(bd.mtx.pair)<-c('s1','s2','b_dist')
bd.mtx.pair<-bd.mtx.pair[!is.na(bd.mtx.pair$b_dist),]

clinic.info<-readRDS('clinic.info.2.rds')
group=clinic.info[,c('New.ID','FISH_1','Diagnosis')]

bd.mtx.pair<-left_join(bd.mtx.pair,group,by=c('s1'='New.ID'))
colnames(bd.mtx.pair)[c(4,5)]<-c('s1_FISH_1','s1_Diagnosis')

bd.mtx.pair<-left_join(bd.mtx.pair,group,by=c('s2'='New.ID'))
colnames(bd.mtx.pair)[c(6,7)]<-c('s2_FISH_1','s2_Diagnosis')


bd.mtx.pair$s1_FISH_1[bd.mtx.pair$s1_FISH_1=='Hyperdiploidy']<-'HPD'
bd.mtx.pair$s1_FISH_1[bd.mtx.pair$s1_FISH_1=='Translocation']<-'Trans'
bd.mtx.pair$s2_FISH_1[bd.mtx.pair$s2_FISH_1=='Hyperdiploidy']<-'HPD'
bd.mtx.pair$s2_FISH_1[bd.mtx.pair$s2_FISH_1=='Translocation']<-'Trans'

bd.mtx.pair$diagnosis_group<-paste0(bd.mtx.pair$s1_Diagnosis,'_',bd.mtx.pair$s2_Diagnosis)
bd.mtx.pair$diagnosis_group[bd.mtx.pair$diagnosis_group=='SMM_MGUS']<-'MGUS_SMM'
bd.mtx.pair$diagnosis_group[bd.mtx.pair$diagnosis_group=='NDMM_MGUS']<-'MGUS_NDMM'
bd.mtx.pair$diagnosis_group<-factor(bd.mtx.pair$diagnosis_group,levels=c('MGUS_MGUS','MGUS_SMM','MGUS_NDMM','SMM_SMM','SMM_NDMM','NDMM_NDMM'))

bd.mtx.pair$FISH_group<-paste0(bd.mtx.pair$s1_FISH_1,'_',bd.mtx.pair$s2_FISH_1)
bd.mtx.pair.use<-bd.mtx.pair[bd.mtx.pair$FISH_group%in%c('HPD_HPD','HPD_Trans','Trans_HPD','Trans_Trans'),]
bd.mtx.pair.use$FISH_group[bd.mtx.pair.use$FISH_group=='Trans_HPD']<-'HPD_Trans'
bd.mtx.pair.use$FISH_group<-factor(bd.mtx.pair.use$FISH_group,levels=c('HPD_HPD','HPD_Trans','Trans_Trans'))

bd.mtx.pair.use$total_group<-paste0(bd.mtx.pair.use$s1_Diagnosis,'_',bd.mtx.pair.use$s1_FISH_1,'-',bd.mtx.pair.use$s2_Diagnosis,'_',bd.mtx.pair.use$s2_FISH_1)
bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='MGUS_Trans-MGUS_HPD']<-'MGUS_HPD-MGUS_Trans'
bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='NDMM_HPD-MGUS_HPD']<-'MGUS_HPD-NDMM_HPD'
bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='NDMM_HPD-MGUS_Trans']<-'MGUS_Trans-NDMM_HPD'
bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='NDMM_Trans-MGUS_HPD']<-'MGUS_HPD-NDMM_Trans'
bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='NDMM_Trans-MGUS_Trans']<-'MGUS_Trans-NDMM_Trans'
bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='SMM_HPD-MGUS_HPD']<-'MGUS_HPD-SMM_HPD'
bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='NDMM_Trans-NDMM_HPD']<-'NDMM_HPD-NDMM_Trans'
bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='SMM_HPD-MGUS_Trans']<-'MGUS_Trans-SMM_HPD'
bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='SMM_Trans-MGUS_HPD']<-'MGUS_HPD-SMM_Trans'

bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='SMM_Trans-MGUS_Trans']<-'MGUS_Trans-SMM_Trans'
bd.mtx.pair.use$total_group[bd.mtx.pair.use$total_group=='SMM_Trans-SMM_HPD']<-'SMM_HPD-SMM_Trans'

bd.mtx.pair.use$total_group<-factor(bd.mtx.pair.use$total_group,levels=c("MGUS_HPD-MGUS_HPD","MGUS_Trans-MGUS_Trans", 
                                                                         "SMM_HPD-SMM_HPD","SMM_Trans-SMM_Trans",
                                                                         "NDMM_HPD-NDMM_HPD","NDMM_Trans-NDMM_Trans",
                                                                         "MGUS_Trans-NDMM_Trans","MGUS_Trans-SMM_Trans","SMM_Trans-NDMM_Trans",
                                                                         "MGUS_HPD-MGUS_Trans","SMM_HPD-SMM_Trans","NDMM_HPD-NDMM_Trans",
                                                                         "MGUS_HPD-SMM_HPD","MGUS_HPD-NDMM_HPD","SMM_HPD-NDMM_HPD",
                                                                         "MGUS_Trans-SMM_HPD" ,
                                                                         "MGUS_Trans-NDMM_HPD",
                                                                         "MGUS_HPD-SMM_Trans",
                                                                         "MGUS_HPD-NDMM_Trans",
                                                                         "SMM_Trans-NDMM_HPD",
                                                                         "SMM_HPD-NDMM_Trans"))

bd.mtx.pair.diagnosis_group.median<-bd.mtx.pair %>%
  group_by(diagnosis_group) %>%
  summarise(median = median(b_dist, na.rm=TRUE)) %>%
  arrange(median)
bd.mtx.pair$diagnosis_group<-factor(bd.mtx.pair$diagnosis_group,levels=bd.mtx.pair.diagnosis_group.median$diagnosis_group)

p1<-ggplot(bd.mtx.pair, aes(x=diagnosis_group, y=b_dist) ) + 
  geom_boxplot(alpha = 0.4, show.legend = T, outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), alpha=1) +
  theme_boxplot() +
  labs(x="", y="Bhattacharyya distance")

bd.mtx.pair.use$FISH_group<-factor(bd.mtx.pair.use$FISH_group,levels=c('Trans_Trans','HPD_HPD','HPD_Trans'))
p2<-ggplot(bd.mtx.pair.use, aes(x=FISH_group, y=b_dist) ) + 
  geom_boxplot(alpha = 0.4, show.legend = T, outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), alpha=1) +
  geom_signif(test="wilcox.test",
              comparisons = list(c('HPD_HPD', 'HPD_Trans'),c('HPD_HPD','Trans_Trans'),c('HPD_Trans', 'Trans_Trans')),
              map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)
              vjust=0.5,
              textsize=4,
              size=0.5,
              step_increase = 0.15)+
  theme_boxplot() +
  labs(x="", y="Bhattacharyya distance")

p3<-ggplot(bd.mtx.pair.use, aes(x=total_group, y=b_dist) ) + 
  geom_boxplot(alpha = 0.4, show.legend = T, outlier.shape=NA) + 
  geom_jitter(position=position_jitter(width=0.1, height=0), alpha=0.9) +
  geom_hline(yintercept=mean(bd.mtx.pair.use$b_dist), linetype="dashed", color = "blue")+
  theme_boxplot() +
  labs(x="", y="Bhattacharyya distance")


bd.mtx.pair.use.group.median<-bd.mtx.pair.use %>%
  group_by(total_group) %>%
  summarise(median = median(b_dist, na.rm=TRUE)) %>%
  arrange(median)
bd.mtx.pair.use$total_group<-factor(bd.mtx.pair.use$total_group,levels=bd.mtx.pair.use.group.median$total_group)

bd.mtx.pair.use$total_group2<-'same_diagnosis_same_genetics'
bd.mtx.pair.use$total_group2[bd.mtx.pair.use$s1_FISH_1==bd.mtx.pair.use$s2_FISH_1 & bd.mtx.pair.use$s1_Diagnosis!=bd.mtx.pair.use$s2_Diagnosis]<-'different_diagnosis_same_genetics'
bd.mtx.pair.use$total_group2[bd.mtx.pair.use$s1_FISH_1!=bd.mtx.pair.use$s2_FISH_1 & bd.mtx.pair.use$s1_Diagnosis==bd.mtx.pair.use$s2_Diagnosis]<-'same_diagnosis_different_genetics'
bd.mtx.pair.use$total_group2[bd.mtx.pair.use$s1_FISH_1!=bd.mtx.pair.use$s2_FISH_1 & bd.mtx.pair.use$s1_Diagnosis!=bd.mtx.pair.use$s2_Diagnosis]<-'different_diagnosis_different_genetics'

bd.mtx.pair.use$diagnosis_group<-factor(bd.mtx.pair.use$diagnosis_group,levels=c('MGUS_MGUS',
                                                                                 'SMM_SMM',
                                                                                 'NDMM_NDMM',
                                                                                 'MGUS_SMM',
                                                                                 'MGUS_NDMM',
                                                                                 'SMM_NDMM'))

pval_MGUS_MGUS<-pairwise.wilcox.test(bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='MGUS_MGUS',]$b_dist,bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='MGUS_MGUS',]$total_group,pool.sd = FALSE,p.adjust.method='none')$p.value
pval_SMM_SMM<-pairwise.wilcox.test(bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='SMM_SMM',]$b_dist,bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='SMM_SMM',]$total_group,pool.sd = FALSE,p.adjust.method='none')$p.value
pval_NDMM_NDMM<-pairwise.wilcox.test(bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='NDMM_NDMM',]$b_dist,bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='NDMM_NDMM',]$total_group,pool.sd = FALSE,p.adjust.method='none')$p.value
pval_MGUS_SMM<-pairwise.wilcox.test(bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='MGUS_SMM',]$b_dist,bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='MGUS_SMM',]$total_group,pool.sd = FALSE,p.adjust.method='none')$p.value
pval_MGUS_NDMM<-pairwise.wilcox.test(bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='MGUS_NDMM',]$b_dist,bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='MGUS_NDMM',]$total_group,pool.sd = FALSE,p.adjust.method='none')$p.value
pval_SMM_NDMM<-pairwise.wilcox.test(bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='SMM_NDMM',]$b_dist,bd.mtx.pair.use[bd.mtx.pair.use$diagnosis_group=='SMM_NDMM',]$total_group,pool.sd = FALSE,p.adjust.method='none')$p.value

p4<-ggplot(bd.mtx.pair.use, aes(x=total_group, y=b_dist) ) + 
  geom_boxplot(alpha = 1, show.legend = T, outlier.shape=NA, fill='#b8e186') + 
  geom_jitter(position=position_jitter(width=0.1, height=0)) +
  geom_tile(data=bd.mtx.pair.use, aes(x = total_group, y = 0, fill = s1_FISH_1), width=0.8, color='black')+
  geom_tile(data=bd.mtx.pair.use, aes(x = total_group, y = -1, fill = s2_FISH_1), width=0.8, color='black')+
  geom_tile(data=bd.mtx.pair.use, aes(x = total_group, y = -2, fill = s1_Diagnosis), width=0.8, color='black')+
  geom_tile(data=bd.mtx.pair.use, aes(x = total_group, y = -3, fill = s2_Diagnosis), width=0.8, color='black')+
  scale_fill_manual(values=c('Trans'='#8dd3c7','HPD'='#ff7f00',
                             'MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b'))+
  facet_grid(.~diagnosis_group,scales = 'free',space='free') +
  stat_compare_means(comparisons=list(c(1,2), c(1,3), c(2,3)), method="wilcox") +
  
  theme_boxplot() +
  labs(x="", y="Bhattacharyya distance")

plotx<-ggarrange(ggarrange(p1,p2, ncol = 2, labels = c("A", "B")),                                                 # First row with scatter plot
                 p4, # Second row with box and dot plots
                 nrow = 2                                      # Labels of the scatter plot
)
ggsave('bhat_distance.pdf',plotx,height=10,width=15)


###################### Fig 2C correlate malignant cell ITH score driver events and diagnosis ######################
data<-as.data.frame(read_excel('malignant.properties.xlsx',sheet = 1))
data$Diagnosis<-factor(data$Diagnosis,levels=c('MGUS','SMM','NDMM'))

P1<-ggplot(subset(data,Diagnosis%in%c('MGUS','SMM','NDMM') & primary_driver!='N'), aes(x=primary_driver, y=ITH.PC.DistED) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, outlier.shape=NA) + # aes(fill=Hipodiploidy),
  geom_jitter(position=position_jitter(width=0.1, height=0), size=3,shape=16, aes(color=Diagnosis)) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = list(c('H','T')),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
                        vjust=0.1,
                        textsize=5,
                        size=0.75,
                        step_increase = 0.15)+
  scale_color_manual(values=c('normal'='#4dac26','MGUS'='#1c75bc','SMM'='#fbb040','NDMM'='#d01c8b'))+
  #facet_wrap(.~cell.status,ncol=12,scales='free') +
  theme_boxplot() +
  labs(x="", y="ITH score")

ggsave('malignant.ITH_score.driver_event.correlation.pdf',useDingbats=F,P1,width = 14,height=8)


###################### Fig 2D correlate malignant cell ITH score with 1q gain ######################
data<-as.data.frame(read_excel('malignant.properties.xlsx',sheet = 1))
data.use<-subset(data,primary_driver!='N')

data.use$group<-'non-HY'
data.use$group[data.use$New.ID%in%c('NDMM05s','NDMM06','SMM29','MGUS13','MGUS10')]<-'non-HY-1q gain' # ,'MGUS09','MGUS01s'
data.use$group[data.use$New.ID%in%c('MGUS02s','MGUS21','MGUS12','MGUS07','SMM18','SMM07s','SMM08','SMM05s','SMM16','SMM24','SMM20','SMM27','SMM15','SMM04s','SMM11','NDMM07')]<-'HY'
data.use$group[data.use$New.ID%in%c('NDMM02s','NDMM04s','SMM25','SMM21','SMM22')]<-'HY-1q gain'

data.use$group<-factor(data.use$group,levels=c('non-HY','HY','non-HY-1q gain','HY-1q gain'))
data.use$Diagnosis<-factor(data.use$Diagnosis,levels=c('MGUS','SMM','NDMM'))

data.use$group2<-'1q_gain_Pos'
data.use$group2[data.use$group%in%c('non-HY','HY')]<-'1q_gain_Neg'

P1<-ggplot(data.use, aes(x=group, y=ITH.PC.DistED) ) + 
  geom_boxplot(alpha = 0.8, show.legend = T,width = 0.75, outlier.shape=NA) + # aes(fill=Hipodiploidy),
  geom_jitter(position=position_jitter(width=0.1, height=0), size=3, aes(shape=Diagnosis,color=FISH_2)) +
  ggsignif::geom_signif(test="wilcox.test", comparisons = list(c('non-HY','non-HY-1q gain'),c('HY','HY-1q gain'),c('non-HY','HY'),c('non-HY-1q gain','HY-1q gain')),
                        map_signif_level = F,   # c("***"=0.001, "**"=0.01, "*"=0.05)       
                        vjust=0.1,
                        textsize=5,
                        size=0.75,
                        step_increase = 0.15)+
  scale_color_manual(values=c('t(4;14)'='#8dd3c7','t(11;14)'='#bebada','t(14;16)'='#80b1d3','t(14;20)'='#fccde5','Hyperdiploidy'='#ff7f00'))+
  #facet_wrap(.~cell.status,ncol=12,scales='free') +
  scale_y_continuous(expand = expansion(mult = c(0.15)))+
  theme_boxplot() +
  labs(x="", y="ITH score")

ggsave('malignant.ITH_score.1qgain.correlation.pdf',useDingbats=F,P1,width = 14,height=8)


###################### Fig 2E CD138+ samples subcluster composition stackbarplot by cytogenetics ######################
pbmc.f.sdc1.malignant<-readRDS('MM.CD138_pos.malignant.filter.scs.rds')

pbmc.f.sdc1.malignant.use<-subset(pbmc.f.sdc1.malignant,FISH_1%in%c('Translocation','Hyperdiploidy'))

data=pbmc.f.sdc1.malignant.use[[c('New.ID','New.Cluster')]]
data$New.Cluster<-gsub('(.*)_(.*)','\\2',data$New.Cluster)


data1=as.data.frame(table(data))
data2 = data1 %>% group_by(New.ID) %>% summarise(count=sum(Freq))
data3<-inner_join(data1,data2,by='New.ID')
data3$percentage<- data3$Freq*100/data3$count;

fish<-unique(pbmc.f.sdc1.malignant.use[[c('New.ID','FISH_2')]])

data3<-left_join(data3,fish,by=c('New.ID'='New.ID'))
data3$FISH_2<-factor(data3$FISH_2,levels=c('t(4;14)','t(11;14)','t(14;16)','t(14;20)',"Hyperdiploidy"))

data4<-subset(data1,Freq>0)
data3$New.ID<-factor(data3$New.ID,levels=names(sort(table(data4$New.ID))))
data3$New.Cluster<-factor(data3$New.Cluster,levels=paste0('C',c(9:0)))

plotx<-ggplot(data3, aes(x=New.ID,fill=New.Cluster,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values = c('C0'='#543005',
                               'C1'='#8c510a',
                               'C2'='#bf812d',
                               'C3'='#dfc27d',
                               'C4'='#f6e8c3',
                               'C5'='#c7eae5',
                               'C6'='#80cdc1',
                               'C7'='#35978f',
                               'C8'='#01665e',
                               'C9'='#003c30')) +
  #scale_x_continuous(breaks=c(0,25,50,75,100)) +
  labs(x='',y = '% (of CD138+ cells)') +
  facet_grid(.~FISH_2,scales = 'free',space = 'free') +
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

ggsave('CD138_pos.intrapatient.composition.stackbarplot.by.cytogenetics.pdf',useDingbats=F,plotx,width = 14,height=8)

###################### Fig 2E CD138+ samples BCR clonotype composition stackbarplot by cytogenetics ######################
pbmc.f.sdc1.malignant<-readRDS('MM.CD138_pos.malignant.filter.scs.rds')

pbmc.f.sdc1.malignant.use<-subset(pbmc.f.sdc1.malignant,FISH_1%in%c('Translocation','Hyperdiploidy'))

data=pbmc.f.sdc1.malignant.use[[c('New.ID','BCR_clonotype_raw')]]
data$BCR_clonotype_raw<-gsub('(.*)_(.*)','\\1',data$BCR_clonotype_raw)

data$BCR_clonotype_raw[data$BCR_clonotype_raw%in%paste0('clonotype',1)]<-'Top1'
data$BCR_clonotype_raw[data$BCR_clonotype_raw%in%paste0('clonotype',2:10)]<-'Top2-10'
data$BCR_clonotype_raw[data$BCR_clonotype_raw%in%paste0('clonotype',11:100)]<-'Top11-100'
data$BCR_clonotype_raw[data$BCR_clonotype_raw%in%paste0('clonotype',101:100000)]<-'>Top100'

data1=as.data.frame(table(data))
data2 = data1 %>% group_by(New.ID) %>% summarise(count=sum(Freq))
data3<-inner_join(data1,data2,by='New.ID')
data3$percentage<- data3$Freq*100/data3$count;

fish<-unique(pbmc.f.sdc1.malignant.use[[c('New.ID','FISH_2')]])

data3<-left_join(data3,fish,by=c('New.ID'='New.ID'))
data3$FISH_2<-factor(data3$FISH_2,levels=c('t(4;14)','t(11;14)','t(14;16)','t(14;20)',"Hyperdiploidy"))

data3$New.ID<-factor(data3$New.ID,levels=names(sort(table(data4$New.ID)))) # data4 see above
data3$BCR_clonotype_raw<-factor(data3$BCR_clonotype_raw,levels=rev(c('Top1','Top2-10','Top11-100','>Top100')))



plotx<-ggplot(data3, aes(x=New.ID,fill=BCR_clonotype_raw,y=percentage)) +
  geom_bar(stat="identity",position = 'stack') +
  scale_fill_manual(values = c('Top1'='#543005',
                               'Top2-10'='#bf812d',
                               'Top11-100'='#f6e8c3',
                               '>Top100'='#80cdc1')) +
  #scale_x_continuous(breaks=c(0,25,50,75,100)) +
  labs(x='',y = '% (of CD138+ cells)') +
  facet_grid(.~FISH_2,scales = 'free',space = 'free') +
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

ggsave('CD138_pos.sample.BCR.clonotype.cytogenetics.stackbarplot.pdf',useDingbats=F,plotx,width = 14,height=8)

