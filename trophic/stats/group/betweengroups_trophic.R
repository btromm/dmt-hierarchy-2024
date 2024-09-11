pacman::p_load(rstatix,tidyverse,scales,reshape2,data.table,cowplot,ggpubr,lme4,merTools,lmerTest,emmeans,ggradar,ragg,svglite)

# Directedness
dat<-read.csv('/Users/myco/human_brains/LSDAyaCannab/bigdirectedness.csv',header=FALSE)
condition<-read.csv('/Users/myco/human_brains/LSDAyaCannab/bigdirectedness_condition.csv',header=FALSE)
condition2<-read.csv('/Users/myco/human_brains/LSDAyaCannab/condition2bigtc.csv',header=FALSE)
dir2<-data.frame(as.vector(dat),as.vector(condition),as.vector(condition2))
colnames(dir2)[1]='Directedness'
colnames(dir2)[2]='Condition'
colnames(dir2)[3]='Condition2'
dir2$Condition<-factor(dir2$Condition,levels=c("Baseline","Ayahuasca","Placebo","DMT","Chronic/Sober","Chronic/Cannabis","Occasional/Sober","Occasional/Cannabis"))
dir2$Condition2<-factor(dir2$Condition2,levels=c("Ayahuasca","DMT","Chronic","Occasional"))

outlier <- dir2 %>% group_by(Condition) %>% identify_outliers(Directedness)
data.frame(outlier)

dir2 <- dir2[-which(dir2$Directedness %in% outlier$Directedness[outlier$is.extreme=='TRUE']),]

dir2 %>% ggplot(aes(Condition,Directedness,fill=Condition2,show.legend=FALSE)) +
  ggdist::stat_halfeye(adjust=.5,width=.6,justification=-.35,.width=0,point_colour=NA) +
  geom_boxplot(width=.3,outlier.color=NA) +
  geom_point(color='black',size=1.5,alpha=.6,position=position_jitter(seed=1,width=.1))+
  geom_signif(comparisons=list(c("Baseline","Ayahuasca")),test='wilcox.test',test.args=list(paired=TRUE),map_signif_level=TRUE,y_position=c(0.022, 0.018, 0.017))+
  geom_signif(comparisons=list(c("Chronic/Sober","Chronic/Cannabis")),test='wilcox.test',map_signif_level=TRUE,y_position=c(0.018))+
  geom_signif(comparisons=list(c("Occasional/Sober","Occasional/Cannabis")),annotation=c('***'),y_position=c(0.017))+
  geom_signif(comparisons=list(c("Placebo","DMT")),test='wilcox.test',map_signif_level=TRUE,test.args = list(paired=FALSE), y_position=0.027)+
  theme_classic()+
  scale_fill_brewer(palette="Pastel1") + scale_color_brewer(palette = "Pastel1")+
  theme(axis.text.x=element_text(angle=45,vjust=0.95,hjust=0.95)) +
  #theme(text=element_text(family="SF Pro"))+
  theme(text=element_text(size=20))+
  ylim(0,0.03)+theme(legend.position="none")+xlab("")

ggsave(filename='DIRCHANGE.png',plot=last_plot(),path='/Users/myco/Desktop/',device=ragg::agg_png(),dpi=320,width=777,height=588,units=c('px'))

# Mean hierarchical levels
change<-read.csv('/Users/myco/human_brains/LSDAyaCannab/changeHL.csv',header=TRUE)
colnames(change)[1] = 'Change'
change$Condition<-factor(change$Condition,levels=c("Ayahuasca","DMT","Chronic","Occasional"))


outliers <- change %>% group_by(Condition) %>% identify_outliers(Change)
data.frame(outliers)

#change<-change[-which(change$Change %in% outliers$Change[outlier$is.extreme=='TRUE']),] # not necessary


change %>% ggplot(aes(Condition,Change,fill=Condition, show.legend=FALSE)) +
  ggdist::stat_halfeye(adjust=.5,width=.6,justification =-.35,.width=0,point_colour=NA) +
  geom_boxplot(width=.3,outlier.color=NA,show.legend=FALSE) +
  geom_point(size=1.5,alpha=.6,aes(colour=Condition),position=position_jitter(seed=1,width=.1),show.legend=FALSE)+#geom_text(hjust=0, vjust=0)+
  theme_classic()+
#  theme(text=element_text(family="SF Pro")) +
  theme(text=element_text(size=20)) +
  scale_fill_brewer(palette="Pastel1") + scale_color_brewer(palette="Pastel1")+
  scale_y_continuous(breaks=pretty_breaks(n=5))+
  theme(legend.position="none")+xlab("")+ylab("Hierarchical level change from baseline")

# Mean HL change update
hlsig<-read.csv('/Users/myco/human_brains/LSDAyaCannab/HLCHANGEBIG.csv',header=TRUE)
hlsig$Condidv<-factor(hlsig$Condidv,levels=c("Baseline","Ayahuasca","Placebo","DMT","Baseline/Chronic","Cannabis/Chronic","Baseline/Occasional","Cannabis/Occasional"))
#hlsig$Condition<-factor(sigdat$Condition,levels=c("Placebo","DMT"))

hlsig %>% ggplot(aes(Condidv,HL,fill=Condition))+
  geom_violin(trim=FALSE,show.legend=FALSE)+
  geom_point(aes(group=Region_group,fill=Condition),position=position_dodge(0.2),size=2,shape=1,alpha=0.6,show.legend = FALSE) +
  geom_line(aes(group=Region_group),position=position_dodge(0.2),size=0.5,alpha=0.5,color='black')+
  theme_classic()+
  scale_fill_brewer(palette="Pastel1")+scale_color_brewer(palette="Pastel1")+
  theme(text=element_text(size=20))+
  scale_y_continuous(breaks=pretty_breaks(n=5))+
  xlab("")+ylab("Hierarchical level")+
  theme(axis.text.x=element_text(angle=45,vjust=0.95,hjust=0.95))
  

# hlsig %>% ggplot(aes(Condition,HL,fill=Condition)) +
#   geom_violin(trim=FALSE,show.legend=FALSE) +
#   geom_point(aes(fill=Condition,group=Region),position = position_dodge(0.2),size=2,shape=1,alpha=0.6,show.legend=FALSE) +
#   geom_line(aes(group=Region),position=position_dodge(0.2),size=0.5,alpha=0.5,color="black") +
#   # geom_signif(comparisons=list(c("Placebo","DMT")), test='wilcox.test', map_signif_level=TRUE, y_position=0.3)+
#   theme_classic()+
#   scale_fill_brewer(palette="Pastel1") + scale_color_brewer(palette="Pastel1") +
#   #scale_x_discrete(limits=c("Placebo","DMT")) +
#   # theme(text=element_text(family="SF Pro")) +
#   theme(text=element_text(size=17)) +
#   ylim(-.15, 0.35) + ylab('Change in hierarchical levels')
# scale_y_continuous(breaks=pretty_breaks(n=3),expand = c(0, 0), limits = c(0, NA))

# # Regional hierarchical levels
# 
# rsnhl<-read.csv('/Users/myco/human_brains/LSDAyaCannab/rsnhlsigchanges.csv',header=FALSE)
# rsnhl<-data.frame(as.vector(rsnhl))
# colnames(rsnhl) <- c("Visual","Somatomotor","Dorsal attention","Salience","Limbic","Frontoparietal","Default mode")
# rownames(rsnhl) <- c("Ayahuasca","DMT","Cannabis (Chronic)","Cannabis (Occasional)")
# rsnhl<-rsnhl %>% rownames_to_column(var="Condition")
# rsnhl$Condition<-factor(rsnhl$Condition,levels=c("Ayahuasca","DMT","Cannabis (Chronic)","Cannabis (Occasional)"))
# 
# 
# rsnhlmelt<-melt(rsnhl,id='Condition')
# 
# rsnhlmelt %>% ggplot(aes(x=variable,y=value,fill=Condition)) +
#   geom_col(position='dodge',color='black',width=.7)+
#   #coord_flip()
#   theme_minimal()+
#   scale_fill_brewer(palette="Pastel1")+scale_color_brewer(palette="Pastel1")+
#   scale_y_continuous(breaks=pretty_breaks(n=5))+
#   theme(text=element_text(size=20))+
#   #theme(text=element_text(family="SF Pro")) +
#   ylab("Change in trophic level")+xlab("")+
#   theme(axis.text.x=element_text(angle=45,vjust=0.95,hjust=0.95))+
#   theme(axis.line = element_line(color="black"),
#         axis.ticks = element_line(color="black"),
#         panel.border = element_blank())
# 
# ggsave(filename='RSNCHANGES4.svg',device='svg',plot=last_plot(),path='/Users/myco/human_brains/INSIDEOUT_THESIS/',width=20,height=10)
# 
#   
#   


