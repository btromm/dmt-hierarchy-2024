# loads, updates, installs packages
pacman::p_load(rstatix,tidyverse,scales,reshape2,data.table,cowplot,ggpubr,lme4,merTools,lmerTest,emmeans)

# tc<-read.csv('/Users/myco/human_brains/DMT/trophiccoherence.csv',header=FALSE)
# condition<-read.csv('/Users/myco/human_brains/LSD/condition.csv',header=FALSE)
# subjlist<-read.csv('/Users/myco/human_brains/LSD/subjlist.csv',header=FALSE)
# dat<-data.frame(as.vector(tc),as.vector(condition),as.vector(t(subjlist)))
# colnames(dat)[1]="Directedness"
# colnames(dat)[2]="Condition"
# colnames(dat)[3]="Subj"
# dat$Condition<-factor(dat$Condition,levels=c("Placebo","LSD"))
# 
# summary<-dat %>% group_by(Condition) %>% get_summary_stats(Directedness, type="mean_sd")
# data.frame(summary)
# dat<-na.omit(dat)
# 
# outlier <- dat %>% group_by(Condition) %>% identify_outliers(Directedness)
# data.frame(outlier)
# 
# #outlier_free<-dat[-which(dat$Directedness %in% outlier$Directedness[outlier$is.outlier=='TRUE']),] # not necessary
# 
# normality <- dat %>% group_by(Condition) %>% shapiro_test(Directedness)
# data.frame(normality) #wilcoxon necessary
# 
# 
# qq1<-ggqqplot(dat$Directedness[dat$Condition=='Placebo'])
# qq2<-ggqqplot(dat$Directedness[dat$Condition=='LSD'])
# plot_grid(qq1,qq2,labels="AUTO")
# 
# # Directedness
# dat %>% ggplot(aes(Condition,Directedness,fill=Condition,show.legend=FALSE)) +
#   ggdist::stat_halfeye(adjust=.5,width=.6,justification=-.35,.width=0,point_colour=NA) +
#   geom_boxplot(width=.3,outlier.color=NA,show.legend=FALSE) +
#   geom_point(size=1.5,alpha=.6,color="black",position=position_jitter(seed=1,width=.1),show.legend=FALSE)+
#   geom_signif(comparisons=list(c("Placebo","LSD")), test='wilcox.test', map_signif_level=TRUE, y_position=0.015,show.legend=FALSE) +
#   theme_classic()+
#   #theme(axis.text.x=element_text(angle=45,vjust=0.95,hjust=0.95)) +
#   theme(text=element_text(family="SF Pro")) +
#   theme(text=element_text(size=20)) +
#   scale_fill_brewer(palette="Pastel1") + scale_color_brewer(palette = "Pastel1") +
#   ylim(0,0.025)+
#   scale_y_continuous(breaks=pretty_breaks(n=3))+
#   theme(legend.position="none")

# Mean hierarchical Levels
hlsig<-read.csv('/Users/myco/human_brains/DMT/hlsig.csv',header=FALSE)
hlsigcond<-read.csv('/Users/myco/human_brains/DMT/condition_hlsig.csv',header=FALSE)
regions<-read.csv('/Users/myco/human_brains/DMT/regions_hlsig.csv',header=FALSE)
sigdat<-data.frame(as.vector(hlsig),as.vector(hlsigcond),as.vector(regions))
colnames(sigdat)[1]='Hierarchical Level'
colnames(sigdat)[2]='Condition'
colnames(sigdat)[3]="Region"
sigdat$Condition<-factor(sigdat$Condition,levels=c("Placebo","DMT"))

sigdat %>% ggplot(aes(Condition,`Hierarchical Level`,fill=Condition)) +
  geom_violin(trim=FALSE,show.legend=FALSE) +
  geom_point(aes(fill=Condition,group=Region),position = position_dodge(0.2),size=2,shape=1,alpha=0.6,show.legend=FALSE) +
  geom_line(aes(group=Region),position=position_dodge(0.2),size=0.5,alpha=0.5,color="black") +
 # geom_signif(comparisons=list(c("Placebo","DMT")), test='wilcox.test', map_signif_level=TRUE, y_position=0.3)+
  theme_classic()+
  scale_fill_brewer(palette="Pastel1") + scale_color_brewer(palette="Pastel1") +
  scale_x_discrete(limits=c("Placebo","DMT")) +
 # theme(text=element_text(family="SF Pro")) +
  theme(text=element_text(size=17)) +
  ylim(-.15, 0.35) + ylab('Change in hierarchical levels')
scale_y_continuous(breaks=pretty_breaks(n=3),expand = c(0, 0), limits = c(0, NA))

# All hierarchical levels
hlfull<-read.csv('/Users/myco/human_brains/dmt_update/hierarchicallevels_full.csv',header=FALSE)
names(hlfull)<-c(1:80)
regions<-read.csv('/Users/myco/human_brains/atlases/dbs80/dbs80symm_labels.csv',header=FALSE)$V1
subjlist<-read.csv('/Users/myco/human_brains/dmt_update/subj.csv',header=FALSE)
condition<-read.csv('/Users/myco/human_brains/dmt_update/condition.csv',header=FALSE)
hlfulldat<-data.frame(as.vector(hlfull),as.vector(condition),as.vector(subjlist))
colnames(hlfulldat) <- paste0("old",1:82)
hlfulldat <- hlfulldat %>% rename("Condition"=81,"Subj"=82)
hlfulldat <- hlfulldat %>% rename_with(~ regions,starts_with("old"))

hlfulldatmelt<-melt(setDT(hlfulldat),id.vars=c("Condition","Subj"))
colnames(hlfulldatmelt)[3]="Region"
colnames(hlfulldatmelt)[4]="Hierarchical Level"
summary<-hlfulldatmelt %>% group_by(Condition,Region) %>% get_summary_stats(`Hierarchical Level`, type="mean_sd")
data.frame(summary)
hlfulldatmelt<-na.omit(hlfulldatmelt)

outlier <- hlfulldatmelt %>% group_by(Condition,Region) %>% identify_outliers(`Hierarchical Level`)
data.frame(outlier)
hlfulldatmelt<-hlfulldatmelt[-which(hlfulldatmelt$`Hierarchical Level` %in% outlier$`Hierarchical Level`[outlier$is.extreme=='TRUE']),] # not necessary

colnames(hlfulldatmelt)[4]="Hierarchical Level"

ggplot(hlfulldatmelt,aes(Region,`Hierarchical Level`)) +
  geom_boxplot(aes(fill=Condition),show.legend=FALSE,outlier.color=NA) +
  facet_grid(rows=vars(Condition)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=70,vjust=1,hjust=1))+
  scale_fill_brewer(palette="Pastel1")+
 # theme(text=element_text(family="SF Pro")) +
  theme(text=element_text(size=17))+
  xlab("Region")
