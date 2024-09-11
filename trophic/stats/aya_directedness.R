# loads, updates, installs packages
pacman::p_load(rstatix,tidyverse,scales,reshape2,data.table,cowplot,ggpubr,lme4,merTools,lmerTest,emmeans, purrr, R.matlab)


tc<-read.csv('/Users/myco/human_brains/ayahuasca/trophiccoherence.csv',header=FALSE)
condition<-read.csv('/Users/myco/human_brains/ayahuasca/condition.csv',header=FALSE)
subjlist<-read.csv('/Users/myco/human_brains/ayahuasca/subjlist.csv',header=FALSE)
dat<-data.frame(as.vector(tc),as.vector(condition),as.vector(subjlist))
colnames(dat)[1]="Directedness"
colnames(dat)[2]="Condition"
colnames(dat)[3]="Subj"
dat$Condition<-factor(dat$Condition,levels=c("Baseline","Ayahuasca"))

summary<-dat %>% group_by(Condition) %>% get_summary_stats(Directedness, type="mean_sd")
data.frame(summary)
dat<-na.omit(dat)

outlier <- dat %>% group_by(Condition) %>% identify_outliers(Directedness)
data.frame(outlier)

#outlier_free<-dat[-which(dat$Directedness %in% outlier$Directedness[outlier$is.outlier=='TRUE']),] # not necessary

normality <- dat %>% group_by(Condition) %>% shapiro_test(Directedness)
data.frame(normality) #wilcoxon necessary


qq1<-ggqqplot(dat$Directedness[dat$Condition=='Baseline'])
qq2<-ggqqplot(dat$Directedness[dat$Condition=='Ayahuasca'])
plot_grid(qq1,qq2,labels="AUTO")

# Directedness
dat %>% ggplot(aes(Condition,Directedness,fill=Condition)) +
  ggdist::stat_halfeye(adjust=.5,width=.6,justification=-.35,.width=0,point_colour=NA) +
  geom_boxplot(width=.3,outlier.color=NA) +
  geom_point(size=1.5,alpha=.6,color="black",position=position_jitter(seed=1,width=.1))+
  geom_signif(comparisons=list(c("Baseline","Ayahuasca")), test='wilcox.test', map_signif_level=TRUE, y_position=0.023) +
  theme_classic()+
  #theme(axis.text.x=element_text(angle=45,vjust=0.95,hjust=0.95)) +
  #theme(text=element_text(family="SF Pro")) +
  theme(text=element_text(size=25)) +
  scale_fill_brewer(palette="Pastel1") + scale_color_brewer(palette = "Dark2") +
  ylim(0,0.025)+
  scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme(legend.position="none")+xlab("")
ggsave(filename='AYADIR.svg',device='svg',plot=last_plot(),path='/Users/myco/human_brains/INSIDEOUT_THESIS/',width=300,height=300,units=c('px'))

# Mean hierarchical Levels
hlsig<-read.csv('/Users/myco/human_brains/ayahuasca/hlsig.csv',header=FALSE)
hlsigcond<-read.csv('/Users/myco/human_brains/ayahuasca/condition_hlsig.csv',header=FALSE)
regions<-read.csv('/Users/myco/human_brains/ayahuasca/regions_hlsig.csv',header=FALSE)
sigdat<-data.frame(as.vector(hlsig),as.vector(hlsigcond),as.vector(regions))
colnames(sigdat)[1]='Hierarchical Level'
colnames(sigdat)[2]='Condition'
colnames(sigdat)[3]="Region"
sigdat$Condition<-factor(sigdat$Condition,levels=c("Baseline","Ayahuasca"))

sigdat %>% ggplot(aes(Condition,`Hierarchical Level`,fill=Condition)) +
  geom_violin(trim=FALSE,show.legend=FALSE) +
  geom_point(aes(fill=Condition,group=Region),position = position_dodge(0.2),size=2,shape=1,alpha=0.6,show.legend=FALSE) +
  geom_line(aes(group=Region),position=position_dodge(0.2),size=0.5,alpha=0.5,color="black") +
  geom_signif(comparisons=list(c("Baseline","Ayahuasca")), test='wilcox.test', map_signif_level=TRUE, y_position=0.7)+
  theme_classic() +
  scale_fill_brewer(palette="Pastel1") + scale_color_brewer(palette="Dark2") +
  scale_x_discrete(limits=c("Baseline","Ayahuasca")) +
  theme(text=element_text(family="SF Pro")) +
  theme(text=element_text(size=17)) +
  ylim(0, 0.75) + ylab('Hierarchical Level')
  #scale_y_continuous(breaks=pretty_breaks(n=3),expand = c(0, 0), limits = c(0, NA))

# All hierarchical levels
hlfull<-read.csv('/Users/myco/human_brains/ayahuasca/hierarchicallevels_full.csv',header=FALSE)
sex<-read.csv('/Users/myco/human_brains/ayahuasca/sex.csv',header=FALSE)
age<-read.csv('/Users/myco/human_brains/ayahuasca/age.csv',header=FALSE)
lifetime<-read.csv('/Users/myco/human_brains/ayahuasca/lifetime.csv',header=FALSE)
recency<-read.csv('/Users/myco/human_brains/ayahuasca/recency.csv',header=FALSE)
dose<-read.csv('/Users/myco/human_brains/ayahuasca/dose.csv',header=FALSE)
names(hlfull)<-c(1:80)
regions<-read.csv('/Users/myco/human_brains/atlases/dbs80/dbs80symm_labels.csv',header=FALSE)$V1
subjlist<-read.csv('/Users/myco/human_brains/ayahuasca/subjlist.csv',header=FALSE)
hlfulldat<-data.frame(as.vector(hlfull),as.vector(condition),as.factor(subjlist$V1),as.factor(sex$V1),as.factor(age$V1),as.factor(lifetime$V1),as.factor(recency$V1),as.vector(dose))
colnames(hlfulldat) <- paste0("old",1:82)
hlfulldat <- hlfulldat %>% rename("Condition"=81,"Subj"=82,"Sex"=83,"Age"=84,"Lifetime"=85,"Recency"=86,"Dose"=87)
hlfulldat <- hlfulldat %>% rename_with(~ regions,starts_with("old"))

hlfulldatmelt<-melt(setDT(hlfulldat),id.vars=c("Condition","Subj","Sex","Age","Lifetime","Recency","Dose"))
colnames(hlfulldatmelt)[8]="Region"
colnames(hlfulldatmelt)[9]="Hierarchical Level"
summary<-hlfulldatmelt %>% group_by(Condition,Region) %>% get_summary_stats(`Hierarchical Level`, type="mean_sd")
data.frame(summary)
hlfulldatmelt<-na.omit(hlfulldatmelt)

outlier <- hlfulldatmelt %>% group_by(Condition,Region) %>% identify_outliers(`Hierarchical Level`)
data.frame(outlier)
hlfulldatmelt<-hlfulldatmelt[-which(hlfulldatmelt$`Hierarchical Level` %in% outlier$`Hierarchical Level`[outlier$is.extreme=='TRUE']),] # not necessary

colnames(hlfulldatmelt)[9]="Hierarchical Level"

ggplot(hlfulldatmelt,aes(Region,`Hierarchical Level`)) +
  geom_boxplot(aes(fill=Condition),show.legend=FALSE,outlier.color=NA) +
  facet_grid(rows=vars(Condition)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=70,vjust=1,hjust=1))+
  scale_fill_brewer(palette="Pastel1")+
  theme(text=element_text(family="SF Pro")) +
  theme(text=element_text(size=17))+
  xlab("")

# Mean hierarchical levels (RSN) NOT USED
rsnhl<-read.csv('/Users/myco/human_brains/ayahuasca/sighlrsn.csv',header=FALSE)
rsncond<-read.csv('/Users/myco/human_brains/ayahuasca/conditionrsn.csv',header=FALSE)
nsigrsn<-read.csv('/Users/myco/human_brains/ayahuasca/nsigrsn.csv',header=FALSE)
rsnsigdat<-data.frame(as.vector(rsnhl),as.vector(rsncond),as.vector(nsigrsn))
colnames(rsnsigdat)[1]='Hierarchical Level'
colnames(rsnsigdat)[2]='Condition'
colnames(rsnsigdat)[3]="Region"
rsnsigdat$Condition<-factor(rsnsigdat$Condition,levels=c("Baseline","Ayahuasca"))

rsnsigdat %>% ggplot(aes(Condition,`Hierarchical Level`,fill=Condition)) +
  geom_violin(trim=FALSE,show.legend=FALSE) +
  geom_point(aes(fill=Condition,group=Region),position = position_dodge(0.2),size=2,shape=1,alpha=0.6,show.legend=FALSE) +
  geom_line(aes(group=Region),position=position_dodge(0.2),size=0.5,alpha=0.5,color="black") +
  geom_signif(comparisons=list(c("Baseline","Ayahuasca")), test='wilcox.test', map_signif_level=TRUE, y_position=0.7)+
  theme_classic() +
  scale_fill_brewer(palette="Pastel1") + scale_color_brewer(palette="Dark2") +
  scale_x_discrete(limits=c("Baseline","Ayahuasca")) +
  theme(text=element_text(family="SF Pro")) +
  theme(text=element_text(size=17)) +
  ylim(0, 0.75) + ylab('Hierarchical Level')
#scale_y_continuous(breaks=pretty_breaks(n=3),expand = c(0, 0), limits = c(0, NA))

# All hierarchical levels (RSN)
hlfull<-read.csv('/Users/myco/human_brains/ayahuasca/hierarchicallevels_full.csv',header=FALSE)
sex<-read.csv('/Users/myco/human_brains/ayahuasca/sex.csv',header=FALSE)
age<-read.csv('/Users/myco/human_brains/ayahuasca/age.csv',header=FALSE)
lifetime<-read.csv('/Users/myco/human_brains/ayahuasca/lifetime.csv',header=FALSE)
recency<-read.csv('/Users/myco/human_brains/ayahuasca/recency.csv',header=FALSE)
dose<-read.csv('/Users/myco/human_brains/ayahuasca/dose.csv',header=FALSE)
names(hlfull)<-c(1:80)
regions<-read.csv('/Users/myco/human_brains/atlases/dbs80/dbs80symm_labels.csv',header=FALSE)$V1
subjlist<-read.csv('/Users/myco/human_brains/ayahuasca/subjlist.csv',header=FALSE)
hlfulldat<-data.frame(as.vector(hlfull),as.vector(condition),as.factor(subjlist$V1),as.factor(sex$V1),as.factor(age$V1),as.factor(lifetime$V1),as.factor(recency$V1),as.vector(dose))
colnames(hlfulldat) <- paste0("old",1:82)
hlfulldat <- hlfulldat %>% rename("Condition"=81,"Subj"=82,"Sex"=83,"Age"=84,"Lifetime"=85,"Recency"=86,"Dose"=87)
hlfulldat <- hlfulldat %>% rename_with(~ regions,starts_with("old"))

hlfulldatmelt<-melt(setDT(hlfulldat),id.vars=c("Condition","Subj","Sex","Age","Lifetime","Recency","Dose"))
colnames(hlfulldatmelt)[8]="Region"
colnames(hlfulldatmelt)[9]="Hierarchical Level"
summary<-hlfulldatmelt %>% group_by(Condition,Region) %>% get_summary_stats(`Hierarchical Level`, type="mean_sd")
data.frame(summary)
hlfulldatmelt<-na.omit(hlfulldatmelt)

outlier <- hlfulldatmelt %>% group_by(Condition,Region) %>% identify_outliers(`Hierarchical Level`)
data.frame(outlier)
hlfulldatmelt<-hlfulldatmelt[-which(hlfulldatmelt$`Hierarchical Level` %in% outlier$`Hierarchical Level`[outlier$is.extreme=='TRUE']),] # not necessary

colnames(hlfulldatmelt)[9]="Hierarchical Level"

ggplot(hlfulldatmelt,aes(Region,`Hierarchical Level`)) +
  geom_boxplot(aes(fill=Condition),show.legend=FALSE,outlier.color=NA) +
  facet_grid(rows=vars(Condition)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=70,vjust=1,hjust=1))+
  scale_fill_brewer(palette="Pastel1")+
  theme(text=element_text(family="SF Pro")) +
  theme(text=element_text(size=17))+
  xlab("")
