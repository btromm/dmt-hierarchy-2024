# loads, updates, installs packages
pacman::p_load(lme4, tidyverse, merTools, rstatix, ggpubr, mice, mitml, cowplot, emmeans)

tc<-read.csv('/Users/myco/human_brains/cannab/trophiccoherence.csv',header=FALSE)
condition<-read.csv('/Users/myco/human_brains/cannab/condition.csv',header=FALSE)
subjlist<-read.csv('/Users/myco/human_brains/cannab/subjlist.csv',header=FALSE)
dat<-data.frame(as.vector(tc),as.vector(condition),as.vector(subjlist))
colnames(dat)[1]='tc'
colnames(dat)[2]='cond'
colnames(dat)[3]='subj'

summary<-dat %>% group_by(cond) %>% get_summary_stats(tc, type="mean_sd")
data.frame(summary)

outlier<-dat %>% group_by(cond) %>% identify_outliers(tc)
data.frame(outlier)

outlier_free<-dat[-which(dat$tc %in% outlier$tc[outlier$is.extreme=='TRUE']),]
#outlier_free$frq<-factor(outlier_free$frq,levels=c("Occasional","Chronic"))

# shapiro-wilk is a bit too sensitive here
normality <- dat %>% group_by(cond) %>% shapiro_test(tc)
data.frame(normality)

# using qqplots instead
qq1<-ggqqplot(outlier_free$tc[outlier_free$cond=='Chronic/Cannabis'])
qq2<-ggqqplot(outlier_free$tc[outlier_free$cond=='Chronic/Sober'])
qq3<-ggqqplot(outlier_free$tc[outlier_free$cond=='Occasional/Cannabis'])
qq4<-ggqqplot(outlier_free$tc[outlier_free$cond=='Occasional/Sober'])
plot_grid(qq1,qq2,qq3,qq4,labels="AUTO")
# all look good

library(mice)
# since there is missing data, we will use multiple imputation
# combined with a mixed effects model to test the differences
ini<-mice(data=outlier_free,maxit=0)
pred<-ini$pred
pred["tc","subj"]<-2
imp<-mice(data=outlier_free,pred=pred,m=100)
summary(imp)
implist<-mids2mitml.list(imp)
fit3<-with(implist,lmer(tc~1+cond+(1|subj)))
testEstimates(fit3)
fit3.reduced<-with(implist,lmer(tc~1+(1|subj)))
testModels(fit3,fit3.reduced,method="D1") #sig

gg<-lmer(tc~cond+(1|subj),data=outlier_free)

# follow-up testing
posthoc<-emmeans(gg,list(pairwise~cond),adjust="bonferroni")

outlier_free %>% ggplot(aes(cond,tc,fill=cond)) +
  ggdist::stat_halfeye(adjust=.5,width=.6,justification=-.35,.width=0,point_colour=NA) +
  geom_boxplot(width=.3,outlier.color=NA,show.legend=FALSE) +
  geom_point(size=1.5,alpha=.6,aes(colour=cond),position=position_jitter(seed=1,width=.1),show.legend=FALSE)+
  geom_signif(comparisons=list(c("Chronic/Cannabis","Occasional/Cannabis"),c("Occasional/Cannabis","Occasional/Sober")), test='wilcox.test', map_signif_level=TRUE,y_position=c(0.018,0.017),show.legend=FALSE) +
  theme_classic()+
  theme(axis.text.x=element_text(angle=45,vjust=0.95,hjust=0.95))+
  theme(text=element_text(family="SF Pro"))+
  theme(text=element_text(size=20)) +
  scale_fill_brewer(palette="Pastel1") + scale_color_brewer(palette = "Pastel1")+
  ylim(0,0.1)+
  scale_y_continuous(breaks=pretty_breaks(n=3))+
  theme(legend.position="none")+xlab("")+ylab("Directedness")

#p + xlab('Condition') + ylab('Directedness') + scale_fill_brewer(palette="Pastel1")
# need to add in significance manually
p <- ggplot(outlier_free,aes(x=cond,y=tc,fill=frq)) + geom_violin(trim=FALSE, position=position_dodge(1)) + theme_minimal()
p + scale_fill_brewer(palette="Pastel1") + ylim(NA,0.04) + ylab('Directedness') + xlab('Condition') + scale_x_discrete(limits=c("Placebo","Cannabis")) + geom_point(position = position_jitter(seed = 1, width = 0.2)) #+ stat_summary(fun.data="mean_sdl",geom="pointrange",color="blue")

# All hierarchical levels
hlfull<-read.csv('/Users/myco/human_brains/cannab/hierarchicallevels_full.csv',header=FALSE)
names(hlfull)<-c(1:80)
regions<-read.csv('/Users/myco/human_brains/atlases/dbs80/dbs80symm_labels.csv',header=FALSE)$V1
condition<-read.csv('/Users/myco/human_brains/cannab/hl_condition.csv',header=FALSE)
#subjlist<-read.csv('/Users/myco/human_brains/cannab/subjlist.csv',header=FALSE)
hlfulldat<-data.frame(as.vector(hlfull),as.vector(condition))
colnames(hlfulldat) <- paste0("old",1:81)
hlfulldat <- hlfulldat %>% rename("Condition"=81)
hlfulldat <- hlfulldat %>% rename_with(~ regions,starts_with("old"))

hlfulldatmelt<-melt(setDT(hlfulldat),id.vars=c("Condition"))
colnames(hlfulldatmelt)[2]="Region"
colnames(hlfulldatmelt)[3]="Hierarchical Level"
summary<-hlfulldatmelt %>% group_by(Condition,Region) %>% get_summary_stats(`Hierarchical Level`, type="mean_sd")
data.frame(summary)
hlfulldatmelt<-na.omit(hlfulldatmelt)

outlier <- hlfulldatmelt %>% group_by(Condition,Region) %>% identify_outliers(`Hierarchical Level`)
data.frame(outlier)
hlfulldatmelt<-hlfulldatmelt[-which(hlfulldatmelt$`Hierarchical Level` %in% outlier$`Hierarchical Level`[outlier$is.outlier=='TRUE']),] # not necessary

colnames(hlfulldatmelt)[3]="Hierarchical Level"

ggplot(hlfulldatmelt,aes(Region,`Hierarchical Level`)) +
  geom_boxplot(aes(fill=Condition),show.legend=FALSE,outlier.color=NA) +
  facet_grid(rows=vars(Condition)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=70,vjust=1,hjust=1))+
  scale_fill_brewer(palette="Pastel1")+
  theme(text=element_text(family="SF Pro")) +
  theme(text=element_text(size=17))+
  xlab("Region")

# Mean hierarchical Levels
hlsig<-read.csv('/Users/myco/human_brains/cannab/hlsig.csv',header=FALSE)
hlsigcond<-read.csv('/Users/myco/human_brains/cannab/condition_hlsig.csv',header=FALSE)
regions<-read.csv('/Users/myco/human_brains/cannab/regions_hlsig.csv',header=FALSE)
sigdat<-data.frame(as.vector(hlsig),as.vector(hlsigcond),as.vector(regions))
colnames(sigdat)[1]='Hierarchical Level'
colnames(sigdat)[2]='Condition'
colnames(sigdat)[3]="Region"
sigdat$Condition<-factor(sigdat$Condition,levels=c("Chronic","Occasional"))

sigdat %>% ggplot(aes(Condition,`Hierarchical Level`,fill=Condition)) +
  geom_violin(trim=FALSE,show.legend=FALSE) +
  geom_point(aes(fill=Condition,group=Region),position = position_dodge(0.2),size=2,shape=1,alpha=0.6,show.legend=FALSE) +
  geom_line(aes(group=Region),position=position_dodge(0.2),size=0.5,alpha=0.5,color="black") +
  geom_signif(comparisons=list(c("Chronic","Occasional")), test='wilcox.test', map_signif_level=TRUE, y_position=0.83)+
  theme_classic() +
  scale_fill_brewer(palette="Pastel1") + scale_color_brewer(palette="Pastel1") +
  scale_x_discrete(limits=c("Occasional","Chronic")) +
  theme(text=element_text(family="SF Pro")) +
  theme(text=element_text(size=17)) +
  #ylim(0, 0.8) + 
  ylab('Hierarchical Levels')
