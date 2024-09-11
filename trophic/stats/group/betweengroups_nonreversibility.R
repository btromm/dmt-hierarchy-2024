pacman::p_load(tidyverse, rstatix, ggdist, ggsignif)

# Nonreversibility
NR <- read.csv("/Users/myco/human_brains/LSDAyaCannab/bigNR.csv", header = FALSE)
condition <- read.csv("/Users/myco/human_brains/LSDAyaCannab/condition_cond.csv", header = FALSE)
condition2 <- read.csv("/Users/myco/human_brains/LSDAyaCannab/condition_indiv.csv", header = FALSE)
dat <- data.frame(as.vector(NR), as.vector(condition), as.vector(condition2))
colnames(dat)[1] <- "Irreversibility"
colnames(dat)[2] <- "Condition1"
colnames(dat)[3] <- "Condition2"


summary <- dat %>%
  group_by(Condition1) %>%
  get_summary_stats(Irreversibility, type = "mean_sd")
data.frame(summary)

outlier <- dat %>%
  group_by(Condition1) %>%
  identify_outliers(Irreversibility)
data.frame(outlier)

outlier_free <- dat[-which(dat$Irreversibility %in% outlier$Irreversibility[outlier$is.extreme == "TRUE"]), ] # not necessary

normality <- dat %>%
  group_by(Condition1) %>%
  shapiro_test(Irreversibility)
data.frame(normality) # wilcoxon necessary
dat$Condition2 <- factor(dat$Condition2, levels = c("Baseline", "Ayahuasca", "Placebo", "DMT", "Sober/Chronic", "Cannabis/Chronic", "Sober/Occasional", "Cannabis/Occasional"))
outlier_free$Condition2 <- factor(outlier_free$Condition2, levels = c("Baseline", "Ayahuasca", "Placebo", "DMT", "Sober/Chronic", "Cannabis/Chronic", "Sober/Occasional", "Cannabis/Occasional"))

outlier_free %>% ggplot(aes(Condition2, Irreversibility, fill = Condition1)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, justification = -.35, .width = 0, point_colour = NA) +
  geom_boxplot(width = .3, outlier.color = NA) +
  geom_point(size = 1.5, alpha = .6, color = "black", position = position_jitter(seed = 1, width = .1)) +
  geom_signif(comparisons = list(c("Placebo", "DMT")), test = "wilcox.test", test.args = list(paired = TRUE), map_signif_level = TRUE, y_position = .0018) +
  geom_signif(comparisons = list(c("Sober/Chronic", "Cannabis/Chronic")), annotations = c("*"), y_position = 0.0017) +
  geom_signif(comparisons = list(c("Sober/Occasional", "Cannabis/Occasional")), test = "wilcox.test", map_signif_level = TRUE, y_position = 0.0014) +
  geom_signif(comparisons = list(c("Baseline", "Ayahuasca")), test = "wilcox.test", test.args = list(paired = FALSE), map_signif_level = TRUE, y_position = 0.0019) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95)) +
  # theme(text=element_text(family="SF Pro")) +
  theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_color_brewer(palette = "Pastel1") +
  theme(legend.position = "none") +
  xlab("")


# Hierarchy
hierarchy <- read.csv("/Users/myco/human_brains/LSDAyaCannab/bigHierarchy.csv", header = FALSE)
# condition<-read.csv('/Users/myco/human_brains/LSDAyaCannab/condition.csv',header=FALSE)
dathl <- data.frame(as.vector(hierarchy), as.vector(condition), as.vector(condition2))
colnames(dathl)[1] <- "Hierarchy"
colnames(dathl)[2] <- "Condition"
colnames(dathl)[3] <- "Condition2"

summary <- dathl %>%
  group_by(Condition) %>%
  get_summary_stats(Hierarchy, type = "mean_sd")
data.frame(summary)

outlier <- dathl %>%
  group_by(Condition) %>%
  identify_outliers(Hierarchy)
data.frame(outlier)

outlier_free <- dathl[-which(dathl$Hierarchy %in% outlier$Hierarchy[outlier$is.extreme == "TRUE"]), ] # not necessary

normality <- dathl %>%
  group_by(Condition) %>%
  shapiro_test(Hierarchy)
data.frame(normality) # wilcoxon necessary
dathl$Condition2 <- factor(dat$Condition2, levels = c("Baseline", "Ayahuasca", "Placebo", "DMT", "Sober/Chronic", "Cannabis/Chronic", "Sober/Occasional", "Cannabis/Occasional"))
# outlier_free$Condition<-factor(outlier_free$Condition,levels=c("Baseline","Ayahuasca","Placebo","LSD","Sober/Chronic","Cannabis/Chronic","Sober/Occasional","Cannabis/Occasional"))

dathl %>% ggplot(aes(Condition2, Hierarchy, fill = Condition)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, justification = -.35, .width = 0, point_colour = NA) +
  geom_boxplot(width = .3, outlier.color = NA) +
  geom_point(size = 1.5, alpha = .6, color = "black", position = position_jitter(seed = 1, width = .1)) +
  geom_signif(comparisons = list(c("Placebo", "DMT")), test = "wilcox.test", map_signif_level = TRUE, test.args = list(paired = TRUE), y_position = 0.05) +
  geom_signif(comparisons = list(c("Sober/Chronic", "Cannabis/Chronic")), test = "wilcox.test", map_signif_level = TRUE, y_position = 0.05) +
  geom_signif(comparisons = list(c("Sober/Occasional", "Cannabis/Occasional")), test = "wilcox.test", map_signif_level = TRUE, y_position = 0.048) +
  geom_signif(comparisons = list(c("Baseline", "Ayahuasca")), test = "wilcox.test", map_signif_level = TRUE, y_position = 0.048) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95)) +
  # theme(text=element_text(family="SF Pro")) +
  theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_color_brewer(palette = "Pastel1") +
  ylim(0, 0.065) +
  theme(legend.position = "none") +
  xlab("")
