## Test magnitude of effect between ADAD, sAD, and TREM2
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)

## Compare effect of ADvsCO and TREM2vsCO relative to effect of ADADvsCO (ANCOVA)
effect_pval <- read.csv("data/05-effect_pval_adad_ca_trem2.csv")
effect_pval_small <- effect_pval[effect_pval$ADADvsCO_padj < .05 & effect_pval$CAvsCO_pval < 0.05 & effect_pval$TREM2vsCO_pval < 0.05,] %>%
  filter(metab_name != "serotonin") # Remove serotonin due to high missingness

effect_small_long <- 
  data.frame(
    metab = rep(effect_pval_small$metab, 2),
    test = as.factor(c(rep("CA", 16), rep("TREM2", 16))),
    other_coef = c(effect_pval_small$CAvsCO_effect, effect_pval_small$TREM2vsCO_effect),
    ADAD_coef = rep(effect_pval_small$ADADvsCO_effect, 2)
  )

# If coefficient other_coef:test is significant, the slopes are significantly different
mod1 <- aov(ADAD_coef~other_coef*test, data=effect_small_long)
summary(mod1)

mod1a <- lm(ADAD_coef~other_coef*test, data=effect_small_long)
car::Anova(mod1a, type="II")

## Plot effects
metabs_16_df <- effect_pval_small %>% select(ADADvsCO_effect, TREM2vsCO_effect, CAvsCO_effect)

metabs_16_melt <- reshape2::melt(metabs_16_df, id.vars = 1)
ggplot(metabs_16_melt, aes(x=ADADvsCO_effect, y = value, color = variable)) + 
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), show.legend = FALSE) +
  scale_color_manual(name = element_blank(), labels = c("TREM2 vs CO effect", "sAD vs CO effect"), values = c("CAvsCO_effect" = "darkred", "TREM2vsCO_effect" = "steelblue")) +
  xlab("ADAD vs CO effect") +
  ylab("Other effect") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))
# ggsave("plots/05-effect_comparison_plot_16.png", device = "png", width = 8, height = 5)



### Re-check ANCOVA with matched results
metabs_16 <- c(
  "1,5-anhydroglucitol (1,5-AG)",
  "2-methylcitrate/homocitrate",
  "3-hydroxy-2-ethylpropionate",
  "alpha-tocopherol",
  "aspartate",
  "beta-citrylglutamate",
  "CDP-ethanolamine",
  "ergothioneine",
  "gamma-glutamylthreonine",
  "glutamate",
  "glutarate (C5-DC)",
  "glycerophosphoinositol*",
  "N-acetylglutamate",
  "N-carbamoylglutamate",
  "nicotinamide",
  "retinol (Vitamin A)")

############ High CDR
effect_pval_highCDR <- read.csv("data/04-highCDR_effect_pvals.csv")

# Get only the 16 chosen metabolites
effect_pval_small_highCDR <- effect_pval_highCDR %>%
  filter(metab %in% metabs_16) # Remove serotonin due to high missingness

# Restructure data frame to match AD and TREM2 effects to ADAD effects
effect_small_long_highCDR <- 
  data.frame(
    metab = rep(effect_pval_small_highCDR$metab, 2),
    test = as.factor(c(rep("CA", 16), rep("TREM2", 16))),
    other_coef = c(effect_pval_small_highCDR$CAvsCO_effect, effect_pval_small_highCDR$TREM2vsCO_effect),
    ADAD_coef = rep(effect_pval_small_highCDR$ADADvsCO_effect, 2)
  )

# Run ANCOVA model
mod_highCDR <- lm(ADAD_coef~other_coef*test, data=effect_small_long_highCDR)
car::Anova(mod_highCDR, type="II")

# Get only effects for plotting
metabs_16_highCDR <- effect_pval_small_highCDR %>% select(ADADvsCO_effect, CAvsCO_effect, TREM2vsCO_effect)

# Melt data frame and plot lines
metabs_16_melt_highCDR <- melt(metabs_16_highCDR, id.vars = 1)
ggplot(metabs_16_melt_highCDR, aes(x=ADADvsCO_effect, y = value, color = variable)) + 
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), show.legend = FALSE) +
  scale_color_manual(name = element_blank(), labels = c("Sporadic AD vs CO effect", "TREM2 vs CO effect"), values = c("CAvsCO_effect" = "darkred", "TREM2vsCO_effect" = "steelblue")) +
  xlab("ADAD vs CO effect") +
  ylab("Other effect") +
  ggtitle("Matched on CDR") +
  theme_bw()



########### High Tau
effect_pval_highTau <- read.csv("data/04-highTau_effect_pvals.csv")

effect_pval_small_highTau <- effect_pval_highTau %>%
  filter(metab %in% metabs_16) # Remove serotonin due to high missingness

effect_small_long_highTau <- 
  data.frame(
    metab = rep(effect_pval_small_highTau$metab, 2),
    test = as.factor(c(rep("CA", 16), rep("TREM2", 16))),
    other_coef = c(effect_pval_small_highTau$CAvsCO_effect, effect_pval_small_highTau$TREM2vsCO_effect),
    ADAD_coef = rep(effect_pval_small_highTau$ADADvsCO_effect, 2)
  )

# Run ANCOVA model
mod_highTau <- lm(ADAD_coef~other_coef*test, data=effect_small_long_highTau)
car::Anova(mod_highTau, type="II")

# Get only effects for plotting
metabs_16_highTau <- effect_pval_small_highTau %>% select(ADADvsCO_effect, CAvsCO_effect, TREM2vsCO_effect)

# Melt data frame to plot lines
metabs_16_melt_highTau <- melt(metabs_16_highTau, id.vars = 1)
ggplot(metabs_16_melt_highTau, aes(x=ADADvsCO_effect, y = value, color = variable)) + 
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), show.legend = FALSE) +
  scale_color_manual(name = element_blank(), labels = c("Sporadic AD vs CO effect", "TREM2 vs CO effect"), values = c("CAvsCO_effect" = "darkred", "TREM2vsCO_effect" = "steelblue")) +
  xlab("ADAD vs CO effect") +
  ylab("Other effect") +
  ggtitle("Matched on Braak Tau") +
  theme_bw()




############ High Abeta
effect_pval_highAbeta <- read.csv("data/04-highAbeta_effect_pvals.csv")

# Get only effects for metabolites in the 16
effect_pval_small_highAbeta <- effect_pval_highAbeta %>%
  filter(metab %in% metabs_16) # Remove serotonin due to high missingness

# Restructure data frame to match AD and TREM2 effects with ADAD effects
effect_small_long_highAbeta <- 
  data.frame(
    metab = rep(effect_pval_small_highAbeta$metab, 2),
    test = as.factor(c(rep("CA", 16), rep("TREM2", 16))),
    other_coef = c(effect_pval_small_highAbeta$CAvsCO_effect, effect_pval_small_highAbeta$TREM2vsCO_effect),
    ADAD_coef = rep(effect_pval_small_highAbeta$ADADvsCO_effect, 2)
  )

# Run ANCOVA
mod_highAbeta <- lm(ADAD_coef~other_coef*test, data=effect_small_long_highAbeta)
car::Anova(mod_highAbeta, type="II")

# Get only effects for plotting
metabs_16_highAbeta <- effect_pval_small_highAbeta %>% select(ADADvsCO_effect, CAvsCO_effect, TREM2vsCO_effect)

# Melt data frame to plot lines
metabs_16_melt_highAbeta <- melt(metabs_16_highAbeta, id.vars = 1)
ggplot(metabs_16_melt_highAbeta, aes(x=ADADvsCO_effect, y = value, color = variable)) + 
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), show.legend = FALSE) +
  scale_color_manual(name = element_blank(), labels = c("Sporadic AD vs CO effect", "TREM2 vs CO effect"), values = c("CAvsCO_effect" = "darkred", "TREM2vsCO_effect" = "steelblue")) +
  xlab("ADAD vs CO effect") +
  ylab("Other effect") +
  ggtitle("Matched on Braak Abeta") +
  theme_bw()














