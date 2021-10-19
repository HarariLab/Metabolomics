library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggsignif)
library(MASS)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

## Ordinal logistic regression
metab_data_impute <- readRDS("data/15-metab_data_impute.rds")
pheno <- readRDS("data/14-ROSMAP_pheno_clean.rds")
pheno_ad <- pheno %>% filter(Status %in% c("AD", "CO"))
metab_data_impute_ad <- metab_data_impute[pheno_ad$individualID,]

metab_names <- c(
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

metab_names[which(!metab_names %in% colnames(metab_data_impute_ad))]
metab_names <- metab_names[-which(!metab_names %in% colnames(metab_data_impute_ad))]

metabs_15 <- metab_data_impute_ad[pheno_ad$individualID,metab_names]

# Get PC1
pca_obj <- PCA(metabs_15)
pca_ind <- pca_obj$ind$coord

# Add PC1 to pheno info and plot
pheno_ad$Status <- factor(pheno_ad$Status, levels = c("CO", "AD"))

all(row.names(pca_ind) == pheno_ad$individualID)
pheno_ad$PC1 <- pca_ind[,1]

ggplot(pheno_ad, aes(x = Status, y = PC1, fill = Status)) +
  geom_boxplot() +
  geom_signif(
    comparisons = list(c("CO", "AD")),
    y_position = c(6)
    , tip_length = 0
  ) +
  theme_minimal() +
  scale_x_discrete(labels = c("CO", "sAD")) + 
  scale_fill_manual(values=c("#F8766D", "#00C19F")) +
  theme(legend.position = "none")

# Test eigengene between groups
eigen_model <- glm(Status ~ PC1 + age_death + msex + pmi, data = pheno_ad, family = binomial)
summary(eigen_model)

# Check association with duration
pheno_ad$age_first_ad_dx <-
  pheno_ad$age_first_ad_dx %>%
  as.character() %>%
  recode(
    "90+" = "90"
  ) %>%
  as.numeric()

pheno_ad$duration <- pheno_ad$age_death - pheno_ad$age_first_ad_dx

# Remove age at onset of 90+ (unknown duration)
pheno_duration <- pheno_ad %>% filter(age_first_ad_dx < 90)
pheno_duration$duration_adjusted <- pheno_duration$duration/mean(pheno_duration$duration, na.rm = TRUE)
mod1 <- glm(PC1 ~ duration_adjusted + age_death + msex + pmi, data = pheno_duration, family = gaussian)
summary(mod1)

# Check association with Braak tau
mod2 <- glm(PC1 ~ braaksc + age_death + msex + pmi, data = pheno_ad, family = gaussian)
summary(mod2)


# Test association with consensus clinical diagnosis
## Cognitive status
pheno$Status2 <-
  pheno$cogdx %>%
  recode(
    "1" = "CO",
    "2" = "MCI",
    "3" = "MCI",
    "4" = "AD",
    "5" = "AD",
    "6" = "Other"
  ) %>%
  as.factor()

pheno_cog_ad <- pheno %>% filter(Status2 %in% c("AD", "CO"))

table(pheno$Status, pheno$Status2)

metabs_15_cog <- metab_data_impute[pheno_cog_ad$individualID,metab_names]

# Get PC1
pca_obj_cog <- PCA(metabs_15_cog, graph = FALSE)
pca_ind_cog <- pca_obj_cog$ind$coord

# Add PC1 to pheno info and plot
pheno_cog_ad$Status2 <- factor(pheno_cog_ad$Status2, levels = c("CO", "AD"))

all(row.names(pca_ind_cog) == pheno_cog_ad$individualID)
pheno_cog_ad$PC1 <- pca_ind_cog[,1]

eigen_model_cog <- glm(Status2 ~ PC1 + age_death + msex + pmi, data = pheno_cog_ad, family = binomial)
summary(eigen_model_cog)