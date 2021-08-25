library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd("~/ROSMAP_metabolon")

metab_data_transformed <- readRDS("data/14-metab_data_transformed.rds")
pheno <- readRDS("data/14-ROSMAP_pheno_clean.rds")

# Put pheno and metab data together
metab_data_transformed$individualID <- row.names(metab_data_transformed)
model_data <- right_join(pheno, metab_data_transformed)

# Make a data frame with effect and p-value (for brain, uses age at death)
get_effect_pval <- function(metab_name, model_data, status) {
  readings <- cbind(model_data[,c("Status", "age_death", "pmi", "msex")], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
    model <- glm(reading ~ Status + age_death + msex + pmi,
                 data = readings, family = gaussian)
    data.frame(
      metab_name = metab_name,
      effect_ROSMAP_metabolon = as.matrix(summary(model)$coefficients)[paste("Status", status, sep = ""),"Estimate"],
      pval_ROSMAP_metabolon = as.matrix(summary(model)$coefficients)[paste("Status", status, sep = ""),"Pr(>|t|)"]
    )
}

model_data_ad <- model_data %>% filter(model_data$Status %in% c("CO", "AD"))
model_data_ad$Status <- factor(model_data_ad$Status, levels = c("CO", "AD"))

effect_pval_ad <- ldply(colnames(model_data_ad[,21:ncol(model_data_ad)]),
                        get_effect_pval,
                        model_data = model_data_ad,
                        status = "AD")

effect_pval_ad$padj_ROSMAP_metabolon <- p.adjust(effect_pval_ad$pval_ROSMAP_metabolon, method = "BH")
effect_pval_ad <- effect_pval_ad[order(effect_pval_ad$padj_ROSMAP_metabolon),]

# write.csv(effect_pval_ad, "data/16-ROSMAP_Metabolon_effect_pval_ad.csv", row.names = FALSE)
