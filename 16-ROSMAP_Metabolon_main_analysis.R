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
all(effect_pval_ad$metab_name == metab_meta$SHORT_NAME)
effect_pval_ad$CHEMICAL.ID <- metab_meta$CHEM_ID
effect_pval_ad$HMDB.ID <- metab_meta$HMDB

# write.csv(effect_pval_ad, "data/16-ROSMAP_Metabolon_effect_pval_ad.csv", row.names = FALSE)

# Check in common with ADRC
all_effect_pval_orig <- read.csv("data/03-all_pval_effect_clean_withage.csv", stringsAsFactors = FALSE)
joined_effect_pval <- full_join(all_effect_pval_orig, effect_pval_ad, by = "CHEMICAL.ID")
# write.csv(joined_effect_pval, "data/16-effect_pval_joined.csv", row.names = FALSE)

### Linear regression for APOE in sAD
model_data_apoe <-  
  model_data %>%
  filter(!is.na(apoe_status)) %>%
  filter(Status == "AD")


# Make a data frame with effect and p-value
get_effect_pval_APOE <- function(metab_name, model_data) {
  readings <- cbind(model_data[,1:20], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  model <- glm(reading ~ apoe_status + age_death + msex + pmi,
               data = readings, family = gaussian)
  data.frame(
    metab_name = metab_name,
    effect = as.matrix(summary(model)$coefficients)["apoe_statuspositive","Estimate"],
    pval = as.matrix(summary(model)$coefficients)["apoe_statuspositive","Pr(>|t|)"]
  )
}

effect_pval_apoe <- ldply(colnames(model_data_apoe[,-1:-20]), 
                          get_effect_pval_APOE,
                          model_data = model_data_apoe)
effect_pval_apoe$padj <- p.adjust(effect_pval_apoe$pval, method = "BH")

all(effect_pval_apoe$metab_name == metab_meta$SHORT_NAME)
effect_pval_apoe$CHEMICAL.ID <- metab_meta$CHEM_ID

colnames(effect_pval_apoe) <- c("metabolite", "effect_ROSMAP", "pval_ROSMAP", "padj_ROSMAP", "CHEMICAL.ID")

effect_pval_apoe_orig <- read.csv("data/04-effect_pval_apoe.csv", stringsAsFactors = F)

effect_pval_apoe_joined <- full_join(effect_pval_apoe_orig, effect_pval_apoe, by = "CHEMICAL.ID")

# write.csv(effect_pval_apoe, "data/16-effect_pval_apoe_ROSMAP.csv")
# write.csv(effect_pval_apoe_joined, "data/16-effect_pval_apoe_joined.csv")
