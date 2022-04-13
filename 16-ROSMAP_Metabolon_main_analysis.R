library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd("~/ROSMAP_metabolon")

metab_data_transformed <- readRDS("data/14-metab_data_transformed.rds")
pheno <- readRDS("data/14-ROSMAP_pheno_clean.rds")
metab_meta <- readRDS("data/01-metab_meta.rds")

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

### Check metabolite associations with age
## Linear regression modeling metabolite values by age, sex, and PMI with only sAD samples
model_data_ADonly <- 
  model_data %>%
  filter(Status == "AD")

get_effect_pval_age <- function(metab_name, model_data) {
  readings <- cbind(model_data[,c("Status", "age_death", "pmi", "msex")], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  model <- glm(reading ~ age_death + msex + pmi,
               data = readings, family = gaussian)
  data.frame(
    metab_name = metab_name,
    effect = as.matrix(summary(model)$coefficients)["age_death","Estimate"],
    pval = as.matrix(summary(model)$coefficients)["age_death","Pr(>|t|)"]
  )
}

effect_pval_age <- ldply(colnames(model_data_ADonly[,-1:-20]), 
                         get_effect_pval_age,
                         model_data = model_data_ADonly)
effect_pval_age$padj <- p.adjust(effect_pval_age$pval, method = "BH")
all(effect_pval_age$metab_name == metab_meta$SHORT_NAME)
effect_pval_age$CHEMICAL.ID <- metab_meta$CHEM_ID
effect_pval_age <- 
  rename(effect_pval_age,
         "effect_ROSMAP_AD" = "effect",
         "pval_ROSMAP_AD" = "pval",
         "padj_ROSMAP_AD" = "padj",
         "metab_id" = "CHEMICAL.ID")
# write.csv(effect_pval_age, "data/03-effect_pval_age.csv", row.names = FALSE)

## Check association with age in only <90yo 
model_data_ADonly_under90 <- 
  model_data %>%
  filter(age_death < 90)

effect_pval_age_under90 <- ldply(colnames(model_data_ADonly_under90[,-1:-20]), 
                                 get_effect_pval_age,
                                 model_data = model_data_ADonly_under90)
effect_pval_age_under90$padj <- p.adjust(effect_pval_age_under90$pval, method = "BH")
all(effect_pval_age_under90$metab_name == metab_meta$SHORT_NAME)
effect_pval_age_under90$CHEMICAL.ID <- metab_meta$CHEM_ID
effect_pval_age_under90 <- 
  rename(effect_pval_age_under90,
         "effect_ROSMAP_AD_under90" = "effect",
         "pval_ROSMAP_AD_under90" = "pval",
         "padj_ROSMAP_AD_under90" = "padj",
         "metab_id" = "CHEMICAL.ID")

## Check assn with age in CO
model_data_COonly <- 
  model_data %>%
  filter(Status == "CO")

effect_pval_age_CO <- ldply(colnames(model_data_COonly[,-1:-20]), 
                            get_effect_pval_age,
                            model_data = model_data_COonly)
effect_pval_age_CO$padj <- p.adjust(effect_pval_age_CO$pval, method = "BH")
all(effect_pval_age_CO$metab_name == metab_meta$SHORT_NAME)
effect_pval_age_CO$CHEMICAL.ID <- metab_meta$CHEM_ID
effect_pval_age_CO <- 
  rename(effect_pval_age_CO,
         "effect_ROSMAP_CO" = "effect",
         "pval_ROSMAP_CO" = "pval",
         "padj_ROSMAP_CO" = "padj",
         "metab_id" = "CHEMICAL.ID")

## Check in only <90 yo
effect_pval_age_CO_under90 <- ldply(colnames(model_data_COonly[,-1:-20]), 
                                    get_effect_pval_age,
                                    model_data = model_data_COonly[model_data_COonly$age_death < 90,])
effect_pval_age_CO_under90$padj <- p.adjust(effect_pval_age_CO_under90$pval, method = "BH")
all(effect_pval_age_CO_under90$metab_name == metab_meta$SHORT_NAME)
effect_pval_age_CO_under90$CHEMICAL.ID <- metab_meta$CHEM_ID
effect_pval_age_CO_under90 <- 
  rename(effect_pval_age_CO_under90,
         "effect_ROSMAP_CO_under90" = "effect",
         "pval_ROSMAP_CO_under90" = "pval",
         "padj_ROSMAP_CO_under90" = "padj",
         "metab_id" = "CHEMICAL.ID")
