### Stratified analysis for EOAD vs LOAD
library(plyr)
library(dplyr)

model_data_duration <- readRDS("data/06-model_data_duration.rds") %>% select(-duration)
metab_meta <- readRDS("data/04-metab_meta.rds")

covariates <- model_data_duration %>% select( "TubeBarcode", "AAO", "Status", "age_death", "Sex", "PMI", "APOE", "BraakAbeta",  "BraakTau", "CDR")

## Define function for linear regression analysis
# Make a data frame with effect and p-value
get_effect_pval <- function(metab_id, model_data, status, EOAD = FALSE) {
  readings <- cbind(model_data[,1:ncol(covariates)], model_data[,metab_id])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  if(EOAD == TRUE) {
    model <- glm(reading ~ Status + Sex + PMI,
                 data = readings, family = gaussian)
  }
  else {
    model <- glm(reading ~ Status + age_death + Sex + PMI,
                 data = readings, family = gaussian)
  }
  data.frame(
    metab_id = metab_id,
    effect = as.matrix(summary(model)$coefficients)[paste("Status", status, sep = ""),"Estimate"],
    pval = as.matrix(summary(model)$coefficients)[paste("Status", status, sep = ""),"Pr(>|t|)"]
  )
}

## Linear regression for EOAD vs CO -- Not correcting for age
model_data_EOAD <- 
  model_data_duration %>% 
  filter((Status == "CA" & AAO < 65) | Status == "CO") %>%
  mutate(Status = factor(Status, levels = c("CO", "CA")))

EOAD_effect_pval_df <- ldply(colnames(model_data_EOAD[,-1:-ncol(covariates)]), 
                             get_effect_pval,
                             model_data = model_data_EOAD, 
                             status = "CA",
                             EOAD = TRUE)
EOAD_effect_pval_df$padj <- p.adjust(EOAD_effect_pval_df$pval, method = "BH")
if (all(EOAD_effect_pval_df$metab_id == metab_meta$CHEMICAL.ID)) EOAD_effect_pval_df$metab_name <- metab_meta$BIOCHEMICAL
sum(EOAD_effect_pval_df$padj < 0.05)

### Check whether any of these 40 are in the age-associated metabolites
effect_pval_age <- read.csv("data/04-effect_pval_age.csv", stringsAsFactors = FALSE)
intersect(effect_pval_age$metab_id[effect_pval_age$padj < 0.05], EOAD_effect_pval_df$metab_id[EOAD_effect_pval_df$padj < 0.05])
intersect(effect_pval_age$metab_name[effect_pval_age$padj < 0.05], EOAD_effect_pval_df$metab_name[EOAD_effect_pval_df$padj < 0.05])

### Check whether any of these 40 are in ADADvsCO
effect_pval_all <- read.csv("data/04-effect_pval_adad_ca_trem2.csv", stringsAsFactors = FALSE)
intersect(effect_pval_all$CHEMICAL.ID[effect_pval_all$ADADvsCO_padj < 0.05], EOAD_effect_pval_df$metab_id[EOAD_effect_pval_df$padj < 0.05])
intersect(effect_pval_all$BIOCHEMICAL[effect_pval_all$ADADvsCO_padj < 0.05], EOAD_effect_pval_df$metab_name[EOAD_effect_pval_df$padj < 0.05])

# write.csv(EOAD_effect_pval_df, "data/18-EOAD_effect_pval_df_noage.csv", row.names = FALSE)

## Linear regression for LOAD vs CO -- correcting for age (only one significant regardless of correcting for age or not)
model_data_LOAD <- 
  model_data_duration %>% 
  filter((Status == "CA" & AAO >= 65) | Status == "CO") %>%
  mutate(Status = factor(Status, levels = c("CO", "CA")))

LOAD_effect_pval_df <- ldply(colnames(model_data_LOAD[,-1:-ncol(covariates)]), 
                             get_effect_pval,
                             model_data = model_data_LOAD, 
                             status = "CA",
                             EOAD = FALSE)
LOAD_effect_pval_df$padj <- p.adjust(LOAD_effect_pval_df$pval, method = "BH")
if (all(LOAD_effect_pval_df$metab_id == metab_meta$CHEMICAL.ID)) LOAD_effect_pval_df$metab_name <- metab_meta$BIOCHEMICAL
sum(LOAD_effect_pval_df$padj < 0.05)


### Check differences in age at death and disease duration between the two groups
covariates_filtered <- filter(covariates, Status %in% c("CA", "CO"))
covariates_filtered$AAO[covariates_filtered$Status == "CO"] <- NA

# Create a new status column specifying EOAD and LOAD
covariates_filtered$Status2 <- as.character(covariates_filtered$Status)
covariates_filtered$Status2[covariates_filtered$AAO < 65] <- "EOAD"
covariates_filtered$Status2[covariates_filtered$AAO >= 65] <- "LOAD"

# Add duration
covariates_filtered$duration <- covariates_filtered$age_death - covariates_filtered$AAO
table(covariates_filtered$Status2)

# Isolate only EOAD and LOAD for testing
covariates_for_test <- covariates_filtered %>% filter(Status2 %in% c("EOAD", "LOAD"))

## Check age at death difference
summary(covariates_filtered$age_death[covariates_filtered$Status2 == "EOAD"])
summary(covariates_filtered$age_death[covariates_filtered$Status2 == "LOAD"])
sd(covariates_filtered$age_death[covariates_filtered$Status2 == "EOAD"])
sd(covariates_filtered$age_death[covariates_filtered$Status2 == "LOAD"])

t.test(age_death ~ Status2, data =covariates_for_test)


## Check duration difference
summary(covariates_filtered$duration[covariates_filtered$Status2 == "EOAD"])
summary(covariates_filtered$duration[covariates_filtered$Status2 == "LOAD"])
sd(covariates_filtered$duration[covariates_filtered$Status2 == "EOAD"])
sd(covariates_filtered$duration[covariates_filtered$Status2 == "LOAD"])

t.test(duration ~ Status2, data =covariates_for_test)

## Check differences in CDR
t.test(CDR ~ Status2, data = covariates_for_test)

# Braak tau
t.test(BraakTau ~ Status2, data = covariates_for_test)



