library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

## Load data
# p180 for brains
ROSMAP_p180_brain_FIA <- read.csv("data/ROSMAP_p180_brain_FIA.csv", stringsAsFactors = FALSE) %>%
  filter(projid != "5236 study pool")
ROSMAP_p180_brain_UPLC <- read.csv("data/ROSMAP_p180_brain_UPLC.csv", stringsAsFactors = FALSE) %>%
  filter(projid != "5236 study pool")

# Phenotypic data for ROSMAP
ROSMAP_pheno <- read.csv("data/ROSMAP_pheno.csv")

## Clean up phenotype data
# Get only relevant columns
ROSMAP_pheno_small <- 
  ROSMAP_pheno %>%
  select(
    "projid",
    "study",
    "cogdx", # 1 = control, 2/3 = MCI, 4/5 = AD, 6 = Other
    "ceradsc",
    "msex", # 1 = Male, 0 = Female
    "age_death",
    "apoe_genotype",
    "pmi", # in hours
    "braaksc", # Braak stage (Tau)
    "age_at_visit_max"
  )

# Recode sex
ROSMAP_pheno_small$msex <- 
  ROSMAP_pheno_small$msex %>%
  recode(
    "1" = "Male",
    "0" = "Female"
  ) %>%
  as.factor()

ROSMAP_pheno_small$Status <-
  ROSMAP_pheno_small$ceradsc %>%
  recode(
    "1" = "AD",
    "2" = "PROB",
    "3" = "POSS",
    "4" = "NOT"
  )

ROSMAP_pheno_small$Status[ROSMAP_pheno_small$Status == "POSS" & ROSMAP_pheno_small$braaksc < 4] <- "CO"
ROSMAP_pheno_small$Status[ROSMAP_pheno_small$Status == "NOT" & ROSMAP_pheno_small$braaksc < 4] <- "CO"
ROSMAP_pheno_small$Status[ROSMAP_pheno_small$Status == "PROB" & ROSMAP_pheno_small$braaksc >= 4] <- "AD"
ROSMAP_pheno_small$Status[ROSMAP_pheno_small$Status == "POSS"] <- "Other"
ROSMAP_pheno_small$Status[ROSMAP_pheno_small$Status == "PROB"] <- "Other"
ROSMAP_pheno_small$Status[ROSMAP_pheno_small$Status == "NOT"] <- "Other"

ROSMAP_pheno_small$Status[ROSMAP_pheno_small$Status == "AD" & ROSMAP_pheno_small$cogdx == 1] <- "Other"
ROSMAP_pheno_small$Status[ROSMAP_pheno_small$Status == "CO" & ROSMAP_pheno_small$cogdx %in% c(4,5,6)] <- "Other"

# Recode age and set as numeric
ROSMAP_pheno_small$age_death <- 
  ROSMAP_pheno_small$age_death %>%
  as.character() %>%
  recode(
    "90+" = "90"
  ) %>%
  as.numeric()

ROSMAP_pheno_small$age_at_visit_max <- 
  ROSMAP_pheno_small$age_at_visit_max %>%
  as.character() %>%
  recode(
    "90+" = "90"
  ) %>%
  as.numeric()


## Clean up metabolomic data
# Remove columns if they contain no values
not_all_na <- function(x) any(!is.na(x))

# Brain
ROSMAP_p180_brain_FIA_small <- ROSMAP_p180_brain_FIA[,c(1,12:294)] %>%
  select_if(not_all_na) %>%
  mutate(projid = as.integer(projid))
ROSMAP_p180_brain_UPLC_small <- ROSMAP_p180_brain_UPLC[,c(1,12:95)] %>%
  select_if(not_all_na) %>%
  mutate(projid = as.integer(projid))

model_data_all_brain <- # 172 metabs
  ROSMAP_pheno_small %>%
  right_join(ROSMAP_p180_brain_FIA_small) %>%
  right_join(ROSMAP_p180_brain_UPLC_small) 


## Missingness, IQR, transform - brain
# Get how many readings are missing for each metabolite
missingvals_all_brain <- 
  apply(
    model_data_all_brain[,10:ncol(model_data_all_brain)], 
    MARGIN = 2, 
    FUN = function(x) sum(is.na(x)))
# Convert to percent
pct_missing_all_brain <- missingvals_all_brain/nrow(model_data_all_brain)

# Data frame of all metabolites with their % missing
pct_missing_all_brain_df <- 
  data.frame(
    metabs = colnames(model_data_all_brain)[10:ncol(model_data_all_brain)],
    pct_missing = pct_missing_all_brain
  )
# Get only metabolites missing >=20% values
missing20_all_brain <- pct_missing_all_brain_df[which(pct_missing_all_brain_df$pct_missing >= 0.20),]
nrow(missing20_all_brain) # 15 analytes

# 157 metabolites missing less than 20%
model_data_all_brain_clean <- model_data_all_brain[,!colnames(model_data_all_brain) %in% missing20_all_brain$metabs]


## Log10 transform, IQR adjust, and adjust mean
model_data_all_brain_clean <- 
  model_data_all_brain_clean %>%
  mutate_at(12:ncol(model_data_all_brain_clean), log10)

model_data_all_brain_transformed <- model_data_all_brain_clean
model_data_all_brain_transformed_mean <- vector()
brain_IQRs <- apply(model_data_all_brain_clean[12:ncol(model_data_all_brain_clean)], 2, FUN = function(x) IQR(x, na.rm = TRUE))
brain_quantiles <- apply(model_data_all_brain_clean[12:ncol(model_data_all_brain_clean)], 2, FUN = function(x) quantile(x, na.rm = TRUE))

j <- 12
while(j <= ncol(model_data_all_brain_clean)) {
  model_data_all_brain_transformed[!is.na(model_data_all_brain_clean[, j]) & model_data_all_brain_clean[, j] > brain_quantiles[4, j-11] + brain_IQRs[j-11] * 1.5, j] <- NA # Set everything 1.5xIQR above 75% to NA
  model_data_all_brain_transformed[!is.na(model_data_all_brain_clean[, j]) & model_data_all_brain_clean[, j] < brain_quantiles[2, j-11] - brain_IQRs[j-11] * 1.5, j] <- NA # Set everything 1.5xIQR below 25% to NA
  # Calculate mean and normalize
  model_data_all_brain_transformed_mean[j-11] <- mean(model_data_all_brain_transformed[, j], na.rm=TRUE) # Get mean for metabolite
  model_data_all_brain_transformed[, j] <- model_data_all_brain_transformed[, j] - model_data_all_brain_transformed_mean[j-11] # Subtract mean from values to make mean 0
  j <- j + 1
}



## Define functions for analysis
# Make a data frame with effect and p-value (for brain, uses age at death)
get_effect_pval_brain <- function(metab_name, model_data, status) {
  readings <- cbind(model_data[,1:11], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  if (length(unique(readings$msex)) > 1  & length(unique(readings$Status)) > 1) {
    model <- glm(reading ~ Status + age_death + msex + pmi,
                 data = readings, family = gaussian)
    data.frame(
      metab_name = metab_name,
      effect = as.matrix(summary(model)$coefficients)[paste("Status", status, sep = ""),"Estimate"],
      pval = as.matrix(summary(model)$coefficients)[paste("Status", status, sep = ""),"Pr(>|t|)"]
    )
  }
  else {
    data.frame(
      metab_name = metab_name,
      effect = NA,
      pval = NA
    )
  }
}

################################################################################
# AD vs CO - brain
model_data_ad_brain <- 
  model_data_all_brain_transformed %>%
  filter(Status %in% c("AD", "CO"))
model_data_ad_brain$Status <- factor(model_data_ad_brain$Status, levels = c("CO", "AD"))

# Check for subject missingness
# Remove phenotypic variables
ROSMAP_only_metabs_brain <- model_data_ad_brain[,-1:-9]

# Find individuals who are missing >20% of readings
missingness_subjects_brain <- c()
for (i in 1:nrow(ROSMAP_only_metabs_brain)) {
  missingness_subjects_brain[i] <- sum(is.na(ROSMAP_only_metabs_brain[i,]))/ncol(ROSMAP_only_metabs_brain)
}
which(missingness_subjects_brain > .2) # 36

model_data_ad_brain$projid[23] # 51864085 is missing 22%
missingness_subjects_brain[23]

model_data_ad_brain$projid[28] # 51864085 is missing 22%
missingness_subjects_brain[28]

# model_data_ad_brain$projid[36] # 51864085 is missing 22%
# missingness_subjects_brain[36]

# model_data_ad_brain <- model_data_ad_brain %>% filter(projid != 51864085)
model_data_ad_brain <- model_data_ad_brain %>% filter(!projid %in% c(48748345, 51864085))
table(model_data_ad_brain$Status)
# saveRDS(model_data_ad_brain, "data/09-model_data_ad_ROSMAP_brain.rds")

# Create data frame with effects, pvals, and padj for each metabolite
ad_brain_effect_pval_df <- ldply(colnames(model_data_ad_brain[,-1:-11]), 
                                 get_effect_pval_brain,
                                 model_data = model_data_ad_brain, 
                                 status = "AD")
ad_brain_effect_pval_df$padj <- p.adjust(ad_brain_effect_pval_df$pval, method = "BH")
# saveRDS(ad_brain_effect_pval_df, "data/09-effect_pval_ROSMAP_brain_cerad.rds")



