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

# p180 for serum
ROSMAP_p180_serum_FIA <- read.csv("data/ROSMAP_p180_serum_FIA.csv", stringsAsFactors = FALSE) %>%
  filter(projid != "4688 study pool")
ROSMAP_p180_serum_UPLC <- read.csv("data/ROSMAP_p180_serum_UPLC.csv", stringsAsFactors = FALSE) %>%
  filter(projid != "4688 study pool")

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

# Recode diagnosis
ROSMAP_pheno_small$cogdx <-
  ROSMAP_pheno_small$cogdx %>%
  recode(
    "1" = "CO",
    "2" = "MCI",
    "3" = "MCI",
    "4" = "AD",
    "5" = "AD",
    "6" = "Other"
  ) %>%
  as.factor()

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

# Serum
ROSMAP_p180_serum_FIA_small <- ROSMAP_p180_serum_FIA[,c(2,14:154)] %>%
  mutate(projid = as.integer(projid))
ROSMAP_p180_serum_UPLC_small <- ROSMAP_p180_serum_UPLC[,c(2,14:54)] %>%
  select_if(not_all_na) %>%
  mutate(projid = as.integer(projid))

# Average replicates
replicates <- unique(ROSMAP_p180_serum_FIA$projid[duplicated(ROSMAP_p180_serum_FIA$projid)])
replicates_only_serum_FIA <- ROSMAP_p180_serum_FIA_small[ROSMAP_p180_serum_FIA_small$projid %in% replicates,]
replicates_only_serum_UPLC <- ROSMAP_p180_serum_UPLC_small[ROSMAP_p180_serum_UPLC_small$projid %in% replicates,]

# Average UPLC
avgs <- c()
uplc_avgs <- matrix(nrow = 0, ncol = 38)
colnames(uplc_avgs) <- colnames(replicates_only_serum_UPLC)

for (rep_id in replicates) {
  for (i in 2:ncol(replicates_only_serum_UPLC)) {
    avgs[1] <- rep_id
    pair <- replicates_only_serum_UPLC[replicates_only_serum_UPLC$projid == rep_id, i]
    if (sum(!is.na(pair)) > 1) {
      avgs[i] <- mean(pair, na.rm = TRUE)
    }
    else if (any(!is.na(pair))) {
      distrib <- quantile(ROSMAP_p180_serum_FIA_small[,i], c(0.1), na.rm = TRUE)
      if (pair[!is.na(pair)] < distrib) avgs[i] <- pair[!is.na(pair)]
      else avgs[i] <- NA
    }
    else avgs[i] <- NA
  }
  uplc_avgs <- rbind(uplc_avgs, avgs)
}


serum_uplc_noreps <- 
  ROSMAP_p180_serum_UPLC_small[!ROSMAP_p180_serum_UPLC_small$projid %in% replicates,]
ROSMAP_p180_serum_UPLC_small_avgreps <- rbind(serum_uplc_noreps, uplc_avgs)

# Average FIA
avgs <- c()
fia_avgs <- matrix(nrow = 0, ncol = 142)
colnames(fia_avgs) <- colnames(replicates_only_serum_FIA)

for (rep_id in replicates) {
  for (i in 2:ncol(replicates_only_serum_FIA)) {
    avgs[1] <- rep_id
    pair <- replicates_only_serum_FIA[replicates_only_serum_FIA$projid == rep_id, i]
    if (sum(!is.na(pair)) > 1) {
      avgs[i] <- mean(pair, na.rm = TRUE)
    }
    else if (any(!is.na(pair))) {
      distrib <- quantile(ROSMAP_p180_serum_FIA_small[,i], c(0.1), na.rm = TRUE)
      if (pair[!is.na(pair)] < distrib) avgs[i] <- pair[!is.na(pair)]
      else avgs[i] <- NA
    }
    else avgs[i] <- NA
  }
  fia_avgs <- rbind(fia_avgs, avgs)
}

serum_fia_noreps <- 
  ROSMAP_p180_serum_FIA_small[!ROSMAP_p180_serum_FIA_small$projid %in% replicates,]
ROSMAP_p180_serum_FIA_small_avgreps <- rbind(serum_fia_noreps, fia_avgs)

# Create data frames for models
model_data_all_serum <- # 178 metabs
  ROSMAP_pheno_small %>%
  right_join(ROSMAP_p180_serum_FIA_small_avgreps) %>%
  right_join(ROSMAP_p180_serum_UPLC_small_avgreps)



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
  mutate_at(10:ncol(model_data_all_brain_clean), log10)

model_data_all_brain_transformed <- model_data_all_brain_clean
model_data_all_brain_transformed_mean <- vector()
brain_IQRs <- apply(model_data_all_brain_clean[10:ncol(model_data_all_brain_clean)], 2, FUN = function(x) IQR(x, na.rm = TRUE))
brain_quantiles <- apply(model_data_all_brain_clean[10:ncol(model_data_all_brain_clean)], 2, FUN = function(x) quantile(x, na.rm = TRUE))

j <- 10
while(j <= ncol(model_data_all_brain_clean)) {
  model_data_all_brain_transformed[!is.na(model_data_all_brain_clean[, j]) & model_data_all_brain_clean[, j] > brain_quantiles[4, j-9] + brain_IQRs[j-9] * 1.5, j] <- NA # Set everything 1.5xIQR above 75% to NA
  model_data_all_brain_transformed[!is.na(model_data_all_brain_clean[, j]) & model_data_all_brain_clean[, j] < brain_quantiles[2, j-9] - brain_IQRs[j-9] * 1.5, j] <- NA # Set everything 1.5xIQR below 25% to NA
  # Calculate mean and normalize
  model_data_all_brain_transformed_mean[j-9] <- mean(model_data_all_brain_transformed[, j], na.rm=TRUE) # Get mean for metabolite
  model_data_all_brain_transformed[, j] <- model_data_all_brain_transformed[, j] - model_data_all_brain_transformed_mean[j-9] # Subtract mean from values to make mean 0
  j <- j + 1
}


## Missingness, IQR, transform - serum
# Get how many readings are missing for each metabolite
missingvals_all_serum <- 
  apply(
    model_data_all_serum[,10:ncol(model_data_all_serum)], 
    MARGIN = 2, 
    FUN = function(x) sum(is.na(x)))
# Convert to percent
pct_missing_all_serum <- missingvals_all_serum/nrow(model_data_all_serum)

# Data frame of all metabolites with their % missing
pct_missing_all_serum_df <- 
  data.frame(
    metabs = colnames(model_data_all_serum)[10:ncol(model_data_all_serum)],
    pct_missing = pct_missing_all_serum
  )
# Get only metabolites missing >=20% values
missing20_all_serum <- pct_missing_all_serum_df[which(pct_missing_all_serum_df$pct_missing >= 0.20),]
nrow(missing20_all_serum) # 19 analytes

# 159 metabolites missing less than 20%
model_data_all_serum_clean <- model_data_all_serum[,!colnames(model_data_all_serum) %in% missing20_all_serum$metabs]


## Log10 transform, IQR adjust, and adjust mean
model_data_all_serum_clean <- 
  model_data_all_serum_clean %>%
  mutate_at(10:ncol(model_data_all_serum_clean), log10)

model_data_all_serum_transformed <- model_data_all_serum_clean
model_data_all_serum_transformed_mean <- vector()
serum_IQRs <- apply(model_data_all_serum_clean[10:ncol(model_data_all_serum_clean)], 2, FUN = function(x) IQR(x, na.rm = TRUE))
serum_quantiles <- apply(model_data_all_serum_clean[10:ncol(model_data_all_serum_clean)], 2, FUN = function(x) quantile(x, na.rm = TRUE))

j <- 10
while(j <= ncol(model_data_all_serum_clean)) {
  model_data_all_serum_transformed[!is.na(model_data_all_serum_clean[, j]) & model_data_all_serum_clean[, j] > serum_quantiles[4, j-9] + serum_IQRs[j-9] * 1.5, j] <- NA # Set everything 1.5xIQR above 75% to NA
  model_data_all_serum_transformed[!is.na(model_data_all_serum_clean[, j]) & model_data_all_serum_clean[, j] < serum_quantiles[2, j-9] - serum_IQRs[j-9] * 1.5, j] <- NA # Set everything 1.5xIQR below 25% to NA
  # Calculate mean and normalize
  model_data_all_serum_transformed_mean[j-9] <- mean(model_data_all_serum_transformed[, j], na.rm=TRUE) # Get mean for metabolite
  model_data_all_serum_transformed[, j] <- model_data_all_serum_transformed[, j] - model_data_all_serum_transformed_mean[j-9] # Subtract mean from values to make mean 0
  j <- j + 1
}

## Define functions for analysis
# Make a data frame with effect and p-value (for brain, uses age at death)
get_effect_pval_brain <- function(metab_name, model_data, status) {
  readings <- cbind(model_data[,1:9], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  if (length(unique(readings$msex)) > 1  & length(unique(readings$cogdx)) > 1) {
    model <- glm(reading ~ cogdx + age_death + msex + pmi,
                 data = readings, family = gaussian)
    data.frame(
      metab_name = metab_name,
      effect = as.matrix(summary(model)$coefficients)[paste("cogdx", status, sep = ""),"Estimate"],
      pval = as.matrix(summary(model)$coefficients)[paste("cogdx", status, sep = ""),"Pr(>|t|)"]
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

# Regression for serum (uses age at visit max)
get_effect_pval_serum <- function(metab_name, model_data, status) {
  readings <- cbind(model_data[,1:9], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  if (length(unique(readings$msex)) > 1  & length(unique(readings$cogdx)) > 1) {
    model <- glm(reading ~ cogdx + age_at_visit_max + msex,
                 data = readings, family = gaussian)
    data.frame(
      metab_name = metab_name,
      effect = as.matrix(summary(model)$coefficients)[paste("cogdx", status, sep = ""),"Estimate"],
      pval = as.matrix(summary(model)$coefficients)[paste("cogdx", status, sep = ""),"Pr(>|t|)"]
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


# AD vs CO - brain
model_data_ad_brain <- 
  model_data_all_brain_transformed %>%
  filter(cogdx %in% c("AD", "CO"))
model_data_ad_brain$cogdx <- relevel(model_data_ad_brain$cogdx, ref = "CO")

# Check for subject missingness
# Remove phenotypic variables
ROSMAP_only_metabs_brain <- model_data_ad_brain[,-1:-9]

# Find individuals who are missing >20% of readings
missingness_subjects_brain <- c()
for (i in 1:nrow(ROSMAP_only_metabs_brain)) {
  missingness_subjects_brain[i] <- sum(is.na(ROSMAP_only_metabs_brain[i,]))/ncol(ROSMAP_only_metabs_brain)
}
which(missingness_subjects_brain > .2) # 36

model_data_ad_brain$projid[36] # 51864085 is missing 22%
missingness_subjects_brain[36]

model_data_ad_brain <- model_data_ad_brain %>% filter(projid != 51864085)

# saveRDS(model_data_ad_brain, "data/09-model_data_ad_ROSMAP_brain.rds")

# Create data frame with effects, pvals, and padj for each metabolite
ad_brain_effect_pval_df <- ldply(colnames(model_data_ad_brain[,-1:-9]), 
                                 get_effect_pval_brain,
                                 model_data = model_data_ad_brain, 
                                 status = "AD")
ad_brain_effect_pval_df$padj <- p.adjust(ad_brain_effect_pval_df$pval, method = "BH")
# saveRDS(ad_brain_effect_pval_df, "data/09-effect_pval_ROSMAP_brain.rds")


## Serum Analysis
# AD vs CO - serum
model_data_ad_serum <- 
  model_data_all_serum_transformed %>%
  filter(cogdx %in% c("AD", "CO"), !is.na(age_at_visit_max))
model_data_ad_serum$cogdx <- relevel(model_data_ad_serum$cogdx, ref = "CO")

# Check for subject missingness
# Remove phenotypic variables
ROSMAP_only_metabs_serum <- model_data_ad_serum[,-1:-9]

# Find individuals who are missing >20% of readings
missingness_subjects_serum <- c()
for (i in 1:nrow(ROSMAP_only_metabs_serum)) {
  missingness_subjects_serum[i] <- sum(is.na(ROSMAP_only_metabs_serum[i,]))/ncol(ROSMAP_only_metabs_serum)
}
which(missingness_subjects_serum > .2) # None

# saveRDS(model_data_ad_serum, "data/09-model_data_ad_ROSMAP_serum.rds")

# Create data frame with effects, pvals, and padj for each metabolite
ad_serum_effect_pval_df <- ldply(colnames(model_data_ad_serum[,-1:-9]), 
                                 get_effect_pval_serum,
                                 model_data = model_data_ad_serum, 
                                 status = "AD")
ad_serum_effect_pval_df$padj <- p.adjust(ad_serum_effect_pval_df$pval, method = "BH")

# saveRDS(ad_serum_effect_pval_df, "data/09-effect_pval_ROSMAP_serum.rds")
















