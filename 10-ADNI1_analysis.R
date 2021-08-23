library(dplyr)
library(tidyr)
library(stringr)
library(ADNIMERGE)

## Load data
# Data from ADNIMERGE
FIA_ADNI1 <- admcdukep180fia
UPLC_ADNI1 <- admcdukep180uplc

# Blood was obtained at baseline visit, viscode == "bl"
pheno <- adnimerge %>%
  filter(RID %in% FIA_ADNI1$RID, COLPROT == "ADNI1", VISCODE == "bl") %>%
  select(
    RID, DX.bl, AGE, PTGENDER, APOE4, EXAMDATE
  ) %>%
  unique() %>%
  mutate(RID = as.numeric(RID))

# Calculate age at blood draw (baseline visit)
demog_small <- ptdemog %>% filter(COLPROT == "ADNI1") %>% select(RID, PTDOB) %>% unique() %>% na.omit()

pheno_demog <- inner_join(pheno, demog_small)
pheno_demog$PTDOB <- as.numeric(pheno_demog$PTDOB)
pheno_demog$AGE <- as.numeric(pheno_demog$AGE)
pheno_demog$exam_year <- sapply(pheno_demog$EXAMDATE, function(x) as.numeric(str_split(x, pattern = "-")[[1]][1]))

pheno_demog$age_at_exam <- pheno_demog$exam_year - pheno_demog$PTDOB

pheno <- pheno_demog %>%
  select(
    RID, DX.bl, age_at_exam, PTGENDER, APOE4
  )

FIA_ADNI1 <- filter(FIA_ADNI1, RID %in% pheno$RID)[,c(2, 10:150)] %>%
  mutate_at(vars(setdiff(colnames(.), "RID")), na_if, "< LOD") %>%
  mutate_if(is.character, as.numeric)
UPLC_ADNI1 <- filter(UPLC_ADNI1, RID %in% pheno$RID)[,c(2, 10:50)] %>%
  mutate_at(vars(setdiff(colnames(.), "RID")), na_if, "< LOD") %>%
  mutate_if(is.character, as.numeric)

# Average the replicates
replicates <- unique(FIA_ADNI1$RID[duplicated(FIA_ADNI1$RID)])

# Average replicates
replicates_only_FIA <- FIA_ADNI1[FIA_ADNI1$RID %in% replicates,]
replicates_only_UPLC <- UPLC_ADNI1[UPLC_ADNI1$RID %in% replicates,]

# Average UPLC
avgs <- c()
uplc_avgs <- matrix(nrow = 0, ncol = 42)
colnames(uplc_avgs) <- colnames(replicates_only_UPLC)


for (rep_id in replicates) {
  for (i in 2:ncol(replicates_only_UPLC)) {
    avgs[1] <- rep_id
    pair <- replicates_only_UPLC[replicates_only_UPLC$RID == rep_id, i]
    if (all(!is.na(pair))) {
      avgs[i] <- mean(pair)
    }
    else if (any(!is.na(pair))) {
      distrib <- quantile(UPLC_ADNI1[,i], c(0.1), na.rm = TRUE)
      if (pair[!is.na(pair)] < distrib) avgs[i] <- pair[!is.na(pair)]
      else avgs[i] <- NA
    }
    else avgs[i] <- NA
  }
  uplc_avgs <- rbind(uplc_avgs, avgs)
}

uplc_noreps <- 
  UPLC_ADNI1[!UPLC_ADNI1$RID %in% replicates,]
UPLC_ADNI1_avgreps <- rbind(uplc_noreps, uplc_avgs)

# Average FIA
avgs <- c()
fia_avgs <- matrix(nrow = 0, ncol = 142)
colnames(fia_avgs) <- colnames(replicates_only_FIA)

for (rep_id in replicates) {
  for (i in 2:ncol(replicates_only_FIA)) {
    avgs[1] <- rep_id
    pair <- replicates_only_FIA[replicates_only_FIA$RID == rep_id, i]
    if (all(!is.na(pair))) {
      avgs[i] <- mean(pair)
    }
    else if (any(!is.na(pair))) {
      distrib <- quantile(FIA_ADNI1[,i], c(0.1), na.rm = TRUE)
      if (pair[!is.na(pair)] < distrib) avgs[i] <- pair[!is.na(pair)]
      else avgs[i] <- NA
    }
    else avgs[i] <- NA
  }
  fia_avgs <- rbind(fia_avgs, avgs)
}

fia_noreps <- 
  FIA_ADNI1[!UPLC_ADNI1$RID %in% replicates,]
FIA_ADNI1_avgreps <- rbind(fia_noreps, fia_avgs)

# Create data frames for models
model_data_all <- # 172 metabs, log transformed and mean adjusted
  pheno %>%
  right_join(FIA_ADNI1_avgreps) %>%
  right_join(UPLC_ADNI1_avgreps)

## Missingness, IQR, transform
# Get how many readings are missing for each metabolite
missingvals_all <- 
  apply(
    model_data_all[,6:ncol(model_data_all)], 
    MARGIN = 2, 
    FUN = function(x) sum(is.na(x)))
# Convert to percent
pct_missing_all <- missingvals_all/nrow(model_data_all)

# Data frame of all metabolites with their % missing
pct_missing_all_df <- 
  data.frame(
    metabs = colnames(model_data_all)[6:ncol(model_data_all)],
    pct_missing = pct_missing_all
  )
# Get only metabolites missing >=20% values
missing20_all <- pct_missing_all_df[which(pct_missing_all_df$pct_missing >= 0.20),]
nrow(missing20_all) # 33 analytes

# 148 metabolites missing less than 20%
model_data_all_clean <- model_data_all[,!colnames(model_data_all) %in% missing20_all$metabs]


## Log10 transform, IQR adjust, and adjust mean
model_data_all_clean <- 
  model_data_all_clean %>%
  mutate_at(6:ncol(model_data_all_clean), log10)

model_data_all_transformed <- model_data_all_clean
model_data_all_transformed_mean <- vector()
IQRs <- apply(model_data_all_clean[6:ncol(model_data_all_clean)], 2, FUN = function(x) IQR(x, na.rm = TRUE))
quantiles <- apply(model_data_all_clean[6:ncol(model_data_all_clean)], 2, FUN = function(x) quantile(x, na.rm = TRUE))

j <- 6
while(j <= ncol(model_data_all_clean)) {
  model_data_all_transformed[!is.na(model_data_all_clean[, j]) & model_data_all_clean[, j] > quantiles[4, j-5] + IQRs[j-5] * 1.5, j] <- NA # Set everything 1.5xIQR above 75% to NA
  model_data_all_transformed[!is.na(model_data_all_clean[, j]) & model_data_all_clean[, j] < quantiles[2, j-5] - IQRs[j-5] * 1.5, j] <- NA # Set everything 1.5xIQR below 25% to NA
  # Calculate mean and normalize
  model_data_all_transformed_mean[j-5] <- mean(model_data_all_transformed[, j], na.rm=TRUE) # Get mean for metabolite
  model_data_all_transformed[, j] <- model_data_all_transformed[, j] - model_data_all_transformed_mean[j-5] # Subtract mean from values to make mean 0
  j <- j + 1
}


## Define functions for analysis
# Make a data frame with effect and p-value
get_effect_pval <- function(metab_name, model_data, status) {
  readings <- cbind(model_data[,1:5], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  if (length(unique(readings$PTGENDER)) > 1  & length(unique(readings$DX.bl)) > 1) {
    model <- glm(reading ~ DX.bl + age_at_exam + PTGENDER,
                 data = readings, family = gaussian)
    data.frame(
      metab_name = metab_name,
      effect = as.matrix(summary(model)$coefficients)[paste("DX.bl", status, sep = ""),"Estimate"],
      pval = as.matrix(summary(model)$coefficients)[paste("DX.bl", status, sep = ""),"Pr(>|t|)"]
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

## AD vs Control
# AD vs CN data
model_data_ad <-
  model_data_all_transformed %>%
  filter(DX.bl %in% c("AD", "CN"))
model_data_ad$DX.bl <- relevel(model_data_ad$DX.bl, ref = "CN")

# Check for subjects with high missingness
ADNI1_only_metabs <- model_data_ad[,-1:-5]

test <- c()
for (i in 1:nrow(ADNI1_only_metabs)) {
  test[i] <- sum(is.na(ADNI1_only_metabs[i,]))/ncol(ADNI1_only_metabs)
}

# Remove subjects missing >20% of readings
model_data_ad$RID[which(test > .2)]
model_data_ad <- model_data_ad[-which(test > .2),]

# saveRDS(model_data_ad, "data/10-model_data_ad_ADNI1.rds")

# Create data frame with effects, pvals, and padj for each metabolite
ad_effect_pval_df <- plyr::ldply(setdiff(colnames(model_data_ad[,-1:-5]), "Spermine"), 
                                 get_effect_pval,
                                 model_data = model_data_ad, 
                                 status = "AD")
ad_effect_pval_df$padj <- p.adjust(ad_effect_pval_df$pval, method = "BH")
# saveRDS(ad_effect_pval_df, "data/10-effect_pval_ADNI1.rds")
