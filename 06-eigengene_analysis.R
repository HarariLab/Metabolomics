library(dplyr)
library(FactoMineR)

# 17 metabolites were significant after correction in ADAD vs CO, and nominally
# significant in AD vs CO and TREM2 vs CO
# Serotonin was part of the 17 but was removed because of high missingness in the ADAD group
metabs_16 <- readRDS("data/05-metabs_16.rds")

# Get scaled, imputed data for PCA
ScaledImp_avgreps <- read.csv("data/03-ScaledImp_avgreps.csv", stringsAsFactors = FALSE, row.names = 1)
ScaledImp_pca <- as.data.frame(t(ScaledImp_avgreps)) %>% select(all_of(as.character(metabs_16)))

# saveRDS(ScaledImp_pca, "data/05-ScaledImp_pca_16.rds")
model_data_all <- readRDS("data/04-model_data_all.rds")
CA_CO_TREM2_ADAD <- row.names(model_data_all)[model_data_all$Status %in% c("CO", "CA", "ADAD", "TREM2")]
ScaledImp_pca_CA_CO_TREM2_ADAD <- ScaledImp_pca[row.names(ScaledImp_pca) %in% CA_CO_TREM2_ADAD,]

# Get first PC
pca_res <- PCA(ScaledImp_pca_CA_CO_TREM2_ADAD, ncp=10, graph = FALSE)
PC1 <- pca_res$ind$coord[,1]
PC1_df <- data.frame(TubeBarcode = names(PC1), PC1 = PC1)

# saveRDS(PC1_df, "data/06-eigengene_df_16.rds")

AAO_df <- readRDS("data/00-age_at_onset_df.rds")

# Distinguish between age at death and age at onset
model_data_duration <- model_data_all %>% rename(age_death = Age)

# Make tubebarcode a variable to join
model_data_duration$TubeBarcode <- row.names(model_data_duration)

# Add age at onset to model data
model_data_duration <- 
  AAO_df %>%
  right_join(model_data_duration)
# Calculate disease duration
model_data_duration$duration <- model_data_duration$age_death - model_data_duration$AAO
# Rearrange data table
model_data_duration <- model_data_duration[,c(1:13, 638, 14:637)]

# saveRDS(model_data_duration, "data/06-model_data_duration.rds")

model_data_duration_small <- 
  model_data_duration %>%
  select(
    TubeBarcode,
    Status,
    duration,
    AAO,
    age_death,
    Sex,
    PMI
  )

# Add PC1 to model_data
model_data_duration_PC1 <- 
  PC1_df %>%
  inner_join(model_data_duration_small)

# saveRDS(model_data_duration_PC1, "data/06-model_data_duration_PC1_16.rds")

model_data_duration_PC1_CA <- model_data_duration_PC1  %>%
  filter(Status == "CA")

mod1 <- glm(PC1 ~ duration + age_death + Sex + PMI, data = model_data_duration_PC1_CA, family = gaussian)
summary(mod1)

# Try testing age at onset as well -- if corrected by AAD it is the same as duration
mod_aao <- glm(PC1~ AAO + Sex + PMI, data = model_data_duration_PC1_CA, family = gaussian)
summary(mod_aao)


## Regression for PC1 ~ CDR for 16 metabs
model_data_CDR_small <- 
  model_data_all %>%
  select(
    Status,
    CDR,
    Age,
    Sex,
    PMI
  ) %>%
  filter(
    Status %in% c("ADAD", "CA", "CO", "TREM2")
  )
model_data_CDR_small$TubeBarcode <- row.names(model_data_CDR_small)

model_data_PC1_CDR <-
  inner_join(
    model_data_CDR_small,
    PC1_df
  )

# Test assoc with CDR for all statuses together, with and without age at death as
# a covariate.
mod_CDR_PC1_withage <- glm(PC1 ~ CDR + Age + Sex + PMI, data = model_data_PC1_CDR, family = gaussian)
summary(mod_CDR_PC1_withage)


## Regression for PC1 ~ BraakTau for 17 metabs
model_data_tau_small <- 
  model_data_all %>%
  select(
    Status,
    BraakTau,
    Age,
    Sex,
    PMI
  ) %>%
  filter(
    Status %in% c("ADAD", "CA", "CO", "TREM2")
  )
model_data_tau_small$TubeBarcode <- row.names(model_data_tau_small)

model_data_PC1_tau <-
  inner_join(
    model_data_tau_small,
    PC1_df
  )

# All statuses together
mod_tau_PC1_withage <- glm(PC1 ~ BraakTau + Age + Sex + PMI, data = model_data_PC1_tau, family = gaussian)
summary(mod_tau_PC1_withage)

## Regression for PC1 ~ BraakAbeta for 16 metabs
model_data_abeta_small <- 
  model_data_all %>%
  select(
    Status,
    BraakAbeta,
    Age,
    Sex
  ) %>%
  filter(
    Status %in% c("ADAD", "CA", "CO", "TREM2")
  )
model_data_abeta_small$TubeBarcode <- row.names(model_data_abeta_small)

model_data_PC1_abeta <-
  inner_join(
    model_data_abeta_small,
    PC1_df
  )

# All statuses together
mod_abeta_PC1_withage <- glm(PC1 ~ BraakAbeta + Age + Sex, data = model_data_PC1_abeta, family = gaussian)
summary(mod_abeta_PC1_withage)


## Regression for PC1 ~ APOE4 status for 16 metabs
model_data_apoe_small <- 
  model_data_all %>%
  select(
    Status,
    APOE,
    Age,
    Sex,
    PMI
  ) %>%
  filter(
    Status %in% c("ADAD", "CA", "CO", "TREM2"),
    APOE != "0"
  )
model_data_apoe_small$TubeBarcode <- row.names(model_data_apoe_small)

model_data_PC1_apoe <-
  inner_join(
    model_data_apoe_small,
    PC1_df
  )

model_data_PC1_apoe$APOE4_status <-
  recode(
    model_data_PC1_apoe$APOE,
    "32" = "Negative",
    "33" = "Negative",
    "42" = "Positive",
    "43" = "Positive",
    "44" = "Positive"
  )  %>%
  droplevels()

# All statuses together
mod_apoe_PC1_withage <- glm(PC1 ~ APOE4_status + Age + Sex + PMI, data = model_data_PC1_apoe, family = gaussian)
summary(mod_apoe_PC1_withage)


## Remove nicotinamide and try associations again
metabs_16_nonic <- setdiff(metabs_16, 432)

# Get scaled, imputed data for PCA
ScaledImp_pca_nonic <- as.data.frame(t(ScaledImp_avgreps)) %>% select(all_of(as.character(metabs_16_nonic)))

# saveRDS(ScaledImp_pca, "data/05-ScaledImp_pca_16.rds")
CA_CO_TREM2_ADAD <- row.names(model_data_all)[model_data_all$Status %in% c("CO", "CA", "ADAD", "TREM2")]
ScaledImp_pca_CA_CO_TREM2_ADAD_nonic <- ScaledImp_pca_nonic[row.names(ScaledImp_pca_nonic) %in% CA_CO_TREM2_ADAD,]

# Get first PC
pca_res_nonic <- PCA(ScaledImp_pca_CA_CO_TREM2_ADAD_nonic, ncp=10, graph = FALSE)
PC1_nonic <- pca_res_nonic$ind$coord[,1]
PC1_df_nonic <- data.frame(TubeBarcode = names(PC1_nonic), PC1 = PC1_nonic)

## CDR no nicotinamide
model_data_PC1_CDR_nonic <-
  inner_join(
    model_data_CDR_small,
    PC1_df_nonic
  )

mod_CDR_PC1_withage_nonic <- glm(PC1 ~ CDR + Age + Sex + PMI, data = model_data_PC1_CDR_nonic, family = gaussian)
summary(mod_CDR_PC1_withage_nonic)

## Tau no nicotinamide 
model_data_PC1_tau_nonic <-
  inner_join(
    model_data_tau_small,
    PC1_df_nonic
  )

# All statuses together
mod_tau_PC1_withage_nonic <- glm(PC1 ~ BraakTau + Age + Sex + PMI, data = model_data_PC1_tau_nonic, family = gaussian)
summary(mod_tau_PC1_withage_nonic)

## Duration no nicotinamide
# Add PC1 to model_data
model_data_duration_PC1_nonic <- 
  PC1_df_nonic %>%
  inner_join(model_data_duration_small)

model_data_duration_PC1_CA_nonic <- model_data_duration_PC1_nonic  %>%
  filter(Status == "CA")
mod1_nonic <- glm(PC1 ~ duration + age_death + Sex + PMI, data = model_data_duration_PC1_CA_nonic, family = gaussian)
summary(mod1_nonic)
