library(dplyr)

setwd("~/ROSMAP_metabolon")

metab_data <- read.csv("ROSMAP_Metabolon_HD4_Brain514_assay_data.csv", stringsAsFactors = FALSE, check.names = FALSE)
pheno <- read.csv("ROSMAP_clinical.csv", stringsAsFactors = FALSE)

# Remove NIST samples, not sure what these are for now
metab_data <- metab_data[-grep("NIST", metab_data$individualID),]

# No replicates
any(duplicated(metab_data$individualID))

# Put metabolites and pheno in same order
row.names(pheno) <- pheno$individualID
pheno <- pheno[metab_data$individualID,]
all(row.names(pheno) == metab_data$individualID)

# Recode sex
pheno$msex <- 
  pheno$msex %>%
  recode(
    "1" = "Male",
    "0" = "Female"
  ) %>%
  as.factor()

# Recode diagnosis
pheno$Status <-
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

# Recode age and set as numeric
pheno$age_death <- 
  pheno$age_death %>%
  as.character() %>%
  recode(
    "90+" = "90"
  ) %>%
  as.numeric()

pheno$age_at_visit_max <- 
  pheno$age_at_visit_max %>%
  as.character() %>%
  recode(
    "90+" = "90"
  ) %>%
  as.numeric()

pheno$apoe_status <- 
  pheno$apoe_genotype %>%
  recode(
    "22" = "negative",
    "23" = "negative",
    "24" = "positive",
    "33" = "negative",
    "34" = "positive",
    "44" = "positive"
  )

pheno <- pheno[!is.na(pheno$Status),]
row.names(metab_data) <- metab_data$individualID
metab_data <- metab_data[pheno$individualID,]

# saveRDS(pheno, "data/13-ROSMAP_pheno_initial_clean.rds")
# saveRDS(metab_data, "data/13-ROSMAP_metab_inital_clean.rds")
