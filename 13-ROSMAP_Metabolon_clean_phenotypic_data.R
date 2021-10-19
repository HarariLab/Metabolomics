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

# Recode diagnosis based on CERAD, Braak, and cogdx
pheno$Status <-
  pheno$ceradsc %>%
  recode(
    "1" = "AD",
    "2" = "PROB",
    "3" = "POSS",
    "4" = "NOT"
  )

pheno$Status[pheno$Status == "POSS" & pheno$braaksc < 4] <- "CO"
pheno$Status[pheno$Status == "NOT" & pheno$braaksc < 4] <- "CO"
pheno$Status[pheno$Status == "PROB" & pheno$braaksc >= 4] <- "AD"
pheno$Status[pheno$Status == "POSS"] <- "Other"
pheno$Status[pheno$Status == "PROB"] <- "Other"
pheno$Status[pheno$Status == "NOT"] <- "Other"

pheno$Status[pheno$Status == "AD" & pheno$cogdx == 1] <- "Other"
pheno$Status[pheno$Status == "CO" & pheno$cogdx %in% c(4,5,6)] <- "Other"


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
