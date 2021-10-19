library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(VennDiagram)
setwd("~/ROSMAP_metabolon")

metab_data <- read.csv("ROSMAP_Metabolon_HD4_Brain514_assay_data.csv", stringsAsFactors = FALSE, check.names = FALSE)
metab_meta <- read.csv("ROSMAP_brain_syn26007830/ROSMAP Metabolon HD4 Data Dictionary.csv", stringsAsFactors = FALSE)
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

# Recode diagnosis (based on consensus clinical diagnosis only)
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

metab_meta$SUPER_PATHWAY[metab_meta$SUPER_PATHWAY == ""] <- "Unidentified"

# Replace metabolite IDs with names
all(colnames(metab_data)[9:1063] == metab_meta$CHEM_ID)
colnames(metab_data)[9:1063] <- metab_meta$SHORT_NAME

# Drop unknown metabolites
metab_data <- metab_data[,-grep("X -", colnames(metab_data))]
metab_meta <- metab_meta[metab_meta$SHORT_NAME %in% colnames(metab_data),]

# Find metabolites with >20% missing
percent_na <- apply(metab_data[,9:ncol(metab_data)], MARGIN = 2, FUN = function(x) sum(is.na(x))/length(x))
sum(percent_na > 0.2)

# Find which metabs we are eliminating
metab_data_missing20 <- metab_data[,which(percent_na >= 0.2) + 8]
metab_meta_missing20 <- metab_meta[metab_meta$SHORT_NAME %in% names(percent_na[percent_na >= 0.2]),]

### Fisher tests to recover metabs
all(row.names(metab_data_missing20) == row.names(pheno))
missing20_status_missing <- apply(metab_data_missing20, 2, FUN = function(x) pheno$Status[is.na(x)])

# Temp contingency table for fisher test
cont_table <- data.frame(matrix(nrow=2, ncol=4))
colnames(cont_table) <- names(table(pheno$Status))
rownames(cont_table) <- c("missing", "not_missing")


fisher_results <- data.frame(matrix(nrow=1, ncol=3)) # Store all results
colnames(fisher_results) <- c("metab", "AD_OR", "AD_pval")

# Run fisher test
for (i in 1:length(missing20_status_missing)){
  # Create contingency table
  cont_table[1,]<- table(missing20_status_missing[i]) # Missing
  cont_table[2,]<- table(pheno$Status) - table(missing20_status_missing[i]) # Not Missing
  fisher_results[i,]<- c(names(missing20_status_missing)[i],
                         fisher.test(cont_table[,c(1,2)], workspace = 2000000, conf.level = 0.95)$estimate,
                         fisher.test(cont_table[,c(1,2)], workspace = 2000000, conf.level = 0.95)$p.value)
}

fisher_results$AD_OR <- as.numeric(fisher_results$AD_OR)
fisher_results$AD_pval <- as.numeric(fisher_results$AD_pval)
fisher_results$AD_padj <- p.adjust(fisher_results$AD_pval, method = "BH")
fisher_results_sig <- fisher_results %>% filter(AD_padj < 0.05)

# Make a data frame with effect and p-value
get_effect_pval <- function(metab_name, model_data, status) {
  readings <- cbind(model_data[,c("Status", "age_death", "pmi", "msex")], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  model <- glm(reading ~ Status + age_death + msex + pmi,
               data = readings, family = gaussian)
  data.frame(
    metab = metab_name,
    effect_ROSMAP_metabolon = as.matrix(summary(model)$coefficients)[paste("Status", status, sep = ""),"Estimate"],
    pval_ROSMAP_metabolon = as.matrix(summary(model)$coefficients)[paste("Status", status, sep = ""),"Pr(>|t|)"]
  )
}

model_data_ad <- inner_join(pheno, metab_data) %>% filter(Status %in% c("AD", "CO")) %>% mutate(Status = factor(Status, levels = c("CO", "AD")))

# Are any of them differentially abundant in the same direction as missingness?
fisher_linear_res <- ldply(fisher_results_sig$metab, get_effect_pval, model_data_ad, "AD")
fisher_linear_res <- inner_join(fisher_results_sig, fisher_linear_res)
fisher_linear_res$keep <- (sign(fisher_linear_res$AD_OR-1) != sign(fisher_linear_res$effect_ROSMAP_metabolon)) &
  fisher_linear_res$pval_ROSMAP_metabolon < 0.05

# Metabolites to recover
recover <- fisher_linear_res$metab[fisher_linear_res$keep]