library(plyr)
library(dplyr)
library(ggplot2)
library(VennDiagram)
setwd("~/ROSMAP_metabolon")

pheno <- readRDS("data/13-ROSMAP_pheno_initial_clean.rds")
metab_data <- readRDS("data/13-ROSMAP_metab_inital_clean.rds")
metab_meta <- read.csv("/40/Public_Data/Metabolomics/01.-Raw_Data/ROSMAP_brain_syn26007830/ROSMAP Metabolon HD4 Data Dictionary.csv", stringsAsFactors = FALSE)

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

# also recover tryptophan betaine, low missingness (~25%) and was sig in sADvsCO
recover <- fisher_linear_res$metab[fisher_linear_res$keep]

#### Remove metabolites with >20% missing
drop <- as.character(metab_meta_missing20$SHORT_NAME[!metab_meta_missing20$SHORT_NAME %in% recover])

metab_data_filtered <- metab_data[,setdiff(colnames(metab_data)[9:ncol(metab_data)], drop)]
metab_meta_filtered <- metab_meta[metab_meta$SHORT_NAME %in% colnames(metab_data_filtered),]

# Log transform readings, median has been scaled to 1 but I don't think it has been log transformed
metab_data_filtered_log <- log10(metab_data_filtered)

brain_IQRs <- apply(metab_data_filtered_log, 2, FUN = function(x) IQR(x, na.rm = TRUE))
brain_quantiles <- apply(metab_data_filtered_log, 2, FUN = function(x) quantile(x, na.rm = TRUE))

# Remove outliers
metab_data_transformed <- metab_data_filtered_log
for(j in 1:ncol(metab_data_transformed)) {
  metab_data_transformed[!is.na(metab_data_transformed[, j]) & metab_data_transformed[, j] > brain_quantiles[4, j] + brain_IQRs[j] * 1.5, j] <- NA # Set everything 1.5xIQR above 75% to NA
  metab_data_transformed[!is.na(metab_data_transformed[, j]) & metab_data_transformed[, j] < brain_quantiles[2, j] - brain_IQRs[j] * 1.5, j] <- NA # Set everything 1.5xIQR below 25% to NA
  # Calculate mean and normalize
  metab_mean <- mean(metab_data_transformed[, j], na.rm=TRUE) # Get mean for metabolite
  metab_data_transformed[, j] <- metab_data_transformed[, j] - metab_mean # Subtract mean from values to make mean 0
}


# Find metabolites with >20% missing
percent_na_transformed <- apply(metab_data_transformed, MARGIN = 2, FUN = function(x) sum(is.na(x))/length(x))
sum(percent_na_transformed > 0.2) # 8 metabolites are now missing >20% of readings
percent_na_transformed[which(percent_na_transformed > 0.2)]

# Subjects missing more than 20%
percent_na_subjects <- apply(metab_data_transformed, MARGIN = 1, FUN = function(x) sum(is.na(x))/length(x))
percent_na_subjects[percent_na_subjects > .2]
keep_subjects <- setdiff(row.names(pheno), c("R7506996", "R2730285")) # Keep R1538032, only 20.17% missing

# Filter to only subjects that have low missingness
pheno_clean <- pheno[keep_subjects,]
metab_data_transformed <- metab_data_transformed[keep_subjects,]

# saveRDS(metab_data_transformed, "data/14-metab_data_transformed.rds")
# saveRDS(pheno_clean, "data/14-ROSMAP_pheno_clean.rds")
# saveRDS(metab_meta_filtered, "data/14-metab_meta_filtered.rds")
