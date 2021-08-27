## Main analysis script

library(plyr)
library(dplyr)

Origscale_clean_transformed <- read.csv("data/02-Origscale_clean_transformed.csv", check.names = FALSE, row.names = 1)
Origscale_clean_avgreps <- read.csv("data/01-Origscale_clean_avgreps.csv")
pheno_avg <- readRDS("data/03-pheno_avg.rds")

metabs <- as.data.frame(Origscale_clean_transformed[row.names(Origscale_clean_transformed) %in% row.names(pheno_avg),])

# Make sure tubes are in the same order
all(row.names(metabs) == row.names(pheno_avg))

## Set up model data

# Make TREM2 their own CACO group
trem2_ind <- pheno_avg$TREM2_all_variants != "0" & !is.na(pheno_avg$TREM2_all_variants)
pheno_avg$Final_CACO[trem2_ind] <- "TREM2"

covariates <- 
  data.frame(
    Status = pheno_avg$Final_CACO,
    Age = pheno_avg$AAD,
    Sex = pheno_avg$SEX,
    PMI = pheno_avg$PMI,
    APOE = as.factor(pheno_avg$APOE),
    BraakAbeta = pheno_avg$BraakAbeta,
    BraakTau = pheno_avg$BraakTau,
    CDR = pheno_avg$CDR_Expire
  )
covariates$BraakAbeta <-
  covariates$BraakAbeta %>%
  as.character() %>%
  recode(
    "A" = "1",
    "B" = "2",
    "C" = "3"
  ) %>%
  na_if("#N/A") %>%
  na_if("Missing") %>%
  na_if("0") %>%
  as.numeric()

covariates$BraakTau <-
  covariates$BraakTau %>%
  na_if(0)

model_data_all <- cbind(covariates, metabs, stringsAsFactors = TRUE)
saveRDS(model_data_all, "data/04-model_data_all.rds")

## Define function for linear regression analysis
# Make a data frame with effect and p-value
get_effect_pval <- function(metab_name, model_data, status, ADAD = FALSE) {
  readings <- cbind(model_data[,1:ncol(covariates)], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  if(ADAD == TRUE) {
    model <- glm(reading ~ Status + Sex + PMI,
                 data = readings, family = gaussian)
  }
  else {
    model <- glm(reading ~ Status + Age + Sex + PMI,
                 data = readings, family = gaussian)
  }
  data.frame(
    metab_name = metab_name,
    effect = as.matrix(summary(model)$coefficients)[paste("Status", status, sep = ""),"Estimate"],
    pval = as.matrix(summary(model)$coefficients)[paste("Status", status, sep = ""),"Pr(>|t|)"]
  )
}

## Linear regression for ADAD vs CO
model_data_adad <- 
  model_data_all %>%
  filter(Status %in% c("ADAD", "CO"))
model_data_adad$Status <- relevel(model_data_adad$Status, ref = "CO")

# Create data frame with effects, pvals, and padj for each metabolite
adad_effect_pval_df <- ldply(colnames(model_data_adad[,-1:-ncol(covariates)]), 
                             get_effect_pval,
                             model_data = model_data_adad, 
                             status = "ADAD", 
                             ADAD = TRUE)
adad_effect_pval_df$padj <- p.adjust(adad_effect_pval_df$pval, method = "BH")

## Linear regression for TREM2 vs CO
model_data_trem2 <- 
  model_data_all %>%
  filter(Status %in% c("CO", "TREM2"))
model_data_trem2$Status <- relevel(model_data_trem2$Status, ref = "CO")

# Create data frame with effects, pvals, and padj for each metabolite
trem2_effect_pval_df <- ldply(colnames(model_data_trem2[,-1:-ncol(covariates)]), 
                              get_effect_pval,
                              model_data = model_data_trem2, 
                              status = "TREM2")
trem2_effect_pval_df$padj <- p.adjust(trem2_effect_pval_df$pval, method = "BH")

## Linear regression for sAD vs CO
model_data_ca <- 
  model_data_all %>%
  filter(Status %in% c("CA", "CO"))
model_data_ca$Status <- relevel(model_data_ca$Status, ref = "CO")

# Create data frame with effects, pvals, and padj for each metabolite
ca_effect_pval_df <- ldply(colnames(model_data_ca[,-1:-ncol(covariates)]), 
                           get_effect_pval,
                           model_data = model_data_ca, 
                           status = "CA")
ca_effect_pval_df$padj <- p.adjust(ca_effect_pval_df$pval, method = "BH")

## Data frame for effect, pval, and padj for each metabolite for ADADvsCO, CAvsCO, TREM2vsCO
# Make distinct column names for each comparison
colnames(adad_effect_pval_df)[2:4] <- c("ADADvsCO_effect", "ADADvsCO_pval", "ADADvsCO_padj")
colnames(ca_effect_pval_df)[2:4] <- c("CAvsCO_effect", "CAvsCO_pval", "CAvsCO_padj")
colnames(trem2_effect_pval_df)[2:4] <- c("TREM2vsCO_effect", "TREM2vsCO_pval", "TREM2vsCO_padj")

# Join into one dataframe
all_effect_pval_df <- 
  adad_effect_pval_df %>%
  inner_join(ca_effect_pval_df) %>%
  inner_join(trem2_effect_pval_df)

# write.csv(all_effect_pval_df, "data/04-effect_pval_adad_ca_trem2.csv", row.names = FALSE)


## Linear regression modeling metabolite values by age, sex, and PMI with only sAD samples
model_data_caonly <- 
  model_data_all %>%
  filter(Status == "CA")

# Make a data frame with effect and p-value
get_effect_pval_age <- function(metab_name, model_data) {
  readings <- cbind(model_data[,1:ncol(covariates)], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  model <- glm(reading ~ Age + Sex,
               data = readings, family = gaussian)
  data.frame(
    metab_name = metab_name,
    effect = as.matrix(summary(model)$coefficients)["Age","Estimate"],
    pval = as.matrix(summary(model)$coefficients)["Age","Pr(>|t|)"]
  )
}

# Create data frame with effects, pvals, and padj for each metabolite
effect_pval_age_ca <- ldply(colnames(model_data_caonly[,-1:-ncol(covariates)]), 
                            get_effect_pval_age,
                            model_data = model_data_caonly)
effect_pval_age_ca$padj <- p.adjust(effect_pval_age_ca$pval, method = "BH")



### Linear regression for APOE
model_data_apoe <-  
  model_data_all %>%
  filter(APOE != "0") %>%
  filter(Status == "CA")
model_data_apoe$APOE4_status <- 
  recode(
    model_data_apoe$APOE,
    "32" = "Negative",
    "33" = "Negative",
    "42" = "Positive",
    "43" = "Positive",
    "44" = "Positive"
  ) %>%
  droplevels() # Drop unused levels ("0")
model_data_apoe <- model_data_apoe[,c(1:8,636,9:635)]
model_data_apoe$APOE4_status <- factor(model_data_apoe$APOE4_status, levels = c("Negative", "Positive"))


# Make a data frame with effect and p-value
get_effect_pval_APOE <- function(metab_name, model_data) {
  readings <- cbind(model_data[,1:(ncol(covariates) + 1)], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  model <- glm(reading ~ APOE4_status + Age + Sex + PMI,
               data = readings, family = gaussian)
  data.frame(
    metab_name = metab_name,
    effect = as.matrix(summary(model)$coefficients)["APOE4_statusPositive","Estimate"],
    pval = as.matrix(summary(model)$coefficients)["APOE4_statusPositive","Pr(>|t|)"]
  )
}

effect_pval_apoe <- ldply(colnames(model_data_apoe[,-1:-(ncol(covariates) + 1)]), 
                            get_effect_pval_APOE,
                            model_data = model_data_apoe)
effect_pval_apoe$padj <- p.adjust(effect_pval_apoe$pval, method = "BH")


### Run regressions again matching for CDR and Braak scores
## Linear regression of ADAD vs CO with only high CDR (CDR = 3)
model_data_ADAD_highCDR <- 
  model_data_all %>%
  filter(Status == "CO" | (Status == "ADAD" & CDR == 3))
model_data_ADAD_highCDR$Status <- relevel(model_data_ADAD_highCDR$Status, ref = "CO")

# Create data frame with effects, pvals, and padj for each metabolite
ADAD_highCDR_effect_pval_df <- ldply(colnames(model_data_ADAD_highCDR[,-c(1:11)]), 
                                     get_effect_pval,
                                     model_data = model_data_ADAD_highCDR, 
                                     status = "ADAD",
                                     ADAD = TRUE)
ADAD_highCDR_effect_pval_df$padj <- p.adjust(ADAD_highCDR_effect_pval_df$pval, method = "BH")



## Linear regression of TREM2 vs CO with only high Abeta (3)
model_data_trem2_highAbeta <- 
  model_data_all %>%
  filter(Status == "CO" | (Status == "TREM2" & BraakAbeta == 3))
model_data_trem2_highAbeta$Status <- relevel(model_data_trem2_highAbeta$Status, ref = "CO")

# Create data frame with effects, pvals, and padj for each metabolite
trem2_highAbeta_effect_pval_df <- ldply(colnames(model_data_trem2_highAbeta[,-c(1:11)]), 
                                        get_effect_pval,
                                        model_data = model_data_trem2_highAbeta, 
                                        status = "TREM2")
trem2_highAbeta_effect_pval_df$padj <- p.adjust(trem2_highAbeta_effect_pval_df$pval, method = "BH")

## Linear regression of TREM2 vs CO with only high Tau (4-6)
model_data_trem2_highTau <- 
  model_data_all %>%
  filter(Status == "CO" | (Status == "TREM2" & BraakTau > 3))
model_data_trem2_highTau$Status <- relevel(model_data_trem2_highTau$Status, ref = "CO")

# Create data frame with effects, pvals, and padj for each metabolite
trem2_highTau_effect_pval_df <- ldply(colnames(model_data_trem2_highTau[,-c(1:11)]), 
                                      get_effect_pval,
                                      model_data = model_data_trem2_highTau, 
                                      status = "TREM2")
trem2_highTau_effect_pval_df$padj <- p.adjust(trem2_highTau_effect_pval_df$pval, method = "BH")

## Linear regression of TREM2 vs CO with only high CDR (CDR = 3)
model_data_trem2_highCDR <- 
  model_data_all %>%
  filter(Status == "CO" | (Status == "TREM2" & CDR == 3))
model_data_trem2_highCDR$Status <- relevel(model_data_trem2_highCDR$Status, ref = "CO")

# Create data frame with effects, pvals, and padj for each metabolite
trem2_highCDR_effect_pval_df <- ldply(colnames(model_data_trem2_highCDR[,-c(1:11)]), 
                                      get_effect_pval,
                                      model_data = model_data_trem2_highCDR, 
                                      status = "TREM2")
trem2_highCDR_effect_pval_df$padj <- p.adjust(trem2_highCDR_effect_pval_df$pval, method = "BH")


## Linear regression of sAD vs CO with only high Abeta
model_data_ca_highAbeta <- 
  model_data_all %>%
  filter(Status == "CO" | (Status == "CA" & BraakAbeta == 3))
model_data_ca_highAbeta$Status <- relevel(model_data_ca_highAbeta$Status, ref = "CO")

# Create data frame with effects, pvals, and padj for each metabolite
ca_highAbeta_effect_pval_df <- ldply(colnames(model_data_ca_highAbeta[,-c(1:11)]), 
                                     get_effect_pval,
                                     model_data = model_data_ca_highAbeta, 
                                     status = "CA")
ca_highAbeta_effect_pval_df$padj <- p.adjust(ca_highAbeta_effect_pval_df$pval, method = "BH")

## Linear regression of sAD vs CO with only high Tau (4-6)
model_data_ca_highTau <- 
  model_data_all %>%
  filter(Status == "CO" | (Status == "CA" & BraakTau > 3))
model_data_ca_highTau$Status <- relevel(model_data_ca_highTau$Status, ref = "CO")

# Create data frame with effects, pvals, and padj for each metabolite
ca_highTau_effect_pval_df <- ldply(colnames(model_data_ca_highTau[,-c(1:11)]), 
                                   get_effect_pval,
                                   model_data = model_data_ca_highTau, 
                                   status = "CA")
ca_highTau_effect_pval_df$padj <- p.adjust(ca_highTau_effect_pval_df$pval, method = "BH")

## Linear regression of sAD vs CO with only high CDR (CDR = 3)
model_data_ca_highCDR <- 
  model_data_all %>%
  filter(Status == "CO" | (Status == "CA" & CDR == 3))
model_data_ca_highCDR$Status <- relevel(model_data_ca_highCDR$Status, ref = "CO")

# Create data frame with effects, pvals, and padj for each metabolite
ca_highCDR_effect_pval_df <- ldply(colnames(model_data_ca_highCDR[,-c(1:11)]), 
                                   get_effect_pval,
                                   model_data = model_data_ca_highCDR, 
                                   status = "CA")
ca_highCDR_effect_pval_df$padj <- p.adjust(ca_highCDR_effect_pval_df$pval, method = "BH")




### Merge matched results
## Abeta
colnames(adad_effect_pval_df) <- c("metab", "ADADvsCO_effect", "ADADvsCO_pval", "ADADvsCO_padj")
colnames(ca_highAbeta_effect_pval_df) <- c("metab", "CAvsCO_effect", "CAvsCO_pval", "CAvsCO_padj")
colnames(trem2_highAbeta_effect_pval_df) <- c("metab", "TREM2vsCO_effect", "TREM2vsCO_pval", "TREM2vsCO_padj")

high_Abeta_effects_pvals <- 
  inner_join(adad_effect_pval_df, ca_highAbeta_effect_pval_df) %>%
  inner_join(trem2_highAbeta_effect_pval_df)
# write.csv(high_Abeta_effects_pvals, "data/04-highAbeta_effect_pvals.csv", row.names = FALSE)

## Tau
colnames(adad_effect_pval_df) <- c("metab", "ADADvsCO_effect", "ADADvsCO_pval", "ADADvsCO_padj")
colnames(ca_highTau_effect_pval_df) <- c("metab", "CAvsCO_effect", "CAvsCO_pval", "CAvsCO_padj")
colnames(trem2_highTau_effect_pval_df) <- c("metab", "TREM2vsCO_effect", "TREM2vsCO_pval", "TREM2vsCO_padj")

high_Tau_effects_pvals <- 
  inner_join(adad_effect_pval_df, ca_highTau_effect_pval_df) %>%
  inner_join(trem2_highTau_effect_pval_df)
# write.csv(high_Tau_effects_pvals, "data/04-highTau_effect_pvals.csv", row.names = FALSE)

## CDR
colnames(ADAD_highCDR_effect_pval_df) <- c("metab", "ADADvsCO_effect", "ADADvsCO_pval", "ADADvsCO_padj")
colnames(ca_highCDR_effect_pval_df) <- c("metab", "CAvsCO_effect", "CAvsCO_pval", "CAvsCO_padj")
colnames(trem2_highCDR_effect_pval_df) <- c("metab", "TREM2vsCO_effect", "TREM2vsCO_pval", "TREM2vsCO_padj")

high_CDR_effects_pvals <- 
  inner_join(ADAD_highCDR_effect_pval_df, ca_highCDR_effect_pval_df) %>%
  inner_join(trem2_highCDR_effect_pval_df)
# write.csv(high_CDR_effects_pvals, "data/04-highCDR_effect_pvals.csv", row.names = FALSE)

















