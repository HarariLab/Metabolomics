library(plyr)
library(dplyr)

# Read in data
model_data_all_ROSMAP <- readRDS("data/09-model_data_ad_ROSMAP_serum.rds")
model_data_all_ADNI1 <- readRDS("data/10-model_data_ad_ADNI1.rds")
model_data_all_ADNI2GO <- readRDS("data/11-model_data_ad_ADNI2GO.rds")

# Set up data to merge datasets
# Set up ADNI1 data to merge
model_data_ad_ADNI1 <- 
  model_data_all_ADNI1 %>% 
  mutate_if(is.factor, as.character) %>% # Make factors characters for merging
  mutate_if(is.character, as.character) %>% # Remove column labels (there might be a better way to to do this)
  mutate_if(is.numeric, as.numeric) # Same as above

# Fix incorrect column names
colnames(model_data_ad_ADNI1)[colnames(model_data_ad_ADNI1) == "Carmosine"] <- "Carnosine"
colnames(model_data_ad_ADNI1)[colnames(model_data_ad_ADNI1) == "Met.So"] <- "Met.SO"

# Make column names consistent between datasets
colnames(model_data_ad_ADNI1) <- 
  colnames(model_data_ad_ADNI1) %>%
  recode(
    "RID" = "ID",
    "DX.bl" = "Status",
    "PTGENDER" = "Sex",
    "age_at_exam" = "Age"
  ) 


# Set up ADNI2GO data to merge
model_data_ad_ADNI2GO <-
  model_data_all_ADNI2GO %>% 
  mutate_if(is.factor, as.character) %>%
  mutate_if(is.character, as.character) %>%
  mutate_if(is.numeric, as.numeric)

# Make column names consistent between datasets
colnames(model_data_ad_ADNI2GO) <- 
  colnames(model_data_ad_ADNI2GO) %>%
  recode(
    "RID" = "ID",
    "DX.bl" = "Status",
    "PTGENDER" = "Sex",
    "age_at_exam" = "Age"
  )


# Set up ROSMAP data to merge
model_data_ad_ROSMAP <- 
  model_data_all_ROSMAP %>% 
  mutate(
    cogdx = recode(cogdx, "CO" = "CN") # Match control to ADNI
  ) %>%
  mutate_if(is.factor, as.character)

# Make column names consistent between datasets
colnames(model_data_ad_ROSMAP) <- 
  colnames(model_data_ad_ROSMAP) %>%
  recode(
    "projid" = "ID",
    "cogdx" = "Status",
    "msex" = "Sex",
    "age_death" = "Age"
  )

## Add dummy variables for study
# Add dummy variables
model_data_ad_ADNI1$ADNI1 <- 1
model_data_ad_ADNI1$ADNI2GO <- 0

model_data_ad_ADNI2GO$ADNI1 <- 0
model_data_ad_ADNI2GO$ADNI2GO <- 1

model_data_ad_ROSMAP$ADNI1 <- 0
model_data_ad_ROSMAP$ADNI2GO <- 0

## Combine
model_data_ad_all <- 
  bind_rows(
    model_data_ad_ADNI1, 
    model_data_ad_ADNI2GO, 
    model_data_ad_ROSMAP) %>%
  # Get rid of unnecessary variables
  select(!c("APOE4", "study", "apoe_genotype", "pmi", "braaksc")) %>%
  # Recreate factors
  mutate(
    Status = factor(Status, levels = c("CN", "AD")),
    Sex = factor(Sex, levels = c("Female", "Male")),
    ADNI1 = factor(ADNI1, levels = c("0","1")),
    ADNI2GO = factor(ADNI2GO, levels = c("0","1"))
  )


# Define functions for analysis
# Make a data frame with effect and p-value
get_effect_pval <- function(metab_name, model_data, status) {
  readings <- cbind(model_data[,c(1:4, 154, 155)], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  # Only run models if there is more than one level for each factor
  if (length(unique(readings$Sex)) > 1  & 
      length(unique(readings$Status)) > 1 & 
      length(unique(readings$ADNI1)) > 1 & 
      length(unique(readings$ADNI2GO)) > 1)  {
    model <- glm(reading ~ Status + Age + Sex + ADNI1 + ADNI2GO,
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

# Create data frame with effects, pvals, and padj for each metabolite
ad_effect_pval_df <- ldply(colnames(model_data_ad_all)[-c(1:4, 154, 155)],
                           get_effect_pval,
                           model_data = model_data_ad_all, 
                           status = "AD")
ad_effect_pval_df$padj <- p.adjust(ad_effect_pval_df$pval, method = "BH")

# write.csv(ad_effect_pval_df, "data/12-serum_meta-analysis_effect_pval_clean.csv", row.names = FALSE)

