library(plyr)
library(dplyr)

# Read in metabolomics and meds data
model_data_all <- readRDS("data/04-model_data_all.rds")
pheno <- read.csv("data/03-pheno_avg.csv", row.names = 1, stringsAsFactors = FALSE)
trem2_ind <- pheno$TREM2_all_variants != "0" & !is.na(pheno$TREM2_all_variants)
pheno$Final_CACO[trem2_ind] <- "TREM2"

meds_data <- read.csv("data/t1727_meds_050421.csv", stringsAsFactors = FALSE, na.strings = c("", "NA"))
meds_data$testyear <- sapply(meds_data$testdate, function(x) as.numeric(strsplit(x, "/")[[1]][3]))
meds_data$testdate_dt <- as.Date(meds_data$testdate, format = "%m/%d/%Y")

# Get only IDs, dates, and drugs
meds_data_small <- meds_data[,c(6, 42, 43, 14:41)]

# Get year of death
map_pheno <- read.csv("data/MAP_pheno.csv", stringsAsFactors = FALSE)
map_pheno$full_mapid <- paste0("MAP_", map_pheno$MAP_ID)
map_pheno_yoe <- map_pheno %>% select(full_mapid, Basic_Demographics__YOE, AAO)

pheno <- left_join(pheno, map_pheno_yoe, by = c("MAPID" = "full_mapid"))


# Add MAPIDs to model data
all(row.names(model_data_all) == pheno$TubeBarcode)
model_data_all$MAPID <- pheno$MAPID
model_data_all$AAO <- pheno$AAO
model_data_all$AAO[model_data_all$AAO == "PC"] <- NA
model_data_all$AAO[model_data_all$AAO == ".m"] <- NA
model_data_all$AAO <- as.numeric(model_data_all$AAO)
model_data_all$duration <- model_data_all$Age - as.numeric(model_data_all$AAO)
model_data_all <- model_data_all[,c(636:638,1:635)]

# Get only statuses of interest that also have meds data
model_data_filtered <-
  model_data_all %>% 
  filter(MAPID %in% meds_data$UniquePhenoID & Status %in% c("CA", "CO", "TREM2", "ADAD", "Presymptomatic"))

# Put meds data and pheno data together, get years between last visit and death
meds_pheno_filtered <- 
  inner_join(meds_data_small, pheno, by = c("UniquePhenoID" = "MAPID")) %>%
  filter(Final_CACO %in% c("CA", "CO", "TREM2", "ADAD", "Presymptomatic"))
meds_pheno_filtered$yeardiff <- meds_pheno_filtered$Basic_Demographics__YOE - meds_pheno_filtered$testyear

num_visits <- table(meds_pheno_filtered$UniquePhenoID)
mean(num_visits)
min(num_visits)
max(num_visits)

# filter for only the most recent visit
meds_pheno_filtered_latest  <- 
  meds_pheno_filtered %>%
  group_by(UniquePhenoID) %>%
  slice_max(order_by = testdate_dt)

mean(meds_pheno_filtered_latest$yeardiff)
sd(meds_pheno_filtered_latest$yeardiff)
min(meds_pheno_filtered_latest$yeardiff)
max(meds_pheno_filtered_latest$yeardiff)

# Define functions
# Make a data frame with effect and p-value
get_effect_pval <- function(metab_name, model_data, status, ADAD = FALSE) {
  readings <- cbind(model_data[,1:12], model_data[,metab_name])
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

# Get data for fluoxetine
fluox_names <- c("PROZAC", "FLUOXETI", "fluoxetine")

# Get all rows including fluox
has_any_fluox <- apply(meds_pheno_filtered, MARGIN = 1, function(x) any(x %in% fluox_names))
sum(has_any_fluox)

# Data for only IDs with fluox
any_fluox_data <- meds_pheno_filtered[has_any_fluox,]
length(unique(any_fluox_data$UniquePhenoID))
table(any_fluox_data$Final_CACO[!duplicated(any_fluox_data$UniquePhenoID)])


## Run metabolite models without individuals who have taken fluox in the past 5 years
# AD vs CO
model_data_ad <- 
  model_data_all %>%
  filter(Status %in% c("CA", "CO") & !MAPID %in% any_fluox_data$UniquePhenoID[any_fluox_data$yeardiff <= 5])
model_data_ad$Status <- relevel(model_data_ad$Status, ref = "CO")

table(model_data_ad$Status)

ad_effect_pval_df <- ldply(colnames(model_data_ad[,-1:-14]), 
                           get_effect_pval,
                           model_data = model_data_ad, 
                           status = "CA")
ad_effect_pval_df$padj <- p.adjust(ad_effect_pval_df$pval, method = "BH")

head(ad_effect_pval_df[order(ad_effect_pval_df$padj),])

# TREM2 vs CO
model_data_trem2 <- 
  model_data_all %>%
  filter(Status %in% c("TREM2", "CO") & !MAPID %in% any_fluox_data$UniquePhenoID[any_fluox_data$yeardiff <= 5])
model_data_trem2$Status <- relevel(model_data_trem2$Status, ref = "CO")

table(model_data_trem2$Status)

trem2_effect_pval_df <- ldply(colnames(model_data_trem2[,-1:-14]), 
                              get_effect_pval,
                              model_data = model_data_trem2, 
                              status = "TREM2")
trem2_effect_pval_df$padj <- p.adjust(trem2_effect_pval_df$pval, method = "BH")

head(trem2_effect_pval_df[order(trem2_effect_pval_df$padj),])

# Check association of fluoxetine usage in past 5 years with beta-citrylglutamate levels
model_data_all$fluox_last_5 <- model_data_all$MAPID %in% any_fluox_data$UniquePhenoID[any_fluox_data$yeardiff <= 5]
mod <- glm(fluox_last_5 ~ `100003271`, data = model_data_all[model_data_all$Status == "CA",], family = binomial)
summary(mod)

# Check association of fluox with other variables
model_data_all$has_any_fluox <- model_data_all$MAPID %in% any_fluox_data$UniquePhenoID

# Binomial regression for age at death
mod <- glm(has_any_fluox ~ Age, data = model_data_all[model_data_all$Status == "CA",], family = binomial)
summary(mod)

# Binomial regression for age at onset
mod <- glm(has_any_fluox ~ AAO, data = model_data_all[model_data_all$Status == "CA",], family = binomial)
summary(mod)

# Binomial regression for duration
mod <- glm(has_any_fluox ~ duration, data = model_data_all[model_data_all$Status == "CA",], family = binomial)
summary(mod)

# Vitamins
## Check for any vitamins
has_vitamin <- apply(meds_pheno_filtered, MARGIN = 1, function(x) length(grep("vit|B COMPLE|B-COMPLE|NIACIN", tolower(x))) > 0)
sum(has_vitamin)

vitamin_data <- meds_pheno_filtered[has_vitamin,]
length(unique(vitamin_data$UniquePhenoID))

vitamins_taken <- apply(meds_pheno_filtered, MARGIN = 1, function(x) x[grep("vit", tolower(x))]) %>% unlist() %>% unique()
length(vitamins_taken)



##### Vitamin E

# Only select rows that contain vitamin E, exclude other vitamins
vitE <- setdiff(vitamins_taken, c("VIT B12", "VIT C", "VIT D", "VIT G", "VIT A &", "VIT B", "VIT B CO", "calcium-vitamin D",
                                  "VIT B1", "VIT-B,C,", "VIT B6", "VIT B6 B", "VIT A&D", "VIT A", "VIT B&", "VIT B&C", "VIT B CO",
                                  "VIT B &", "VIT B&C", "VIT D &", "VIT C +", "VIT D +", "VIT A/D", "VIT B 12", "VIT", "NIACIN",
                                  "B COMPLE", "B-COMPLE"))
vitE

# This works for regex purposes
vitE[vitE == "CA++ VIT"] <- " VIT"
vitE[vitE == "CA + VIT"] <- " VIT"
length(vitE)

# Get only people who have taken vitamin E at some point
has_vitE <- apply(meds_pheno_filtered, MARGIN = 1, 
                  function(x) length(grep(paste0(vitE, collapse = "|"), x)) > 0)
vitE_data <- meds_pheno_filtered[has_vitE,]

## Re-run regressions without subjects who took vitamin E in past 5 yrs
# AD vs CO
model_data_ad <- 
  model_data_all %>%
  filter(Status %in% c("CA", "CO") & !MAPID %in% vitE_data$UniquePhenoID[vitE_data$yeardiff <= 5])
model_data_ad$Status <- relevel(model_data_ad$Status, ref = "CO")

table(model_data_ad$Status)

ad_effect_pval_df_novitE <- ldply(colnames(model_data_ad[,-1:-17]), 
                                  get_effect_pval,
                                  model_data = model_data_ad, 
                                  status = "CA")
ad_effect_pval_df_novitE$padj <- p.adjust(ad_effect_pval_df_novitE$pval, method = "BH")

head(ad_effect_pval_df_novitE[order(ad_effect_pval_df_novitE$padj),])


# TREM2 vs CO
model_data_trem2 <- 
  model_data_all %>%
  filter(Status %in% c("TREM2", "CO") & !MAPID %in% vitE_data$UniquePhenoID[vitE_data$yeardiff <= 5])
model_data_trem2$Status <- relevel(model_data_trem2$Status, ref = "CO")

table(model_data_trem2$Status)

trem2_effect_pval_df_novitE <- ldply(colnames(model_data_trem2[,-1:-17]), 
                                     get_effect_pval,
                                     model_data = model_data_trem2, 
                                     status = "TREM2")
trem2_effect_pval_df_novitE$padj <- p.adjust(trem2_effect_pval_df_novitE$pval, method = "BH")

head(trem2_effect_pval_df_novitE[order(trem2_effect_pval_df_novitE$padj),])

# ADAD vs CO
model_data_adad <- 
  model_data_all %>%
  filter(Status %in% c("ADAD", "CO") & !MAPID %in% vitE_data$UniquePhenoID[vitE_data$yeardiff <= 5])
model_data_adad$Status <- relevel(model_data_adad$Status, ref = "CO")

table(model_data_adad$Status)

adad_effect_pval_df_novitE <- ldply(colnames(model_data_adad[,-1:-17]), 
                                    get_effect_pval,
                                    model_data = model_data_adad, 
                                    status = "ADAD",
                                    ADAD = TRUE)
adad_effect_pval_df_novitE$padj <- p.adjust(adad_effect_pval_df_novitE$pval, method = "BH")

head(adad_effect_pval_df_novitE[order(adad_effect_pval_df_novitE$padj),])



##### Retinol (Vitamin A)
vitA <- setdiff(vitamins_taken, c("VIT E", "VIT B12", "VIT C", "VIT D", "VIT G", "VIT B", "VIT B CO", "calcium-vitamin D",
                                  "VIT B1", "VIT-B,C,", "VIT B6", "VIT B6 B", "VIT B&", "VIT B&C", "VIT B CO", "VIT",
                                  "VIT B &", "VIT B&C", "VIT D &", "VIT C +", "VIT D +", "VIT B 12", "vitamin E", "NIACIN",
                                  "B COMPLE", "B-COMPLE"))
vitA

# This works for regex purposes
vitA[vitA == "CA++ VIT"] <- " VIT"
vitA[vitA == "CA + VIT"] <- " VIT"
length(vitA)

# Get only people who have taken vitamin E at some point
has_vitA <- apply(meds_pheno_filtered, MARGIN = 1, 
                  function(x) length(grep(paste0(vitA, collapse = "|"), x)) > 0)
vitA_data <- meds_pheno_filtered[has_vitA,]

## Re-run regressions without subjects who took vitamin A in past 5 yrs
# AD vs CO
model_data_ad <- 
  model_data_all %>%
  filter(Status %in% c("CA", "CO") & !MAPID %in% vitA_data$UniquePhenoID[vitA_data$yeardiff <= 5])
model_data_ad$Status <- relevel(model_data_ad$Status, ref = "CO")

table(model_data_ad$Status)

ad_effect_pval_df_novitA <- ldply(colnames(model_data_ad[,-1:-17]), 
                                  get_effect_pval,
                                  model_data = model_data_ad, 
                                  status = "CA")
ad_effect_pval_df_novitA$padj <- p.adjust(ad_effect_pval_df_novitA$pval, method = "BH")

head(ad_effect_pval_df_novitA[order(ad_effect_pval_df_novitA$padj),])


# TREM2 vs CO
model_data_trem2 <- 
  model_data_all %>%
  filter(Status %in% c("TREM2", "CO") & !MAPID %in% vitA_data$UniquePhenoID[vitA_data$yeardiff <= 5])
model_data_trem2$Status <- relevel(model_data_trem2$Status, ref = "CO")

table(model_data_trem2$Status)

trem2_effect_pval_df_novitA <- ldply(colnames(model_data_trem2[,-1:-17]), 
                                     get_effect_pval,
                                     model_data = model_data_trem2, 
                                     status = "TREM2")
trem2_effect_pval_df_novitA$padj <- p.adjust(trem2_effect_pval_df_novitA$pval, method = "BH")

head(trem2_effect_pval_df_novitA[order(trem2_effect_pval_df_novitA$padj),])

# ADAD vs CO
model_data_adad <- 
  model_data_all %>%
  filter(Status %in% c("ADAD", "CO") & !MAPID %in% vitA_data$UniquePhenoID[vitA_data$yeardiff <= 5])
model_data_adad$Status <- relevel(model_data_adad$Status, ref = "CO")

table(model_data_adad$Status)

adad_effect_pval_df_novitA <- ldply(colnames(model_data_adad[,-1:-17]), 
                                    get_effect_pval,
                                    model_data = model_data_adad, 
                                    status = "ADAD",
                                    ADAD = TRUE)
adad_effect_pval_df_novitA$padj <- p.adjust(adad_effect_pval_df_novitA$pval, method = "BH")

head(adad_effect_pval_df_novitA[order(adad_effect_pval_df_novitA$padj),])



##### Nicotinamide (Vitamin B3)
# Only select rows that contain vitamin B3, exclude other vitamins
vitB <- setdiff(vitamins_taken, c("VIT E", "VIT B12", "VIT C", "VIT D", "VIT G", "VIT A &", "calcium-vitamin D",
                                  "VIT B1", "VIT B6", "VIT A&D", "VIT A", "VIT D &", "VIT C +", "VIT D +", "VIT A/D", "VIT B 12",
                                  "vitamin E", "EYEVITES", "IVITE", "VIT"))


# This works for regex purposes
vitB[vitB == "CA++ VIT"] <- " VIT"
vitB[vitB == "CA + VIT"] <- " VIT"
vitB <- c(vitB, "B COMPLE", "B-COMPLE", "NIACIN")
vitB
length(vitB)

# Get only people who have taken vitamin E at some point
has_vitB <- apply(meds_pheno_filtered, MARGIN = 1, 
                  function(x) length(grep(paste0(vitB, collapse = "|"), x)) > 0)
vitB_data <- meds_pheno_filtered[has_vitB,]

## Re-run regressions without subjects who took vitamin B in past 5 yrs
# AD vs CO
model_data_ad <- 
  model_data_all %>%
  filter(Status %in% c("CA", "CO") & !MAPID %in% vitB_data$UniquePhenoID[vitB_data$yeardiff <= 5])
model_data_ad$Status <- relevel(model_data_ad$Status, ref = "CO")

table(model_data_ad$Status)

ad_effect_pval_df_novitB <- ldply(colnames(model_data_ad[,-1:-17]), 
                                  get_effect_pval,
                                  model_data = model_data_ad, 
                                  status = "CA")
ad_effect_pval_df_novitB$padj <- p.adjust(ad_effect_pval_df_novitB$pval, method = "BH")

head(ad_effect_pval_df_novitB[order(ad_effect_pval_df_novitB$padj),])

# TREM2 vs CO
model_data_trem2 <- 
  model_data_all %>%
  filter(Status %in% c("TREM2", "CO") & !MAPID %in% vitB_data$UniquePhenoID[vitB_data$yeardiff <= 5])
model_data_trem2$Status <- relevel(model_data_trem2$Status, ref = "CO")

table(model_data_trem2$Status)

trem2_effect_pval_df_novitB <- ldply(colnames(model_data_trem2[,-1:-17]), 
                                     get_effect_pval,
                                     model_data = model_data_trem2, 
                                     status = "TREM2")
trem2_effect_pval_df_novitB$padj <- p.adjust(trem2_effect_pval_df_novitB$pval, method = "BH")

head(trem2_effect_pval_df_novitB[order(trem2_effect_pval_df_novitB$padj),])

# ADAD vs CO
model_data_adad <- 
  model_data_all %>%
  filter(Status %in% c("ADAD", "CO") & !MAPID %in% vitB_data$UniquePhenoID[vitB_data$yeardiff <= 5])
model_data_adad$Status <- relevel(model_data_adad$Status, ref = "CO")

table(model_data_adad$Status)

adad_effect_pval_df_novitB <- ldply(colnames(model_data_adad[,-1:-17]), 
                                    get_effect_pval,
                                    model_data = model_data_adad, 
                                    status = "ADAD",
                                    ADAD = TRUE)
adad_effect_pval_df_novitB$padj <- p.adjust(adad_effect_pval_df_novitB$pval, method = "BH")

head(adad_effect_pval_df_novitB[order(adad_effect_pval_df_novitB$padj),])