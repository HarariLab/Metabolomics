library(plyr)
library(dplyr)
library(tidyr)

# Metabolon data with unknown metabolites and metadata removed
Origscale <- read.csv("data/OrigScale_metabolites_known.csv", stringsAsFactors = FALSE, na.strings = c(""))
replicate_pairs <- read.csv("data/replicates.csv", head=T, stringsAsFactors = FALSE)
pheno <- read.csv("data/00-pheno_extract.csv", stringsAsFactors = FALSE)

# Remove individual with mismatched status
drop_ids <- pheno$TubeBarcode[pheno$Final_Status == "Neuro_CO" & pheno$CDR_Expire == 1 | is.na(pheno$Final_Status)]
pheno <- pheno %>% filter(!TubeBarcode %in% drop_ids)
Origscale <- Origscale[,setdiff(colnames(Origscale),  drop_ids)]

# Make sure phenofile and metabolite data are in the same order
ordered_tubes <- data.frame(TubeBarcode = colnames(Origscale[14:484]))
pheno <- inner_join(ordered_tubes, pheno)
pheno$Final_CACO <- as.factor(pheno$Final_CACO)

# Get number of missing values for each metabolite
missingvals <- apply(Origscale[,-(1:13)], MARGIN = 1, FUN = function(x) sum(is.na(x)))
Origscale_no_missing <- Origscale[missingvals == 0,]

# Get correlation for log transformed values with high CV
Origscale_no_missing_log <- log10(Origscale_no_missing[,-1:-13])

# Get coefficient of variation for each metabolite
cv = c()
corr <- c()
for( i in 1:nrow(Origscale_no_missing_log)) {
  analyte = t(Origscale_no_missing_log[i,]) # transpose 
  cv = c( cv, sd(analyte)/ mean(analyte) )
}
hist( cv, breaks=50)
with.cv = cv>=0.03& cv<=0.1 # only get results with cv >=0.3, remove one outlier at right tail
sum( with.cv)
for( i in 1:10) {
  data.rep1 = Origscale_no_missing_log[with.cv,replicate_pairs[i,]$Replicate1]
  data.rep2 = Origscale_no_missing_log[with.cv,replicate_pairs[i,]$Replicate2]
  head( data.rep1)
  head( data.rep2)
  c = cor.test(data.rep1, data.rep2)$estimate
  corr[i]<- c
  m = sprintf(" LOG Pair %i: %s and  %s  corr = %0.3f", i, replicate_pairs[i,]$Replicate1, replicate_pairs[i,]$Replicate2, c)
  plot( data.rep1, data.rep2, main= m )
}

### Average replicates
# empty data frame to store replicate info
avg_reps <- as.data.frame(matrix(nrow=nrow(Origscale), ncol=13, NA))
colnames(avg_reps) <- replicate_pairs[,2]

# counters to keep track of how many values are kept/dropped
counter_keep <- 0
counter_NA <- 0

# initialize empty vectors to get metabolite names of values dropped because in middle of reading distribution
metab_drop <- data.frame()
metab_keep <- data.frame()

for (j in 1:13){ # this loops across replicate values (1:13)
  set2 <- data.frame(Rep1 = Origscale[,replicate_pairs[j,2]], # extract replicates for each metabolite
                     Rep2 = Origscale[,replicate_pairs[j,3]]) # extract 1 pair of subjects, all metabs
  for (i in 1:nrow(set2)) { # this loops across metabolites
    couple <- c(set2[i,1], set2[i,2]) # vector of replicate values for 1 metabolite
    test <- sum(!is.na(couple)) # 2 = both values present, 1 = one value present, 0 = both values missing
    if (test == 2) { # average values & keep average if both values present
      avg_reps[i,j] <- mean(couple)
    } 
    else if (test == 1) { # if one value present, check a few conditions to keep or drop
      distrib <- quantile(Origscale[i,14:484], c(0.1,0.9), na.rm = TRUE)
      val <- couple[which(!is.na(couple))]
      if (val < distrib[1]) { # if second value is <10th percentile, keep
        avg_reps[i,j] <- val # keep value
        counter_keep <- counter_keep + 1
        metab_keep <- rbind(metab_keep, Origscale[i,]) # get metab name kept
      } 
      else if (val > distrib[1]) {#if second value is above 10th percentile, drop
        avg_reps[i,j] <- NA # drop both values
        counter_NA <- counter_NA + 1
        metab_drop <- rbind(metab_drop, Origscale[i,]) # get metab name dropped
      } 
      else print("error")
    } 
    else if (test == 0){ # if both values missing, make NA
      avg_reps[i,j] <- NA
    } 
    else print("error")
  }
}

Origscale_noreps1 <- Origscale[,!names(Origscale) %in% replicate_pairs[,2]]
Origscale_noreps2 <- Origscale_noreps1[,!names(Origscale_noreps1) %in% replicate_pairs[,3]]

Origscale_avgreps <- cbind(Origscale_noreps2, avg_reps)


# Count missing values by metabolite, transform to percentage, create df
pct_missing <- apply(Origscale_avgreps[,-(1:13)], MARGIN = 1, FUN = function(x) sum(is.na(x))/460)

# Drop irrelevent columns
drop_fields <- c("PATHWAY.SORTORDER","PLATFORM", "COMP.ID", "RI", "MASS", "PUBCHEM", "CAS", "KEGG", "Group.HMDB")
missingvals_full_df <- Origscale_avgreps[,!(names(Origscale_avgreps) %in% drop_fields)]
missingvals_full_df$pct_missing <- pct_missing

missing20_full <- missingvals_full_df %>% filter(pct_missing >= 0.2)

###### Get pheno for only averaged values
pheno_avg <- pheno[pheno$TubeBarcode %in% colnames(Origscale_avgreps[14:471]),]
ordered_tubes <- as.data.frame(colnames(Origscale_avgreps[14:471]))
colnames(ordered_tubes) <- "TubeBarcode"
pheno_avg <- inner_join(ordered_tubes, pheno_avg)

all(pheno_avg$TubeBarcode == colnames(Origscale_avgreps[,-c(1:13)]))

missing20_status_missing <- apply(missing20_full[,5:463], 1, FUN = function(x) pheno_avg$Final_CACO[is.na(x)])
names(missing20_status_missing) <- missing20_full$BIOCHEMICAL

## Fisher Tests
# Temp contingency table for fisher test
cont_table <- data.frame(matrix(nrow=2, ncol=6))
colnames(cont_table) <- names(table(pheno_avg$Final_CACO))
rownames(cont_table) <- c("Missing", "Not.Missing")


fisher_results_caco <- data.frame(matrix(nrow=1, ncol=13)) # Store all results
fisher_results_caco_nottested <- data.frame(matrix(nrow=1,ncol=2)) # metabs not tested
colnames(fisher_results_caco) <- c("metab_name", "ADAD_OR", "ADAD_fisher_pval", 
                                   "CA_OR", "CA_fisher_pval", "FTD_OR", "FTD_fisher_pval", 
                                   "OT_OR", "OT_fisher_pval", "Presym_OR", "Presym_fisher_pval",
                                   "ADADvsCA_OR", "ADADvsCA_fisher_pval")
colnames(fisher_results_caco_nottested)<- c("position", "metab_name")

# Run fisher test
for (i in 1:length(missing20_status_missing)){
  # Create contingency table
  cont_table[1,]<- table(missing20_status_missing[i]) # Missing
  cont_table[2,]<- table(pheno_avg$Final_CACO) - table(missing20_status_missing[i]) # Not Missing
  fisher_results_caco[i,]<- c(missing20_full$BIOCHEMICAL[i],
                              fisher.test(cont_table[,c(1,3)], workspace = 2000000, conf.level = 0.95)$estimate,
                              fisher.test(cont_table[,c(1,3)], workspace = 2000000, conf.level = 0.95)$p.value,
                              fisher.test(cont_table[,c(2,3)], workspace = 2000000, conf.level = 0.95)$estimate,
                              fisher.test(cont_table[,c(2,3)], workspace = 2000000, conf.level = 0.95)$p.value,
                              fisher.test(cont_table[,c(4,3)], workspace = 2000000, conf.level = 0.95)$estimate,
                              fisher.test(cont_table[,c(4,3)], workspace = 2000000, conf.level = 0.95)$p.value,
                              fisher.test(cont_table[,c(5,3)], workspace = 2000000, conf.level = 0.95)$estimate,
                              fisher.test(cont_table[,c(5,3)], workspace = 2000000, conf.level = 0.95)$p.value,
                              fisher.test(cont_table[,c(6,3)], workspace = 2000000, conf.level = 0.95)$estimate,
                              fisher.test(cont_table[,c(6,3)], workspace = 2000000, conf.level = 0.95)$p.value,
                              fisher.test(cont_table[,c(1,2)], workspace = 2000000, conf.level = 0.95)$estimate,
                              fisher.test(cont_table[,c(1,2)], workspace = 2000000, conf.level = 0.95)$p.value)
}

# Remove NAs from results
fisher_results_caco <- na.omit(fisher_results_caco)
fisher_results_caco[,2:13] <- lapply(fisher_results_caco[,2:13], as.numeric)


fisher_results_sig <- fisher_results_caco %>% filter(ADAD_fisher_pval < 0.05 | CA_fisher_pval < 0.05 | ADADvsCA_fisher_pval < 0.05)

get_effect_pval <- function(metab_name, model_data, status, ADAD = FALSE) {
  readings <- cbind(model_data[,1:18], model_data[,metab_name])
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),]
  if (!all(table(readings$Final_CACO, readings$reading)[1,] == 0) & !all(table(readings$Final_CACO, readings$reading)[2,] == 0) & nrow(table(readings$SEX, readings$reading)) == 2) {
    
    if(ADAD == TRUE) {
      model <- lm(reading ~ Final_CACO + SEX,
                  data = readings)
    }
    else {
      model <- lm(reading ~ Final_CACO + SEX + AAD,
                  data = readings)
    }
    fstat <- summary(model)$fstatistic
    fpval <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
    if (!is.na(fpval) & fpval < 0.05) {
      
      df <- data.frame(
        metab_name = metab_name,
        effect = as.matrix(summary(model)$coefficients)[paste("Final_CACO", status, sep = ""),"Estimate"],
        pval = as.matrix(summary(model)$coefficients)[paste("Final_CACO", status, sep = ""),"Pr(>|t|)"]
      )
      colnames(df)[2:3] <- c(paste0(status, "_effect"), paste0(status, "_pval"))
      df
    }
    else {
      df <- data.frame(metab_name = metab_name, effect = NA, pval = NA)
      colnames(df)[2:3] <- c(paste0(status, "_effect"), paste0(status, "_pval"))
      df
    }
  }
  else {
    df <- data.frame(metab_name = metab_name, effect = NA, pval = NA)
    colnames(df)[2:3] <- c(paste0(status, "_effect"), paste0(status, "_pval"))
    df
  }
}

row.names(Origscale_avgreps) <- Origscale_avgreps$BIOCHEMICAL
metab_data <- t(Origscale_avgreps[,-c(1:13)]) %>% as.data.frame
metab_data$TubeBarcode <- row.names(metab_data)

model_data <- inner_join(pheno, metab_data) %>% mutate(Final_CACO = factor(Final_CACO, levels = c("CO", "ADAD", "CA", "FTD", "OT", "Presymptomatic")))

### Are any of them differentially abundant in the same direction as missingness?

# CA vs CO
model_data_CA <- 
  model_data %>%
  filter(Final_CACO %in% c("CA", "CO")) %>%
  droplevels()
model_data_CA$Final_CACO <- relevel(model_data_CA$Final_CACO, ref = "CO")

fisher_linear_res_CA <- ldply(fisher_results_sig$metab_name, get_effect_pval, model_data_CA, "CA")

# FTD vs CO
model_data_FTD <- 
  model_data %>%
  filter(Final_CACO %in% c("FTD", "CO")) %>%
  droplevels()
model_data_FTD$Final_CACO <- relevel(model_data_FTD$Final_CACO, ref = "CO")

fisher_linear_res_FTD <- ldply(fisher_results_sig$metab_name, get_effect_pval, model_data_FTD, "FTD")

# OT vs CO
model_data_OT <- 
  model_data %>%
  filter(Final_CACO %in% c("OT", "CO")) %>%
  droplevels()
model_data_OT$Final_CACO <- relevel(model_data_OT$Final_CACO, ref = "CO")

fisher_linear_res_OT <- ldply(fisher_results_sig$metab_name, get_effect_pval, model_data_OT, "OT")

# Presym vs CO
model_data_Presym <- 
  model_data %>%
  filter(Final_CACO %in% c("Presymptomatic", "CO")) %>%
  droplevels()
model_data_Presym$Final_CACO <- relevel(model_data_Presym$Final_CACO, ref = "CO")

fisher_linear_res_Presym <- ldply(fisher_results_sig$metab_name, get_effect_pval, model_data_Presym, "Presymptomatic")


# ADAD vs CO
model_data_ADAD <- 
  model_data %>%
  filter(Final_CACO %in% c("ADAD", "CO")) %>%
  droplevels()
model_data_ADAD$Final_CACO <- relevel(model_data_ADAD$Final_CACO, ref = "CO")

fisher_linear_res_ADAD <- ldply(fisher_results_sig$metab_name, get_effect_pval, model_data_ADAD, "ADAD", ADAD = TRUE)

# ADAD vs CA
model_data_ADADvsCA <- 
  model_data %>%
  filter(Final_CACO %in% c("ADAD", "CA")) %>%
  droplevels()
model_data_ADADvsCA$Final_CACO <- relevel(model_data_ADADvsCA$Final_CACO, ref = "CA")

fisher_linear_res_ADADvsCA <- ldply(fisher_results_sig$metab_name, get_effect_pval, model_data_ADADvsCA, "ADAD", ADAD = TRUE)
colnames(fisher_linear_res_ADADvsCA)[2:3] <- c("ADADvsCA_effect", "ADADvsCA_pval")

fisher_linear_res_all <-
  fisher_results_sig %>%
  inner_join(fisher_linear_res_CA) %>%
  inner_join(fisher_linear_res_ADAD) %>%
  inner_join(fisher_linear_res_ADADvsCA) %>%
  inner_join(fisher_linear_res_OT) %>%
  inner_join(fisher_linear_res_Presym) %>%
  inner_join(fisher_linear_res_FTD)


fisher_linear_res_all$keep <- 
  ((sign(fisher_linear_res_all$CA_OR-1) != sign(fisher_linear_res_all$CA_effect)) & fisher_linear_res_all$CA_pval < 0.05 & fisher_linear_res_all$CA_fisher_pval < 0.05) |
  ((sign(fisher_linear_res_all$ADAD_OR-1) != sign(fisher_linear_res_all$ADAD_effect)) & fisher_linear_res_all$ADAD_pval < 0.05 & fisher_linear_res_all$ADAD_fisher_pval < 0.05) |
  ((sign(fisher_linear_res_all$ADADvsCA_OR-1) != sign(fisher_linear_res_all$ADADvsCA_effect)) & fisher_linear_res_all$ADADvsCA_pval < 0.05 & fisher_linear_res_all$ADADvsCA_fisher_pval < 0.05) |
  ((sign(fisher_linear_res_all$OT_OR-1) != sign(fisher_linear_res_all$OT_effect)) & fisher_linear_res_all$OT_pval < 0.05 & fisher_linear_res_all$OT_fisher_pval < 0.05) |
  ((sign(fisher_linear_res_all$Presym_OR-1) != sign(fisher_linear_res_all$Presymptomatic_effect)) & fisher_linear_res_all$Presymptomatic_pval < 0.05 & fisher_linear_res_all$Presym_fisher_pval < 0.05) |
  ((sign(fisher_linear_res_all$FTD_OR-1) != sign(fisher_linear_res_all$FTD_effect)) & fisher_linear_res_all$FTD_pval < 0.05 & fisher_linear_res_all$FTD_fisher_pval < 0.05)

## Remove metabolites with high missingness
fisher_linear_res_all$metab_name[fisher_linear_res_all$keep == TRUE & !is.na(fisher_linear_res_all$keep)]

# Restore metabolites identified in Fisher tests except pregnenetriol sulfate,
# It has too many NAs (>80%) and therefore the IQR adjustment fails.
restore <- c("3-methyl-2-oxobutyrate","4-hydroxyphenylpyruvate","acetylcholine",
             "androsterone sulfate","cysteinylglycine disulfide*", 
             "gamma-glutamyl-epsilon-lysine", "gamma-glutamylphenylalanine",
             "pregnenediol sulfate (C21H34O5S)*", "serotonin","tryptophan betaine")

drop <- as.character(missing20_full$BIOCHEMICAL[!missing20_full$BIOCHEMICAL %in% restore])

Origscale_clean_avgreps <- Origscale_avgreps[!Origscale_avgreps$BIOCHEMICAL %in% drop,]

# Replace missing values for restored metabolites with minimum
for (metab in restore) {
  Origscale_clean_avgreps[Origscale_clean_avgreps$BIOCHEMICAL == metab,-c(1:13)] <-
    as.numeric(Origscale_clean_avgreps[Origscale_clean_avgreps$BIOCHEMICAL == metab,-c(1:13)]) %>%
    replace_na(min(Origscale_clean_avgreps[Origscale_clean_avgreps$BIOCHEMICAL == metab,-c(1:13)], na.rm = TRUE))
}

# Check whether any individuals have missingness >20%
subject_missingness <- apply(Origscale_clean_avgreps[,14:471], MARGIN = 2, FUN = function(x) sum(is.na(x))/627)
any(subject_missingness > 0.2) # FALSE
max(subject_missingness)

quantile(subject_missingness, probs=0.95)

# write.csv(Origscale_clean_avgreps, "data/01-Origscale_clean_avgreps.csv", row.names = FALSE)








