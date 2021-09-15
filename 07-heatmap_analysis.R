# Create heatmap of eigengene metabolites among participants
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(FactoMineR)
library(ggplot2)
library(ggpubr)

# Get scaled, imputed data to plot heatmap and make sure the rows are in the
# Same order as the model data so they can be merged
ScaledImp_avgreps <- read.csv("data/03-ScaledImp_avgreps.csv", stringsAsFactors = FALSE, row.names = 1)
model_data_all <- readRDS("data/04-model_data_all.rds")
all(row.names(model_data_all) == colnames(ScaledImp_avgreps))

metab_meta <- readRDS("data/04-metab_meta.rds")

# Get only individuals from the four status groups of interest
CA_CO_TREM2_ADAD <- row.names(model_data_all)[model_data_all$Status %in% c("CO", "CA", "ADAD", "TREM2")]
ScaledImp_avgreps <- ScaledImp_avgreps[,which(colnames(ScaledImp_avgreps) %in% CA_CO_TREM2_ADAD),]

# Order the individuals by status group in order to split heatmap
CA <- row.names(model_data_all)[model_data_all$Status == "CA"]
ADAD <- row.names(model_data_all)[model_data_all$Status == "ADAD"]
TREM2 <- row.names(model_data_all)[model_data_all$Status == "TREM2"]
CO <- row.names(model_data_all)[model_data_all$Status == "CO"]
ordered_barcodes <- c(ADAD, TREM2, CA, CO)
ScaledImp_avgreps <- ScaledImp_avgreps[ordered_barcodes]

# 16 metabolites chosen for profile, not
# including serotonin due to high missingness
metabs_16 <- readRDS("data/05-metabs_16.rds")
metab_meta_16 <- metab_meta %>% filter(CHEMICAL.ID %in% metabs_16)

# Get only data for those 16 metabolites
ScaledImp_16 <- ScaledImp_avgreps[row.names(ScaledImp_avgreps) %in% metabs_16,]
all(row.names(ScaledImp_16) == metab_meta_16$CHEMICAL.ID)
row.names(ScaledImp_16) <- metab_meta_16$BIOCHEMICAL
ScaledImp_16_mat <- as.matrix(ScaledImp_16)

# Load age at onset data and order samples
age_at_onset_df <- readRDS("data/00-age_at_onset_df.rds")
row.names(age_at_onset_df) <- age_at_onset_df$TubeBarcode
age_at_onset_df <- age_at_onset_df[ordered_barcodes,]

# Get annotation data
anno <- model_data_all[ordered_barcodes, c("CDR", "BraakTau", "BraakAbeta", "Sex", "Age", "PMI")]

# Add age at onset data to annotation data
all(row.names(anno) == row.names(age_at_onset_df))
anno <- cbind(anno, age_at_onset_df[ordered_barcodes,]$AAO)
colnames(anno)[7] <- "AAO"

# Controls should not have age at onset or duration
anno$AAO[row.names(anno) %in% CO]
anno$Duration <- anno$Age - anno$AAO

anno_small <- anno[,c("CDR","BraakTau")]

# Update ADAD CDR
anno_small$CDR[row.names(anno_small) %in% c("FB06002985", "FB06002999", "FB06002991", "FB06002785", "FB06002677", "FB06002779")] <- 3
anno_small$CDR[row.names(anno_small) == "FB06003001"] <- NA

anno_small$CDR <- factor(anno_small$CDR, levels = c("3", "2", "1", "0.5", "0"))
anno_small$BraakTau <- factor(anno_small$BraakTau, levels = c("6", "5", "4", "3", "2", "1", "0"))

# Specify which individuals correspond to which groups
# for splitting the heatmap
group <- factor(c(rep(c("ADAD", "TREM2", "AD", "CO"), times = c(25, 21, 305, 26))), levels = c("ADAD", "TREM2", "AD", "CO"))

# Specify colors for heatmap
col_fun = colorRamp2(c(-2, 0, 1.5), c("deepskyblue1", "black", "yellow"))

# Specify annotation colors
CDR_colors <- c("0" = "#0D0887FF", "0.5" = "#7E03A8FF", "1" = "#CC4678FF", "2" = "#F89441FF", "3" = "#F0F921FF")
BraakTau_colors <- viridis::plasma(7)
names(BraakTau_colors) <- c("0", "1", "2", "3", "4", "5", "6")

# Draw Heatmap and use clustering algorithm
# to order columns, this is the starting point for
# the later manual ordering
hm <- Heatmap(
  t(scale(t(ScaledImp_16_mat))),
  col = col_fun,
  # cluster_columns = FALSE,
  clustering_method_columns = "ward.D2",
  # show_column_names = FALSE,
  # show_column_dend = FALSE,
  column_names_gp = gpar(fontsize = 5),
  column_dend_height = unit(40, "mm"),
  show_row_dend = FALSE,
  column_split = group,
  cluster_column_slices = FALSE,
  column_title_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 8),
  show_heatmap_legend = FALSE,
  row_order = c(
    "1,5-anhydroglucitol (1,5-AG)",
    "glutamate",
    "gamma-glutamylthreonine",
    "glutarate (C5-DC)",
    "3-hydroxy-2-ethylpropionate",
    "aspartate",
    "alpha-tocopherol",
    "retinol (Vitamin A)",
    "ergothioneine",
    "2-methylcitrate/homocitrate",
    "N-acetylglutamate",
    "CDP-ethanolamine",
    "CDP-choline",
    "glycerophosphoinositol*",
    "nicotinamide",
    "beta-citrylglutamate")
)
hm

#### Include presymptomatic ***This overwrites data frames**

# Get scaled, imputed data to plot heatmap and make sure the rows are in the
# Same order as the model data so they can be merged
ScaledImp_avgreps <- read.csv("data/03-ScaledImp_avgreps.csv", stringsAsFactors = FALSE, row.names = 1)
model_data_all <- readRDS("data/04-model_data_all.rds")
all(row.names(model_data_all) == colnames(ScaledImp_avgreps))

# Get only individuals from the four status groups of interest
CA_CO_TREM2_ADAD_Presym <- row.names(model_data_all)[model_data_all$Status %in% c("CO", "CA", "ADAD", "TREM2", "Presymptomatic")]
ScaledImp_avgreps <- ScaledImp_avgreps[,which(colnames(ScaledImp_avgreps) %in% CA_CO_TREM2_ADAD_Presym),]

# Order the individuals by status group in order to split heatmap
Presymptomatic <- row.names(model_data_all)[model_data_all$Status == "Presymptomatic"]
ordered_barcodes <- c(ADAD, TREM2, CA, Presymptomatic, CO)
ScaledImp_avgreps <- ScaledImp_avgreps[ordered_barcodes]

# Get only data for those 16 metabolites
ScaledImp_16 <- ScaledImp_avgreps[row.names(ScaledImp_avgreps) %in% c(metabs_16),] #, "serotonin"),]
all(row.names(ScaledImp_16) == metab_meta_16$CHEMICAL.ID)
row.names(ScaledImp_16) <- metab_meta_16$BIOCHEMICAL
ScaledImp_16_mat <- as.matrix(ScaledImp_16)

# Set groups to split up heatmap by status
group <- factor(c(rep(c("ADAD", "TREM2", "sAD", "Presymptomatic", "CO"), times = c(25, 21, 305, 15, 26))), 
                levels = c("CO", "Presymptomatic", "sAD", "TREM2", "ADAD"))

# Annotations for heatmap
# Get annotation data
anno <- model_data_all[ordered_barcodes, c("CDR", "BraakTau", "BraakAbeta", "Sex", "Age", "PMI")]

# Add age at onset data to annotation data
# Load age at onset data and order samples
age_at_onset_df <- readRDS("data/00-age_at_onset_df.rds")
row.names(age_at_onset_df) <- age_at_onset_df$TubeBarcode
age_at_onset_df <- age_at_onset_df[ordered_barcodes,]
all(row.names(anno) == row.names(age_at_onset_df))
anno <- cbind(anno, age_at_onset_df[ordered_barcodes,]$AAO)
colnames(anno)[7] <- "AAO"

# Controls should not have age at onset or duration
anno$AAO[row.names(anno) %in% CO]
anno$Duration <- anno$Age - anno$AAO
anno_small <- model_data_all[ordered_barcodes, c("CDR", "BraakTau")]

# Fixing ADAD CDR based on file from Fengxian
anno_small$CDR[row.names(anno_small) %in% c("FB06002985", "FB06002999", "FB06002991", "FB06002785", "FB06002677", "FB06002779")] <- 3
anno_small$CDR[row.names(anno_small) == "FB06003001"] <- NA

# Relevel CDR and BraakTau so that the legends are in the right order
anno_small$CDR <- factor(anno_small$CDR, levels = c("3", "2", "1", "0.5", "0"))
anno_small$BraakTau <- factor(anno_small$BraakTau, levels = c("6", "5", "4", "3", "2", "1", "0"))

hm4 <- Heatmap(
  t(scale(t(ScaledImp_16_mat))),
  col = col_fun,
  # clustering_method_columns = "ward.D2",
  show_column_names = FALSE,
  column_title_rot = 90,
  column_dend_height = unit(40, "mm"),
  show_row_dend = FALSE,
  column_split = group,
  height = unit(8, "cm"),
  cluster_column_slices = FALSE,
  column_title_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(title = "Metabolite Levels", title_position = "leftcenter-rot"),
  row_order = c(
    "1,5-anhydroglucitol (1,5-AG)",
    "glutamate",
    "gamma-glutamylthreonine",
    "glutarate (C5-DC)",
    "3-hydroxy-2-ethylpropionate",
    "aspartate",
    "alpha-tocopherol",
    "retinol (Vitamin A)",
    "ergothioneine",
    "2-methylcitrate/homocitrate",
    "N-acetylglutamate",
    "CDP-choline",
    "CDP-ethanolamine",
    "glycerophosphoinositol*",
    "nicotinamide",
    "beta-citrylglutamate"),
  # Keep original column orders the same and just add presym
  column_order = c(column_order(hm)$CO, 352:366, column_order(hm)$AD, column_order(hm)$TREM2, column_order(hm)$ADAD),
  bottom_annotation = HeatmapAnnotation(df = anno_small,
                                        col = list(
                                          BraakTau = BraakTau_colors,
                                          CDR = CDR_colors
                                        ),
                                        na_col = "white",
                                        annotation_name_gp = gpar(fontsize = 10),
                                        annotation_legend_param = list(title_position = "leftcenter-rot")
  )
)
draw(hm4, heatmap_legend_side = "left", annotation_legend_side = "right")

AD_barcode_ESAD <- colnames(ScaledImp_16_mat)[column_order(hm)$AD[1:30]]
# saveRDS(AD_barcode_ESAD, "data/07-ESAD_barcodes.rds")
dat <- anno[CA,]
dat$ESAD <- row.names(dat) %in% AD_barcode_ESAD

##### CDR
cdr_model <- glm(ESAD ~ CDR + Age + Sex + PMI, data = dat, family = binomial)
summary(cdr_model)
mean(dat$CDR[dat$ESAD])
sd(dat$CDR[dat$ESAD])
mean(dat$CDR[!dat$ESAD])
sd(dat$CDR[!dat$ESAD])

##### Braak Tau
braaktau_model <- glm(ESAD ~ BraakTau + Age + Sex + PMI, data = dat, family = binomial)
summary(braaktau_model)
mean(dat$BraakTau[dat$ESAD], na.rm = TRUE)
sd(dat$BraakTau[dat$ESAD], na.rm = TRUE)
mean(dat$BraakTau[!dat$ESAD], na.rm = TRUE)
sd(dat$BraakTau[!dat$ESAD], na.rm = TRUE)

##### Age at Onset
aao_model <- glm(ESAD ~ AAO + Sex + PMI, data = dat, family = binomial)
summary(aao_model)

##### Duration
duration_model <- glm(ESAD ~ Duration + Age + Sex + PMI, data = dat, family = binomial)
summary(duration_model)

mean(dat$Duration[dat$ESAD], na.rm = TRUE)
sd(dat$Duration[dat$ESAD], na.rm = TRUE)
mean(dat$Duration[!dat$ESAD], na.rm = TRUE)
sd(dat$Duration[!dat$ESAD], na.rm = TRUE)

##### AAD
aad_model <- glm(ESAD ~ Age + Sex + PMI, data = dat, family = binomial)
summary(duration_model)

## Plot eigengene between status groups
model_data_all <- readRDS("data/04-model_data_all.rds")

AD_barcode_sAD <- setdiff(CA, AD_barcode_ESAD)

model_data_all$status2 <- as.character(model_data_all$Status)
model_data_all$status2[row.names(model_data_all) %in% AD_barcode_ESAD] <- "Early Stage AD"
model_data_all$status2[row.names(model_data_all) %in% AD_barcode_sAD] <- "AD"
model_data_all$status2 <- as.factor(model_data_all$status2)

# Get scaled, imputed data for PCA
ScaledImp_avgreps <- read.csv("data/03-ScaledImp_avgreps.csv", stringsAsFactors = FALSE, row.names = 1)
ScaledImp_pca <- as.data.frame(t(ScaledImp_avgreps)) %>% dplyr::select(all_of(as.character(metabs_16)))

# Get tube barcodes in the desired status groups
CA_CO_TREM2_ADAD_Presym <- row.names(model_data_all)[model_data_all$Status %in% c("CO", "CA", "ADAD", "TREM2", "Presymptomatic")]
ScaledImp_pca_withpresym <- ScaledImp_pca[CA_CO_TREM2_ADAD_Presym ,]

pca_res_withpresym <- PCA(ScaledImp_pca_withpresym, ncp=10, graph = FALSE)

## Get first PC
PC1_withpresym <- pca_res_withpresym$ind$coord[,1]
PC1_df_withpresym <- data.frame(TubeBarcode = names(PC1_withpresym), PC1 = PC1_withpresym)

# Join statuses to first PC
status2_df <- data.frame(TubeBarcode = row.names(model_data_all), Status = model_data_all$status2)
PC1_df_withpresym <- inner_join(status2_df, PC1_df_withpresym)

# Relevel to order plots
PC1_df_withpresym$Status <- 
  factor(
    PC1_df_withpresym$Status, 
    levels = c("CO", "Presymptomatic", "Early Stage AD", "AD", "TREM2", "ADAD")
  )

ggboxplot(PC1_df_withpresym, x = "Status", y = "PC1", fill = "Status") +
  ggsignif::stat_signif(
    comparisons = list(
      c("Presymptomatic", "CO"),
      c("Presymptomatic", "Early Stage AD"),
      c("AD", "Early Stage AD"),
      c("Early Stage AD", "CO"),
      c("AD", "CO"),
      c("CO", "TREM2"),
      c("CO", "ADAD"),
      c("AD", "ADAD")
    ),
    y_position = c(6.25, 7, 6.25, 8, 9, 10, 11, 12) + 1,
    tip_length = 0
  ) +
  scale_x_discrete(labels = rev(c("ADAD", "TREM2", "sAD", "ESAD", "Presymptomatic", "CO"))) +
  scale_fill_viridis_d(name = NULL, alpha = 0.5) +
  theme(legend.position = "none")


### Checking if PC1 is different between groups
anno$TubeBarcode <- row.names(anno)

# CO vs presym
PC1_df_CO_presym <- filter(PC1_df_withpresym, Status %in% c("CO", "Presymptomatic")) %>% inner_join(anno)

eigen_model_CO_presym <- glm(Status ~ PC1 + Age + Sex + PMI, data = PC1_df_CO_presym, family = binomial)
summary(eigen_model_CO_presym)

# ESAD vs presym
PC1_df_ESAD_presym <- filter(PC1_df_withpresym, Status %in% c("Early Stage AD", "Presymptomatic")) %>% inner_join(anno)

eigen_model_ESAD_presym <- glm(Status ~ PC1 + Age + Sex + PMI, data = PC1_df_ESAD_presym, family = binomial)
summary(eigen_model_ESAD_presym)

# ESAD vs sAD
PC1_df_ESAD_sAD <- filter(PC1_df_withpresym, Status %in% c("Early Stage AD", "AD")) %>% inner_join(anno)

eigen_model_ESAD_sAD <- glm(Status ~ PC1 + Age + Sex + PMI, data = PC1_df_ESAD_sAD, family = binomial)
summary(eigen_model_ESAD_sAD)

# CO vs ESAD
PC1_df_CO_ESAD <- filter(PC1_df_withpresym, Status %in% c("CO", "Early Stage AD")) %>% inner_join(anno)

eigen_model_CO_ESAD <- glm(Status ~ PC1 + Age + Sex + PMI, data = PC1_df_CO_ESAD, family = binomial)
summary(eigen_model_CO_ESAD)

# CO vs sAD
PC1_df_CO_sAD <- filter(PC1_df_withpresym, Status %in% c("CO", "AD")) %>% inner_join(anno)

eigen_model_CO_sAD <- glm(Status ~ PC1 + Sex + PMI, data = PC1_df_CO_sAD, family = binomial)
summary(eigen_model_CO_sAD)

# CO vs TREM2
PC1_df_CO_TREM2 <- filter(PC1_df_withpresym, Status %in% c("CO", "TREM2")) %>% inner_join(anno)

eigen_model_CO_TREM2 <- glm(Status ~ PC1 + Sex + PMI, data = PC1_df_CO_TREM2, family = binomial)
summary(eigen_model_CO_TREM2)

# CO vs ADAD
PC1_df_CO_ADAD <- filter(PC1_df_withpresym, Status %in% c("CO", "ADAD")) %>% inner_join(anno)

eigen_model_CO_ADAD <- glm(Status ~ PC1 + Sex + PMI, data = PC1_df_CO_ADAD, family = binomial)
summary(eigen_model_CO_ADAD)

# ADAD vs sAD
PC1_df_ADAD_sAD <- filter(PC1_df_withpresym, Status %in% c("ADAD", "AD")) %>% inner_join(anno)

eigen_model_ADAD_sAD <- glm(Status ~ PC1 + Sex + PMI, data = PC1_df_ADAD_sAD, family = binomial)
summary(eigen_model_ADAD_sAD)



### CDR in presym vs ESAD
PC1_df_ESAD_presym$Status <- factor(PC1_df_ESAD_presym$Status, levels = c("Early Stage AD", "Presymptomatic"))
cdr_model_ESAD_presym <- glm(Status ~ CDR + Age + Sex + PMI, data = PC1_df_ESAD_presym, family = binomial)
summary(cdr_model_ESAD_presym)

tau_model_ESAD_presym <- glm(Status ~ BraakTau + Age + Sex + PMI, data = PC1_df_ESAD_presym, family = binomial)
summary(tau_model_ESAD_presym)


pvals <- 
  rbind(
    c("Presymptomatic", "CO"),
    c("Presymptomatic", "Early Stage AD"),
    c("AD", "Early Stage AD"),
    c("Early Stage AD", "CO"),
    c("AD", "CO"),
    c("CO", "TREM2"),
    c("CO", "ADAD"),
    c("AD", "ADAD")
  ) %>% as.data.frame()

colnames(pvals) <- c("group1", "group2")

pvals$p.adj <- c(
  "0.176",
  "0.057",
  "2.38x10\U207B\U2070\U2078",
  "0.661",
  "1.26x10\U207B\U2070\U2077",
  "0.006",
  "0.012",
  "3.21x10\U207B\U2070\U2075"
)

# With custom p-values
ggboxplot(PC1_df_withpresym, x = "Status", y = "PC1", fill = "Status") +
  stat_pvalue_manual(
    pvals,
    y.position = c(6.25, 7, 6.25, 8, 9, 10, 11, 12) + 1,
    label = "{p.adj}",
    tip.length = 0
  ) +
  scale_x_discrete(labels = rev(c("ADAD", "TREM2", "sAD", "ESAD", "Presymptomatic", "CO"))) +
  scale_fill_viridis_d(name = NULL, alpha = 0.5) +
  theme(legend.position = "none")
