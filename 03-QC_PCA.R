library(dplyr)
library(FactoMineR)
library(factoextra)



metabs_clean <- read.csv("data/01-Origscale_clean_avgreps.csv", comment.char="#", stringsAsFactors=FALSE)
ScaledImpData <- read.csv("data/ScaledImp_metabolites_known.csv", stringsAsFactors=FALSE)
pheno <- read.csv("data/00-pheno_extract.csv", stringsAsFactors = FALSE)
pheno$Final_CACO[pheno$TREM2_all_variants != 0 & !is.na(pheno$TREM2_all_variants)] <- "TREM2"
pheno$Final_CACO[pheno$Final_CACO == "CA"] <- "sAD"
pheno <- pheno%>% filter(Final_CACO %in% c("ADAD", "sAD", "CO", "Presymptomatic", "TREM2"))
ordered_tubes <- as.data.frame(colnames(ScaledImpData[14:487]))
colnames(ordered_tubes) <- "TubeBarcode"
pheno <- inner_join(ordered_tubes, pheno)

replicates <- read.csv("data/replicates.csv", quote="", stringsAsFactors=FALSE)

# Arrange metabs_clean
# make rownames ad Biochemical ID
rownames(metabs_clean) <- metabs_clean$CHEMICAL.ID

# Drop subjects with missing PC values
drop_ids <- pheno$TubeBarcode[pheno$Final_Status == "Neuro_CO" & pheno$CDR_Expire == 1 | is.na(pheno$Final_Status)]
pheno_clean <- pheno %>% filter(!TubeBarcode %in% drop_ids)

# Only extract metabolites of interest from scaled, imputed data
ScaledImp_subset <- as.data.frame(ScaledImpData[ScaledImpData$CHEMICAL.ID %in% metabs_clean$CHEMICAL.ID,])
rownames(ScaledImp_subset) <- ScaledImp_subset$CHEMICAL.ID

# Only keep subject columns
metabolome_red <- subset(ScaledImp_subset, select = -c(PATHWAY.SORTORDER, BIOCHEMICAL, SUPER.PATHWAY, SUB.PATHWAY, COMP.ID, PLATFORM, CHEMICAL.ID, RI, MASS, PUBCHEM, CAS, KEGG, Group.HMDB))

# Average replicates
avg_reps <- (metabolome_red[,replicates$Replicate1] + metabolome_red[,replicates$Replicate2])/2
colnames(avg_reps) <- replicates$Replicate1

# Drop replicates, append added replicates back
ScaledImp_noreps_temp <- metabolome_red[,!colnames(metabolome_red) %in% replicates[,2]]
ScaledImp_noreps <- ScaledImp_noreps_temp[,!colnames(ScaledImp_noreps_temp) %in% replicates[,3]]
ScaledImp_avgreps <- cbind(ScaledImp_noreps, avg_reps)

# Subset pheno_clean to only get clean (averaged) tubes
rownames(pheno_clean) <- pheno_clean$TubeBarcode
pheno_avg <- pheno_clean[colnames(ScaledImp_avgreps),] # subject x phenotype
pheno_avg <- pheno_avg[!is.na(pheno_avg$MAPID),]

# Transpose dataframe
ScaledImp_avgreps <- ScaledImp_avgreps[,colnames(ScaledImp_avgreps) %in% pheno_avg$TubeBarcode]
ScaledImp_pca <- as.data.frame(t(ScaledImp_avgreps)) # data is now: subject x metabolite

pca_res_10 <- PCA(ScaledImp_pca, ncp=10, graph = FALSE) # Keep 10 pcs


# Color by Status (no ellipses)
fviz_pca_ind(pca_res_10, geom.ind = "point", col.ind = pheno_avg$Final_CACO, 
             addEllipses = FALSE, legend.title = "Status", mean.point = FALSE)


# Remove outliers
pheno_avg <- pheno_avg %>% filter(!TubeBarcode %in% c("FB06002749", "FB06002881", "FB06003057", "FB06003026", "FB06002744"))

ScaledImp_avgreps <- ScaledImp_avgreps[,colnames(ScaledImp_avgreps) %in% pheno_avg$TubeBarcode]

# write.csv(ScaledImp_avgreps, "data/03-ScaledImp_avgreps.csv")
# write.csv(pheno_avg, "data/03-pheno_avg.csv")
# saveRDS(pheno_avg, "data/03-pheno_avg.rds")


ScaledImp_pca_new <- as.data.frame(t(ScaledImp_avgreps)) # data is now: subject x metabolite

pca_res_10_new <- PCA(ScaledImp_pca_new, ncp=10, graph = FALSE) # Keep 10 pcs


# Color by Status (no ellipses)
fviz_pca_ind(pca_res_10_new, geom.ind = "point", col.ind = pheno_avg$Final_CACO, 
             addEllipses = FALSE, legend.title = "Status", mean.point = FALSE)

# Testing PMI association with PC1
ind <- get_pca_ind(pca_res_10)
pc10 <- as.data.frame(ind$coord)
pc10_PMI <- cbind(pc10, PMI = pheno_avg$PMI, stringsAsFactors = TRUE)

model_pmi<- glm(PMI ~ Dim.1, data = pc10_PMI, family = gaussian)
summary(model_pmi)

