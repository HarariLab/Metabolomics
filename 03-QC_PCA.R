library(dplyr)
library(FactoMineR)
library(factoextra)


metabs_clean <- read.csv("data/01-Origscale_clean_avgreps.csv", comment.char="#", stringsAsFactors=FALSE)
ScaledImpData <- read.csv("data/ScaledImp_metabolites_known.csv", stringsAsFactors=FALSE)
pheno <- read.csv("data/00-pheno_extract.csv", stringsAsFactors = FALSE)
ordered_tubes <- as.data.frame(colnames(ScaledImpData[14:487]))
colnames(ordered_tubes) <- "TubeBarcode"
pheno <- inner_join(ordered_tubes, pheno)

replicates <- read.csv("data/replicates.csv", quote="", stringsAsFactors=FALSE)

# Arrange metabs_clean
# make rownames ad Biochemical ID
rownames(metabs_clean) <- metabs_clean$BIOCHEMICAL

# Drop subjects with missing PC values
pheno_clean <- pheno[complete.cases(pheno["PC1"]),]

# Only extract metabolites of interest from scaled, imputed data
ScaledImp_subset <- as.data.frame(ScaledImpData[ScaledImpData$BIOCHEMICAL %in% metabs_clean$BIOCHEMICAL,])
rownames(ScaledImp_subset) <- ScaledImp_subset$BIOCHEMICAL

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
# 931 was duplicated, 60810 is removed)

# Transpose dataframe

ScaledImp_avgreps <- ScaledImp_avgreps[,colnames(ScaledImp_avgreps) %in% pheno_avg$TubeBarcode]
ScaledImp_pca <- as.data.frame(t(ScaledImp_avgreps)) # data is now: subject x metabolite

pca_res_10 <- PCA(ScaledImp_pca, ncp=10, graph = FALSE) # Keep 10 pcs


# Color by Status (no ellipses)
fviz_pca_ind(pca_res_10, geom.ind = "point", col.ind = pheno_avg$Final_CACO, 
             addEllipses = FALSE, legend.title = "CACO Status", mean.point = FALSE)


# Remove outliers
pheno_avg <- pheno_avg %>% filter(!TubeBarcode %in% c("FB06002749", "FB06002881", "FB06003057", "FB06003026", "FB06002744"))

ScaledImp_avgreps <- ScaledImp_avgreps[,colnames(ScaledImp_avgreps) %in% pheno_avg$TubeBarcode]


# write.csv(ScaledImp_avgreps, "data/03-ScaledImp_avgreps.csv")
# write.csv(pheno_avg, "data/03-pheno_avg.csv")







