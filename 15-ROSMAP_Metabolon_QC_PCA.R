library(FactoMineR)
library(factoextra)

metab_data_transformed <- readRDS("data/14-metab_data_transformed.rds")
pheno <- readRDS("data/14-ROSMAP_pheno_clean.rds")

metab_data_impute <- metab_data_transformed

# Set NAs to minimum
for(i in 1:ncol(metab_data_impute)) {
  metab_data_impute[which(is.na(metab_data_impute[,i])),i] <- min(metab_data_impute[,i], na.rm = TRUE)
}

# No obvious outliers
pca_obj <- PCA(metab_data_impute, graph = FALSE)

pca_ind_plot <- fviz_pca_ind(pca_obj, geom = "point", col.ind = pheno$Status, addEllipses = TRUE)
# ggsave(pca_ind_plot, filename = "plots/02-pca_ind_plot.png", device = "png")

saveRDS(metab_data_impute, "data/15-metab_data_impute.rds")
