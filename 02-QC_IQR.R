library(dplyr)

## Log-transform data, remove outliers, and adjust mean to 0

Origscale_clean_avgreps <- read.csv("data/01-Origscale_clean_avgreps.csv")

Origscale_log = cbind(Origscale_clean_avgreps[,2:3], 
                      log10(Origscale_clean_avgreps[,-1:-13]))

log10_t <- as.matrix(t(Origscale_log[,-1:-2])) # Rows = tubes, cols = metabs
transformed <- log10_t
transformed_mean <- vector()

metab_IQRs <- apply(log10_t, 2, FUN = function(x) IQR(x, na.rm = TRUE))
metab_quantiles <- apply(log10_t, 2, FUN = function(x) quantile(x, na.rm = TRUE))

j <- 1
while(j <= nrow(Origscale_log)){
  transformed[log10_t[, j] > metab_quantiles[4, j] + metab_IQRs[j] * 1.5, j] <- NA # Set everything 1.5xIQR above 75% to NA
  transformed[log10_t[, j] < metab_quantiles[2, j] - metab_IQRs[j] * 1.5, j] <- NA # Set everything 1.5xIQR below 25% to NA
  # Calculate mean and normalize
  transformed_mean[j] <- mean(transformed[, j], na.rm=TRUE) # Get mean for metabolite
  transformed[, j] <- transformed[, j] - transformed_mean[j] # Subtract mean from values to make mean 0
  j <- j + 1
}

# Export IQR adjusted, rescaled, and log-transformed data
Origscale_clean_transformed <- as.data.frame(transformed)
colnames(Origscale_clean_transformed) <- Origscale_clean_avgreps$CHEMICAL.ID
# write.csv(Origscale_clean_transformed, "data/02-Origscale_clean_transformed.csv")
                      
                      