# Metabolomic and lipidomic signatures in autosomal dominant and late-onset Alzheimer disease brains

### Scripts used to QC and analyze metabolomics data

---

## Scripts
#### **Knight ADRC dataset**
##### 01-QC_missingness.R (Knight ADRC Data)
- Check technical replicate correlation and average
- Remove metabolites missing >20% readings
- Recover metabolites with high missingness whose missingness between groups corresponds with differential abundance between groups
- Impute missing values in recovered metabolites with minimum
- Check missingness by individual

##### 02-QC_IQR.R
- Log-transform data, remove outliers (>1.5xIQR), adjust mean

##### 03-QC_PCA.R
- Average replicates in the scaled/imputed data provided by Metabolon
- Calculate PCA and identify outlier samples
- Remove outlier samples

##### 04-main_analysis.R
- Main analysis including linear regressions for each status group (also matched by CDR, BraakTau, and BraakAbeta), association with age at death in sAD, and association with APOE genotype

##### 05-effect_comparison.R
- Compare effects of 16 common metabolites between ADAD, TREM2, and sAD. ANCOVA tests with and without matching individuals by CDR, BraakTau, and BraakAbeta

##### 06-eigengene_analysis.R
- Calculate eigengene (PC1) for 16 metabolites in CO, sAD, ADAD, and TREM2 participants
- Test association of eigengene with phenotypic variables

##### 07-heatmap_analysis.R
- Create heatmaps representing the 16 metabolite-profile and test association of early-stage AD (ESAD) group with phenotypic variables

##### 08-pharmacological_analysis.R
- Check for effects of fluoxetine on beta-citrylglutamate associations and vitamin supplements on vitamin associations

#### **ROSMAP and ADNI p180 datasets**

##### 09.1-ROSMAP_p180_analysis.R
- Full cleaning and analysis of ROSMAP serum samples

##### 09.2-ROSMAP_p180_analysis.R
- Full cleaning and analysis of ROSMAP DLPFC samples (p180 only)

##### 10-ADNI1_analysis.R
- Full cleaning and analysis of ADNI1 serum samples

##### 11-ADNI2GO_analysis.R
- Full cleaning and analysis of ADNIGO/2 serum samples

##### 12-serum_meta-analysis.R
- Meta-analysis of ROSMAP, ADNI1, and ADNI2/GO serum data

#### **ROSMAP Metabolon dataset**

##### 13-ROSMAP_Metabolon_clean_phenotypic_data.R
- Create phenotype data frame to be used with ROSMAP Metabolon data

##### 14.1-ROSMAP_Metabolon_QC_missingness_cogdx_recovery.R
- Recovery of metabolites using consensus clinical diagnosis for ROSMAP Metabolon data, same process as 01-QC_missingness.R

##### 14.2-ROSMAP_Metabolon_QC_missingness.R
- Missingness QC for ROSMAP Metabolon data, same process as 01-QC_missingness.R

##### 15-ROSMAP_Metabolon_QC_PCA.R
- Calculate and visualize PCA for individuals; no outliers removed

##### 16-ROSMAP_Metabolon_main_analysis.R
- Differential abundance analysis for ROSMAP Metabolon data (linear regressions)

##### 17-ROSMAP_Metabolon_eigengene_analysis.R
- Calculate eigengene for the 15 available metbolites in ROSMAP Metabolon data, check association with disease duration and Braak tau

---

## Metabolomics Browser

##### app.R
- Main app file, load data and call modules

#### **Modules Folder**

##### boxplots.R
- Module for "Reading Distributions" tab, display results table, boxplots, and specific metabolite information

##### home.R
- Homepage, display phenotypic summary information

##### volcano.R
- Module for "Volcano Plots" tab, plot DE results, show boxplots, and specific metabolite information




