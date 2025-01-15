# Load necessary libraries
library(gwasrapidd)
library(readr)
library(dplyr)

# Set the directory where you want to save the file
output_directory <- "C:\\Users\\gqu\\Desktop\\projects\\Genevic\\data"  # Update this path to your specific directory

# Retrieve GWAS catalog data for Alzheimer's disease
AD_VariantAssociations <- get_associations(efo_id = "MONDO_0004975")

# Extract associations and risk alleles from the retrieved data
AD_associations <- data.frame(
  association_id = AD_VariantAssociations@associations$association_id,
  pvalue = AD_VariantAssociations@associations$pvalue,
  or_per_copy_number = AD_VariantAssociations@associations$or_per_copy_number
)

AD_riskAlleles <- data.frame(
  association_id = AD_VariantAssociations@risk_alleles$association_id,
  variant_id = AD_VariantAssociations@risk_alleles$variant_id
)

# Merge association data with risk allele data to link variant IDs
AD_GWAS_Summary <- left_join(AD_associations, AD_riskAlleles, by = "association_id") %>%
  group_by(variant_id) %>%
  summarise(
    GWAS_significant_count = sum(pvalue <= 5e-08, na.rm = TRUE),
    mean_log_or = mean(abs(log(or_per_copy_number)), na.rm = TRUE),
    priority_score = mean_log_or * GWAS_significant_count
  ) %>%
  filter(GWAS_significant_count > 0)  # Keep only variants with significant associations

# Save the results
write.csv(AD_GWAS_Summary, file = paste0(output_directory, "/AD_GWAS_Priority_Scores.csv"), row.names = FALSE)
