#' Process data
#'
#' This function annotates a merged dataset with genomic information.
#'
#' @param full_dataset_path The path to the full dataset CSV file.
#' @param output_file_prefix The prefix for the output file name.
#' @return The annotated dataset.
#'
#' @import tidyr
#' @importFrom dplyr select mutate group_by arrange distinct left_join
#' @importFrom GenomicRanges read_regions annotate_regions
#' @importFrom annotatr build_annotations
#' @export

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# BiocManager::install("annotatr")
BiocManager::install("org.Hs.eg.db", force=TRUE)
# ****        process_data function to annotate merged dataset with genomic information
process_data <- function(full_dataset_path, output_file_prefix) {
  options(max.print = 100)
  library(DBI)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(annotatr)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  # Read full dataset
  full_dataset <- read.csv(full_dataset_path)
  
  # Reshape AD_data to wide format (missing weights will be NA)
  full_dataset_wide <- tidyr::spread(full_dataset, key = ID, value = effect_weight, fill = NA)
  
  # Extract chromosome coordinates
  chrom_coord <- full_dataset_wide %>%
    dplyr::select(hm_chr, hm_pos, hm_rsID)
  
  # Generate bed file
  chrom_coord <- chrom_coord %>%
    dplyr::rename(chrom = hm_chr, chromStart = hm_pos) %>%
    dplyr::mutate(chromEnd = chromStart + 1) %>%
    dplyr::select(chrom, chromStart, chromEnd, hm_rsID)
  
  chrom_coord$chrom <- paste("chr", chrom_coord$chrom, sep = "")
  
  # Remove rows where "chrom_coord" column is "NA"
  chrom_coord <- subset(chrom_coord, chromEnd != "NA")
  
  # Save DataFrame as a tab-delimited file without header names
  require(dplyr)
  chrom_coord <- mutate_if(chrom_coord, is.numeric, as.integer)
  write.table(chrom_coord, file = "./chrom_coord.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")
  
  # Read bed file created above for annotations 
  dm_regions <- read_regions(
    con = "./chrom_coord.bed",
    genome = "hg38",
    format = "bed",
    extraCols = c(rsID = "character")
  )
  
  annots <- c("hg38_genes_1to5kb", "hg38_genes_promoters", "hg38_genes_cds", "hg38_genes_5UTRs", "hg38_genes_exons", 
              "hg38_genes_firstexons", "hg38_genes_introns", "hg38_genes_intronexonboundaries", "hg38_genes_exonintronboundaries", 
              "hg38_genes_3UTRs", "hg38_genes_intergenic", "hg38_enhancers_fantom", "hg38_basicgenes", "hg38_cpgs" )
  annotations <- build_annotations(genome = "hg38", annotations = annots)
  
  # Intersect regions with annotations
  dm_annotated <- annotate_regions(
    regions = dm_regions,
    annotations = annotations,
    minoverlap = 1L,
    ignore.strand = TRUE,
    quiet = TRUE
  )
  
  # Convert GRanges object to data frame
  df_dm_annotated <- data.frame(dm_annotated)
  print(names(df_dm_annotated))
  print(full_dataset_path)
  # # See the GRanges dataframe
  # print(head(df_dm_annotated))
  # Select required columns
  df_dm_annotated_select <- df_dm_annotated %>%
    dplyr::select(seqnames, start, annot.symbol, annot.type, annot.width)
  # # See the df_dm_annotated_select dataframe
  # print(head(df_dm_annotated_select))
  
  # Create new column "SNP_coord" with unique SNP name
  full_dataset <- full_dataset %>%
    dplyr::mutate(SNP_coord = paste(hm_chr, hm_pos, sep = "_"))
  # Create new column "SNP_coord" with unique SNP name "start-1" because in annotation dataframe, it is 1 nucleotide ahead
  
  
  df_dm_annotated_select <- df_dm_annotated_select %>%
    mutate(SNP_coord = paste(sub("^chr", "", seqnames), (start-1), sep = "_"))
  # Remove duplicate SNP_coord that had multiple annotations
  # df_dm_annotated_filter <- df_dm_annotated_select %>%
  # group_by(SNP_coord) %>%
  # filter(row_number() == 1)
  
  
  df_dm_annotated_filter <- df_dm_annotated_select %>%
    dplyr::group_by(SNP_coord) %>%
    dplyr::filter(!is.na(annot.symbol)) %>%
    dplyr::arrange(annot.width, .by_group = TRUE) %>%
    dplyr::slice(1) 
  
  # # See the df_dm_annotated_select dataframe
  print(head(df_dm_annotated_select))
  # Merge datasets
  annotated_df <- dplyr::left_join(full_dataset, df_dm_annotated_filter, by = "SNP_coord")
  annotated_df <-  annotated_df %>% dplyr::select("SNP_coord", "hm_rsID","hm_chr", "hm_pos", "effect_allele", "effect_weight",  
                                                  "ID", "annot.symbol", "annot.type")
  # Filter out rows with "NA_NA" in SNP_coord column
  annotated_df <- annotated_df %>%
    filter(!is.na(SNP_coord) & SNP_coord != "NA_NA")
  # Filter for rows with the highest "weight" within each "ID" group, removing duplicate SNPs
  annotated_df <- annotated_df %>%
    dplyr::group_by(ID) %>%
    dplyr::arrange(desc(effect_weight)) %>%
    dplyr::distinct(SNP_coord, .keep_all = TRUE) %>%
    ungroup()
  
  # Add ranks per study: highest weight in each study is ranked 1, 2nd highest is ranked 2, and so on...
  annotated_df <- annotated_df %>%
    dplyr::arrange(ID, effect_weight) %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(ranks = rank(-abs(as.double(effect_weight)), ties.method = "average")) %>%
    ungroup()
  
  # Save the full dataset as a CSV file
  output_file <- paste0(output_file_prefix, "_annotated_dataset.csv")
  write.csv(annotated_df, file = output_file, row.names = FALSE)
  
  return(annotated_df)
}

# Set the input and output directories
input_directory <- 'F:/projects/GENEVIC/data/AlzheimerPGS/Out/'
output_directory <- 'F:/projects/GENEVIC/data/AlzheimerPGS/Annotated/'

# List all CSV files in the input directory
file_list <- list.files(path = input_directory, pattern = "\\.csv$", full.names = TRUE)

# Loop through each file and process it
for (file_path in file_list) {
  file_name <- basename(file_path)
  output_file_prefix <- file.path(output_directory, tools::file_path_sans_ext(file_name))
  
  # Call the existing process_data function on each file
  print(paste("Processing:", file_name))
  process_data(file_path, output_file_prefix)
}

