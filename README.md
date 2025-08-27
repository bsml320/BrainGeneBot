# 🧬 BrainGeneBot

*BrainGeneBot: A Framework for Variant Prioritization and GPT-informed Interpretation across Polygenic Risk Score Studies*

---

## Repository Structure
```text
├── R/
│   ├── annotate_data.R                 # Annotate SNP datasets with genomic features (hg38)
│   └── Calculate_GWAS_priority_score.R # Compute GWAS-based variant priority scores for AD
├── configs/                            # Configuration files (hyperparameters, settings, etc.)
├── models/
│   ├── dataset.py      # Feature engineering from SNP-level data into numeric feature sets
│   ├── sup_models.py   # Supervised models (ML + NN) for prediction and SNP ranking
│   ├── Data Preprocess.ipynb   # Jupyter notebook for data preprocessing
│   └── unsup_models.py # Unsupervised rank aggregation methods (Borda, Markov Chain, Bayesian, etc.)

├── Analysis_Chat.py                    # Query ClinVar, extract gene/variant info, generate structured prompts & analysis
├── GeneChat.py                         # Generate prompts for gene relevance, STRING network, or enrichment analysis
├── Literature_Chat.py                  # Literature mining and retrieval
├── PGSChat.py                          # Polygenic score prompt generation and analysis of ranked variants
├── build_dbSNP.py                      # Tool to construct SNP reference from dbSNP
├── experiments_sup.py                  # Supervised training/testing experiments
├── experiments_unsup.py                # Unsupervised training/testing experiments
├── impute_rsID.py                      # rsID imputation and mapping tool
├── main_sup.py                         # Entry point for supervised tasks
├── main_unsup.py                       # Entry point for unsupervised tasks
├── utils.py                            # Utility functions and helpers
└── README.md                           # Project documentation
```

## Installation
```text
# Clone the repository
git clone <repo_url>
cd <repo_name>
```
## Requirements
### Python (main)
- `pandas`
- `numpy`
- `tqdm`
- `scikit-learn`
- `torch`
- `matplotlib`
- `seaborn`
- `joblib`
- `Bio` (Biopython)
- `google-generativeai`
- `requests`
- `enrichrpy` (for enrichment analysis)

### R (main)
- `dplyr`
- `tidyr`
- `annotatr`
- `GenomicRanges`
- `TxDb.Hsapiens.UCSC.hg38.knownGene`
- `org.Hs.eg.db`
- `gwasrapidd`




## Data Preparation

The data preparation pipeline processes Alzheimer’s polygenic score (PGS) files from the PGS Catalog into standardized datasets suitable for supervised and unsupervised analysis.

### 1. Download & Decompress
- Read Alzheimer’s-specific entries from `pgs_scores_Alzheimer.csv`.
- Fetch harmonized scoring files (`*_hmPOS_GRCh38.txt.gz`) from FTP.
- Prevent re-downloading existing files.
- Decompress all `.gz` archives into plain text for downstream processing.

### 2. Annotation & Cleaning
- Annotate SNPs using `annotate_data.R` (hg38 gene models).
- Add functional context (promoter, exon, intron, CpG, enhancer, etc.).
- Remove unknown or unmapped rsIDs.
- Exclude problematic files (e.g., `PGS003957`, `PGS003958`).
- Generate `annotation_map.csv` with unique SNP coordinates and alleles.

### 3. Duplicate Detection & Merging
- Identify duplicate PGS files by comparing SNP sets.
- Convert effect sizes into **Z-scores** (absolute, allele-count aware).
- Merge duplicates by summing Z-scores and aggregating ranks into lists.

### 4. Top-N Filtering & Overlap Optimization
- Apply quantile-based filtering to retain top variants by effect weight.
- Optimize N (variants per study) to maximize **cross-study SNP overlap** (~90% of maximum overlap).

### 5. Standardization & Reranking
- Standardize effect weights within each study.
- Replace list of ranks with averaged rank values.
- Re-rank SNPs consistently across studies.

### 6. Final Outputs
- **`data_dict.pkl`** → processed SNP-level study data.
- **`annotation_map.csv`** → SNP-to-annotation mapping.
- **`AD_GWAS_Priority_Scores.csv`** → GWAS-based priority scores for Alzheimer’s disease variants.




## Usage
### 1. Supervised / Unsupervised Ranking
```text
# Ranking aggregation
# Run supervised experiments
python main_sup.py --config configs/sup_config.yaml

# Run unsupervised experiments
python main_unsup.py --config configs/unsup_config.yaml
```
### 2. Using Chat Modules
The chat modules integrate AI-driven reasoning with genetic datasets, GWAS scores, and enrichment tools.
They generate structured prompts, query external resources (ClinVar, STRING, Enrichr), and return JSON-formatted outputs.
These modules require loading ranking results (from data_dict.pkl + aggregated ranks).
```text
# Variant-level ClinVar analysis
python Analysis_Chat.py

# Gene relevance and STRING network
python GeneChat.py

# Literature retrieval
python Literature_Chat.py

# Polygenic score reasoning (load ranks, summarize top genes/variants)
python PGSChat.py
```

Note: Before running chat modules, ensure that processed results (data_dict.pkl, annotation_map.csv, aggregated ranks) are available in data/.

Example:

Supervised pipeline produces ranked variants → load ranks into PGSChat.py.

Unsupervised pipeline produces aggregated ranks → use as input for gene-level prompts.
## Notes
- data_dict.pkl and annotation_map.csv must be generated before running supervised/unsupervised models.
- GWAS-based priority scores (AD_GWAS_Priority_Scores.csv) are produced separately in R/Calculate_GWAS_priority_score.R.
- Example workflows and visualizations are provided in Data Preprocess.ipynb.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use this code or framework in your research, please cite:
```text
@article{your2024braingenebot,
  title={BrainGeneBot: A Framework for Variant Prioritization and GPT-informed Interpretation across Polygenic Risk Score Studies},
  author={Your Name and Collaborators},
  journal={Journal Name},
  year={2024},
  publisher={Publisher}
}
```
