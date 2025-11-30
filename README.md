# Salivary Metagenome Analysis in Sjögren Syndrome

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Kraken2](https://img.shields.io/badge/Kraken2-2.1.3-green)](https://ccb.jhu.edu/software/kraken2/)
[![R](https://img.shields.io/badge/R-%3E%3D4.3-blue)](https://www.r-project.org/)

Analysis pipeline for shotgun metagenomic sequencing of salivary samples from patients with Sjögren syndrome and healthy controls.

> **Note**: This is a work-in-progress project. The pipeline covers preprocessing through microbial diversity analysis.

## Overview

This repository contains a complete bioinformatics workflow for whole-genome shotgun (WGS) metagenomic analysis:

- **Preprocessing pipeline** with FastQC, fastp, and Bowtie2 (host removal)
- **Taxonomic classification** with Kraken2/Bracken using oral microbiome database
- **Statistical analysis** in R for diversity metrics and differential abundance

This project uses whole-genome shotgun sequencing rather than 16S rRNA amplicon sequencing, providing higher resolution (species/strain-level identification), no primer bias, and more accurate abundance estimates.

## Study Design

| Parameter | Description |
|-----------|-------------|
| Disease | Sjögren syndrome (autoimmune disorder affecting salivary glands) |
| Cohort | 144 saliva samples |
| Groups | Sjögren patients vs. healthy controls |
| Platform | Illumina (paired-end, 151 bp) |
| Sequencing type | Whole-genome shotgun (WGS) metagenomics |
| Average depth | ~40-58 million reads per sample |

Kraken2/Bracken was selected for taxonomic classification due to its speed (~1M reads/min) and species-level resolution. A custom oral microbiome database was used to reduce false positives from misclassification to gut/environmental organisms.

Host removal with Bowtie2 is critical for saliva samples, which typically contain 50-90% human DNA contamination.

## Repository Structure

```
shotgun-metagenomics-saliva-sjogren/
├── README.md
├── LICENSE
├── scripts/
│   ├── preprocessing/
│   │   ├── 01_setup_environment.sh       # Project setup and conda env
│   │   ├── 02_quality_control.sh         # FastQC + MultiQC
│   │   ├── 03_fastp_trimming.sh          # Quality trimming
│   │   ├── 04_host_removal.sh            # Human DNA removal
│   │   └── 05_final_qc.sh                # Post-processing QC
│   ├── 01_taxonomic_classification.sh    # Kraken2 + Bracken
│   └── 02_statistical_analysis.Rmd       # R diversity analysis
└── data/
    └── README_data.md                    # Data access instructions
```

## Pipeline Overview

```
Raw FASTQ files (Illumina paired-end)
        │
        ▼
┌───────────────────────────────────────┐
│  PREPROCESSING                        │
│  ├── FastQC (quality assessment)      │
│  ├── fastp (trimming, Q20, 50bp min)  │
│  └── Bowtie2 (human DNA removal)      │
└───────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────┐
│  TAXONOMIC CLASSIFICATION             │
│  ├── Kraken2 (k-mer classification)   │
│  └── Bracken (abundance estimation)   │
└───────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────┐
│  STATISTICAL ANALYSIS (R)             │
│  ├── Alpha diversity (Shannon, etc.)  │
│  ├── Beta diversity (PCoA, PERMANOVA) │
│  └── Differential abundance           │
└───────────────────────────────────────┘
```

## Methods

### 1. Preprocessing

| Step | Tool | Parameters |
|------|------|------------|
| Quality assessment | FastQC v0.12.1 | Default |
| Report aggregation | MultiQC v1.15 | Default |
| Quality trimming | fastp v0.23.4 | Q20, min-length 50bp, adapter auto-detection |
| Host removal | Bowtie2 v2.5.1 | --very-sensitive-local, GRCh38 reference |

### 2. Taxonomic Classification

| Step | Tool | Parameters |
|------|------|------------|
| Classification | Kraken2 v2.1.3 | Oral microbiome database, paired-end |
| Abundance estimation | Bracken v2.8 | Read length 151bp, species level, min 10 reads |

### 3. Statistical Analysis (R)

| Analysis | Method | Purpose |
|----------|--------|---------|
| Alpha diversity | Shannon, Simpson, Richness | Within-sample diversity |
| Beta diversity | Bray-Curtis, PCoA | Between-sample differences |
| Statistical tests | Wilcoxon, PERMANOVA | Group comparisons |
| Differential abundance | Wilcoxon + BH correction | Identify disease-associated taxa |

## Requirements

### Key Dependencies

**Preprocessing:**
- FastQC v0.12.1
- MultiQC v1.15
- fastp v0.23.4
- Bowtie2 v2.5.1
- samtools v1.18
- seqkit v2.5.1

**Taxonomic Analysis:**
- Kraken2 v2.1.3
- Bracken v2.8

**Statistical Analysis (R):**
```r
library(vegan)
library(ape)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
```

### Database Requirements

1. **Human genome index** (for host removal):
   ```bash
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
   bowtie2-build --threads 8 GCF_000001405.40_GRCh38.p14_genomic.fna.gz GRCh38
   ```

2. **Kraken2 oral microbiome database**:
   ```bash
   kraken2-build --download-taxonomy --db oral_microbiome_db
   kraken2-build --download-library bacteria --db oral_microbiome_db
   kraken2-build --build --db oral_microbiome_db --threads 8
   bracken-build -d oral_microbiome_db -t 8 -k 35 -l 151
   ```

## Quick Start

```bash
# Run preprocessing
bash scripts/preprocessing/02_quality_control.sh
bash scripts/preprocessing/03_fastp_trimming.sh
bash scripts/preprocessing/04_host_removal.sh

# Run taxonomic classification
bash scripts/01_taxonomic_classification.sh
```

```r
# Run R analysis
rmarkdown::render("scripts/02_statistical_analysis.Rmd")
```

## Output Files

| Output | Description |
|--------|-------------|
| `species_abundance.tsv` | Species abundance matrix (samples × species) |
| `classification_summary.tsv` | Per-sample classification statistics |
| `Alpha_Diversity_*.tsv` | Alpha diversity metrics and statistics |
| `Beta_Diversity_*.tsv` | Beta diversity statistics (PERMANOVA) |
| `Differential_Abundance_Results.tsv` | Differentially abundant species |
| `*.pdf` | Publication-ready figures |

## Data Access

Raw sequencing data will be deposited in NCBI SRA upon publication.

See `data/README_data.md` for detailed data access instructions.

## Citation

*Manuscript in preparation*

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Maria J. Rus**
[![ORCID](https://img.shields.io/badge/ORCID-0000--0003--3659--2821-green)](https://orcid.org/0000-0003-3659-2821)
Email: marjimrus@gmail.com

## Acknowledgments

- Universidad de Sevilla
- Virgen Macarena University Hospital (sample collection)
