#!/bin/bash
#==============================================================================
# TAXONOMIC CLASSIFICATION WITH KRAKEN2/BRACKEN
# Sjögren Syndrome Salivary Metagenome Analysis
# Input: Host-removed reads from preprocessing/03_host_removal/
# Output: Species abundance tables for diversity analysis
#==============================================================================

set -euo pipefail

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Kraken2/Bracken parameters
READ_LENGTH=151      # Read length after trimming (adjust for your data)
LEVEL="S"            # Taxonomic level: S = Species
MIN_READS=10         # Minimum reads for Bracken abundance estimation

echo -e "${BLUE}=== TAXONOMIC CLASSIFICATION WITH KRAKEN2/BRACKEN ===${NC}"
echo -e "${BLUE}Date: $(date)${NC}"
echo ""

#------------------------------------------------------------------------------
# Check architecture (for Apple Silicon compatibility)
#------------------------------------------------------------------------------
ARCH=$(uname -m)
if [[ "$ARCH" == "arm64" ]]; then
    echo -e "${YELLOW}Detected ARM64 architecture.${NC}"
    echo -e "${YELLOW}Note: Kraken2 requires x86_64. Run under Rosetta if needed.${NC}"
    # Uncomment below to auto-switch to Rosetta:
    # if [[ "${RUNNING_UNDER_ROSETTA:-}" != "true" ]]; then
    #     export RUNNING_UNDER_ROSETTA=true
    #     exec arch -x86_64 /bin/bash "$0" "$@"
    # fi
fi

#------------------------------------------------------------------------------
# Check requirements
#------------------------------------------------------------------------------
if ! command -v kraken2 &> /dev/null; then
    echo -e "${RED}Error: kraken2 not found.${NC}"
    exit 1
fi

if ! command -v bracken &> /dev/null; then
    echo -e "${RED}Error: bracken not found.${NC}"
    exit 1
fi

# Check Kraken2 database
if [[ ! -d "$KRAKEN_DB" ]]; then
    echo -e "${RED}Error: Kraken2 database not found at $KRAKEN_DB${NC}"
    echo ""
    echo -e "${YELLOW}To build an oral microbiome database:${NC}"
    echo "kraken2-build --download-taxonomy --db oral_microbiome_db"
    echo "kraken2-build --download-library bacteria --db oral_microbiome_db"
    echo "kraken2-build --build --db oral_microbiome_db --threads 8"
    exit 1
fi

# Check Bracken database for read length
BRACKEN_DB_FILE="${KRAKEN_DB}/database${READ_LENGTH}mers.kmer_distrib"
if [[ ! -f "$BRACKEN_DB_FILE" ]]; then
    echo -e "${RED}Error: Bracken database not found for ${READ_LENGTH}bp reads${NC}"
    echo ""
    echo -e "${YELLOW}Build Bracken database:${NC}"
    echo "bracken-build -d \"$KRAKEN_DB\" -t $THREADS -k 35 -l $READ_LENGTH"
    exit 1
fi

echo -e "${GREEN}✓ Kraken2 and Bracken databases found${NC}"

#------------------------------------------------------------------------------
# Create output directories
#------------------------------------------------------------------------------
KRAKEN_OUTPUT="$TAXONOMIC_DIR/kraken2_output"
BRACKEN_OUTPUT="$TAXONOMIC_DIR/bracken_output"
ABUNDANCE_TABLES="$TAXONOMIC_DIR/abundance_tables"

mkdir -p "$KRAKEN_OUTPUT"
mkdir -p "$BRACKEN_OUTPUT"
mkdir -p "$ABUNDANCE_TABLES"
mkdir -p "$LOGS_DIR/taxonomic"

#------------------------------------------------------------------------------
# Process samples
#------------------------------------------------------------------------------
TOTAL_SAMPLES=$(wc -l < "$SAMPLE_LIST")
echo -e "${BLUE}Processing $TOTAL_SAMPLES samples${NC}"
echo -e "${BLUE}Database: $(basename $KRAKEN_DB)${NC}"
echo -e "${BLUE}Read length: ${READ_LENGTH}bp${NC}"
echo -e "${BLUE}Taxonomic level: Species${NC}"
echo ""

# Initialize results summary
echo -e "Sample\tTotal_Reads\tClassified_Reads\tPercent_Classified\tTop_Species" > "$TAXONOMIC_DIR/classification_summary.tsv"

CURRENT=0
while IFS= read -r sample; do
    [[ -z "$sample" || "$sample" =~ ^[[:space:]]*# ]] && continue

    CURRENT=$((CURRENT + 1))
    echo -e "${BLUE}[$CURRENT/$TOTAL_SAMPLES] Processing: $sample${NC}"

    # Input files
    R1="${HOST_REMOVAL_DIR}/${sample}_nonhost_1.fastq.gz"
    R2="${HOST_REMOVAL_DIR}/${sample}_nonhost_2.fastq.gz"

    # Check input files
    if [[ ! -f "$R1" ]] || [[ ! -f "$R2" ]]; then
        echo -e "${YELLOW}  Warning: Input files not found, skipping${NC}"
        continue
    fi

    # Output files
    KRAKEN_REPORT="${KRAKEN_OUTPUT}/${sample}.kreport"
    KRAKEN_OUT="${KRAKEN_OUTPUT}/${sample}.kraken"
    BRACKEN_OUT="${BRACKEN_OUTPUT}/${sample}.bracken"
    BRACKEN_REPORT="${BRACKEN_OUTPUT}/${sample}.breport"
    LOG_FILE="$LOGS_DIR/taxonomic/${sample}.log"

    # Skip if already processed
    if [[ -f "$BRACKEN_OUT" ]]; then
        echo -e "${YELLOW}  Already processed, skipping...${NC}"
        continue
    fi

    # Run Kraken2
    echo "  Running Kraken2..."
    kraken2 \
        --db "$KRAKEN_DB" \
        --threads $THREADS \
        --paired \
        --gzip-compressed \
        --report "$KRAKEN_REPORT" \
        --output "$KRAKEN_OUT" \
        "$R1" "$R2" \
        2> "$LOG_FILE"

    # Extract classification statistics
    TOTAL_READS=$(grep -o '[0-9]\+ sequences processed' "$LOG_FILE" | grep -o '[0-9]\+' || echo "0")
    CLASSIFIED_READS=$(grep -o '[0-9]\+ sequences classified' "$LOG_FILE" | grep -o '[0-9]\+' || echo "0")

    if [[ "$TOTAL_READS" -gt 0 ]]; then
        PERCENT_CLASSIFIED=$(echo "scale=2; $CLASSIFIED_READS * 100 / $TOTAL_READS" | bc -l)
    else
        PERCENT_CLASSIFIED="0.00"
    fi

    # Run Bracken for species-level abundance
    echo "  Running Bracken..."
    bracken \
        -d "$KRAKEN_DB" \
        -i "$KRAKEN_REPORT" \
        -o "$BRACKEN_OUT" \
        -w "$BRACKEN_REPORT" \
        -r $READ_LENGTH \
        -l $LEVEL \
        -t $MIN_READS \
        2>> "$LOG_FILE" || {
            echo -e "${YELLOW}  Warning: Bracken failed for $sample${NC}"
            continue
        }

    # Get top species
    if [[ -f "$BRACKEN_OUT" ]] && [[ -s "$BRACKEN_OUT" ]]; then
        TOP_SPECIES=$(awk 'NR>1 && NF>0 {print $1; exit}' "$BRACKEN_OUT" 2>/dev/null || echo "N/A")
        [[ -z "$TOP_SPECIES" ]] && TOP_SPECIES="N/A"
        echo -e "${sample}\t${TOTAL_READS}\t${CLASSIFIED_READS}\t${PERCENT_CLASSIFIED}\t${TOP_SPECIES}" >> "$TAXONOMIC_DIR/classification_summary.tsv"
        echo -e "${GREEN}  ✓ Completed${NC}"
    else
        echo -e "${sample}\t${TOTAL_READS}\t${CLASSIFIED_READS}\t${PERCENT_CLASSIFIED}\tN/A" >> "$TAXONOMIC_DIR/classification_summary.tsv"
        echo -e "${YELLOW}  ⚠ Completed but no Bracken output${NC}"
    fi

done < "$SAMPLE_LIST"

#------------------------------------------------------------------------------
# Create combined abundance table
#------------------------------------------------------------------------------
echo ""
echo -e "${YELLOW}Creating combined abundance table...${NC}"

python3 << 'PYTHON_SCRIPT'
import os
import pandas as pd

bracken_dir = os.environ.get('BRACKEN_OUTPUT', '')
sample_list = os.environ.get('SAMPLE_LIST', '')
output_dir = os.environ.get('ABUNDANCE_TABLES', '')

# Read sample list
with open(sample_list, 'r') as f:
    samples = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]

# Initialize dictionary for abundances
species_dict = {}

# Read each Bracken output
for sample in samples:
    bracken_file = f"{bracken_dir}/{sample}.bracken"
    if os.path.exists(bracken_file):
        df = pd.read_csv(bracken_file, sep='\t')
        for _, row in df.iterrows():
            species = row['name']
            if species not in species_dict:
                species_dict[species] = {}
            species_dict[species][sample] = row['fraction_total_reads']
    else:
        print(f"Warning: No Bracken output for {sample}")

# Create DataFrame
abundance_df = pd.DataFrame.from_dict(species_dict, orient='index')
abundance_df = abundance_df.fillna(0)

# Ensure column order matches sample list
available_samples = [s for s in samples if s in abundance_df.columns]
abundance_df = abundance_df[available_samples]

# Save to file
output_file = f"{output_dir}/species_abundance.tsv"
abundance_df.to_csv(output_file, sep='\t')

print(f"Species abundance table saved to: {output_file}")
print(f"Total species detected: {len(abundance_df)}")
print(f"Samples included: {len(available_samples)}")
print(f"\nTop 10 most abundant species:")
print(abundance_df.sum(axis=1).nlargest(10))
PYTHON_SCRIPT

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------
echo ""
echo -e "${GREEN}=== TAXONOMIC CLASSIFICATION COMPLETED ===${NC}"
echo ""
echo -e "${YELLOW}Results:${NC}"
echo -e "  • Kraken2 reports: $KRAKEN_OUTPUT/"
echo -e "  • Bracken reports: $BRACKEN_OUTPUT/"
echo -e "  • Species abundance table: $ABUNDANCE_TABLES/species_abundance.tsv"
echo -e "  • Classification summary: $TAXONOMIC_DIR/classification_summary.tsv"
echo ""
echo -e "${BLUE}Next steps:${NC}"
echo -e "  1. Review species abundance table"
echo -e "  2. Run statistical analysis with R (02_statistical_analysis.Rmd)"
