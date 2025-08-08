#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=480000
#SBATCH --job-name=provean_analysis
#SBATCH --output=logs/provean_analysis.out
#SBATCH --error=logs/provean_analysis.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@institution.edu

# HBV MTCT Analysis Pipeline - PROVEAN Functional Impact Prediction
#
# Purpose: Predict functional impact of HBV protein mutations using PROVEAN
# Method: PROVEAN algorithm for variant effect prediction
# Input: Protein sequences and variant files
# Output: Functional impact scores and predictions
#
# Citation: Choi, Y. & Chan, A.P. (2015). PROVEAN web server: a tool to predict 
#           the functional effect of amino acid substitutions and indels. 
#           Bioinformatics, 31(16), 2745-2747.
#
# Author: Your Name
# Date: 2024
# Version: 1.0

set -euo pipefail

# Script configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Load configuration
CONFIG_FILE="$PROJECT_ROOT/config/pipeline_config.yaml"
if [ -f "$CONFIG_FILE" ]; then
    OUTPUT_DIR=$(grep -A 10 "paths:" "$CONFIG_FILE" | grep "output_dir:" | awk '{print $2}' | tr -d '"')
    PROVEAN_THRESHOLD=$(grep -A 10 "functional:" "$CONFIG_FILE" | grep "provean_threshold:" | awk '{print $2}')
    THREADS=$(grep -A 10 "resources:" "$CONFIG_FILE" | grep "threads:" | awk '{print $2}')
else
    OUTPUT_DIR="results"
    PROVEAN_THRESHOLD=-2.5
    THREADS=8
fi

# Directories
PROTEIN_DIR="$PROJECT_ROOT/results/protein_analysis"
FUNCTIONAL_DIR="$OUTPUT_DIR/functional_analysis"
LOGS_DIR="$OUTPUT_DIR/logs"
DATABASE_DIR="$PROJECT_ROOT/data/databases"

# PROVEAN paths
PROVEAN_BIN="$CONDA_PREFIX/bin/provean.sh"
BLAST_DB="$DATABASE_DIR/swissprot"

# Create output directories
mkdir -p "$FUNCTIONAL_DIR/provean_results" "$LOGS_DIR"

# Logging functions
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOGS_DIR/provean_analysis.log"
}

log_info() { log "INFO: $1"; }
log_warn() { log "WARN: $1"; }
log_error() { log "ERROR: $1"; }
log_success() { log "SUCCESS: $1"; }

# Validate requirements
validate_requirements() {
    log_info "Validating PROVEAN requirements..."
    
    # Check PROVEAN installation
    if [ ! -f "$PROVEAN_BIN" ]; then
        log_error "PROVEAN not found at: $PROVEAN_BIN"
        exit 1
    fi
    
    # Check BLAST database
    if [ ! -f "$BLAST_DB.phr" ]; then
        log_error "BLAST database not found: $BLAST_DB"
        exit 1
    fi
    
    # Check protein directory
    if [ ! -d "$PROTEIN_DIR" ]; then
        log_error "Protein directory not found: $PROTEIN_DIR"
        exit 1
    fi
    
    log_success "All requirements validated"
}

# Generate variants from protein alignments
generate_variants() {
    local sample=$1
    local gene=$2
    local ref_protein=$3
    local sample_protein=$4
    local variants_file=$5
    
    log_info "Generating variants for $sample $gene..."
    
    # Python script to generate variants from aligned sequences
    python3 << EOF
import sys
from Bio import SeqIO
from Bio.Align import PairwiseAligner

def generate_variants(ref_file, sample_file, output_file):
    try:
        # Read sequences
        ref_seq = str(SeqIO.read(ref_file, "fasta").seq)
        sample_seq = str(SeqIO.read(sample_file, "fasta").seq)
        
        # Simple position-by-position comparison
        variants = []
        min_len = min(len(ref_seq), len(sample_seq))
        
        for i in range(min_len):
            ref_aa = ref_seq[i]
            sample_aa = sample_seq[i]
            
            # Skip gaps and unknown amino acids
            if ref_aa in ['-', 'X', '*'] or sample_aa in ['-', 'X', '*']:
                continue
            
            # Record substitutions
            if ref_aa != sample_aa:
                variants.append(f"{ref_aa}{i+1}{sample_aa}")
        
        # Handle length differences (indels)
        if len(sample_seq) < len(ref_seq):
            # Deletions
            for i in range(len(sample_seq), len(ref_seq)):
                if ref_seq[i] not in ['-', 'X', '*']:
                    variants.append(f"{ref_seq[i]}{i+1}del")
        
        # Write variants file
        with open(output_file, 'w') as f:
            for variant in variants:
                f.write(f"{variant}\n")
        
        return len(variants)
        
    except Exception as e:
        print(f"Error generating variants: {e}", file=sys.stderr)
        return 0

# Generate variants
variant_count = generate_variants("$ref_protein", "$sample_protein", "$variants_file")
print(f"Generated {variant_count} variants")
EOF

    local variant_count=$(tail -1 "$LOGS_DIR/provean_analysis.log" | grep -o "[0-9]\+" || echo "0")
    if [ "$variant_count" -gt 0 ]; then
        log_success "Generated $variant_count variants for $sample $gene"
        return 0
    else
        log_warn "No variants generated for $sample $gene"
        return 1
    fi
}

# Run PROVEAN analysis
run_provean() {
    local query_file=$1
    local variants_file=$2
    local output_file=$3
    local gene_name=$4
    
    log_info "Running PROVEAN for $gene_name..."
    
    # Check input files
    if [ ! -f "$query_file" ] || [ ! -f "$variants_file" ]; then
        log_error "Input files missing for PROVEAN"
        return 1
    fi
    
    # Check if variants file is empty
    if [ ! -s "$variants_file" ]; then
        log_warn "No variants to analyze for $gene_name"
        echo "No variants found" > "$output_file"
        return 0
    fi
    
    # Run PROVEAN with timeout
    local cmd="timeout 1800 $PROVEAN_BIN -q \"$query_file\" -v \"$variants_file\" --num_threads $THREADS"
    
    if eval "$cmd" > "$output_file" 2>&1; then
        # Check if PROVEAN completed successfully
        if grep -q "PROVEAN scores" "$output_file"; then
            log_success "PROVEAN completed for $gene_name"
            return 0
        else
            log_error "PROVEAN failed to generate scores for $gene_name"
            return 1
        fi
    else
        log_error "PROVEAN execution failed for $gene_name (timeout or error)"
        return 1
    fi
}

# Parse PROVEAN results
parse_provean_results() {
    local provean_output=$1
    local summary_file=$2
    local gene_name=$3
    local sample_name=$4
    
    if [ ! -f "$provean_output" ] || [ ! -s "$provean_output" ]; then
        echo "$gene_name,$sample_name,0,0,0,0.0" >> "$summary_file"
        return
    fi
    
    # Extract PROVEAN scores using awk
    awk -v gene="$gene_name" -v sample="$sample_name" -v threshold="$PROVEAN_THRESHOLD" '
    BEGIN {
        total = 0
        deleterious = 0
        neutral = 0
        sum_scores = 0
    }
    /^[A-Z][0-9]+[A-Z]/ || /^[A-Z][0-9]+del/ {
        if (NF >= 2 && $2 ~ /^-?[0-9]+\.?[0-9]*$/) {
            total++
            sum_scores += $2
            if ($2 <= threshold) {
                deleterious++
            } else {
                neutral++
            }
        }
    }
    END {
        avg_score = (total > 0) ? sum_scores / total : 0
        printf "%s,%s,%d,%d,%d,%.3f\n", gene, sample, total, deleterious, neutral, avg_score
    }' "$provean_output" >> "$summary_file"
}

# Process samples for PROVEAN analysis
process_samples() {
    log_info "Processing samples for PROVEAN analysis..."
    
    # Genes suitable for PROVEAN (smaller proteins)
    local genes=("PreC_core" "X_antigen")
    
    # Reference mapping
    declare -A ref_mapping=(
        ["PreC_core"]="HBV_D_PreC_core"
        ["X_antigen"]="HBV_D_X_antigen"
    )
    
    # Sample lists
    local samples=("B2" "B3" "B5" "B12" "B13" "B14" "B15" "M2" "M3" "M5" "M12" "M13" "M14" "M15")
    
    # Initialize summary file
    local summary_file="$FUNCTIONAL_DIR/provean_summary.csv"
    echo "Gene,Sample,Total_Variants,Deleterious,Neutral,Avg_Score" > "$summary_file"
    
    local total_analyses=0
    local successful_analyses=0
    
    # Process each gene and sample combination
    for gene in "${genes[@]}"; do
        local ref_file="$PROJECT_ROOT/results/protein_analysis/reference_proteins/${ref_mapping[$gene]}.fasta"
        
        # Check if reference exists
        if [ ! -f "$ref_file" ]; then
            log_warn "Reference protein not found: $ref_file"
            continue
        fi
        
        for sample in "${samples[@]}"; do
            local sample_file="$PROTEIN_DIR/clean_proteins/${sample}_${gene}_clean.fasta"
            
            # Check if sample protein exists
            if [ ! -f "$sample_file" ]; then
                log_warn "Sample protein not found: $sample_file"
                continue
            fi
            
            total_analyses=$((total_analyses + 1))
            
            # Generate variants
            local variants_file="$FUNCTIONAL_DIR/provean_results/${sample}_${gene}_variants.txt"
            local provean_output="$FUNCTIONAL_DIR/provean_results/${sample}_${gene}_provean.txt"
            
            log_info "[$total_analyses] Processing $sample $gene..."
            
            if generate_variants "$sample" "$gene" "$ref_file" "$sample_file" "$variants_file"; then
                # Run PROVEAN
                if run_provean "$ref_file" "$variants_file" "$provean_output" "${sample}_${gene}"; then
                    successful_analyses=$((successful_analyses + 1))
                    
                    # Parse results
                    parse_provean_results "$provean_output" "$summary_file" "$gene" "$sample"
                else
                    log_error "PROVEAN failed for $sample $gene"
                    echo "$gene,$sample,0,0,0,0.0" >> "$summary_file"
                fi
            else
                log_warn "No variants to analyze for $sample $gene"
                echo "$gene,$sample,0,0,0,0.0" >> "$summary_file"
            fi
        done
    done
    
    log_info "PROVEAN analysis completed: $successful_analyses/$total_analyses successful"
    return $((total_analyses - successful_analyses))
}

# Generate comprehensive report
generate_report() {
    log_info "Generating PROVEAN analysis report..."
    
    local report_file="$FUNCTIONAL_DIR/provean_report.txt"
    local summary_file="$FUNCTIONAL_DIR/provean_summary.csv"
    
    # Calculate overall statistics
    local total_variants=$(awk -F',' 'NR>1 {sum+=$3} END {print sum+0}' "$summary_file")
    local total_deleterious=$(awk -F',' 'NR>1 {sum+=$4} END {print sum+0}' "$summary_file")
    local total_neutral=$(awk -F',' 'NR>1 {sum+=$5} END {print sum+0}' "$summary_file")
    
    local deleterious_pct=0
    if [ "$total_variants" -gt 0 ]; then
        deleterious_pct=$(echo "scale=1; $total_deleterious * 100 / $total_variants" | bc -l)
    fi
    
    # Generate report
    cat > "$report_file" << EOF
HBV MTCT Analysis - PROVEAN Functional Impact Prediction Report
===============================================================

Analysis Date: $(date)
PROVEAN Version: $(provean.sh --version 2>/dev/null | head -1 || echo "Unknown")
PROVEAN Threshold: $PROVEAN_THRESHOLD (deleterious ≤ threshold)
Database: Swiss-Prot

OVERALL RESULTS
===============
Total Variants Analyzed: $total_variants
Deleterious Variants: $total_deleterious (${deleterious_pct}%)
Neutral/Tolerated Variants: $total_neutral

GENE-SPECIFIC RESULTS
=====================
EOF

    # Add gene-specific statistics
    for gene in "PreC_core" "X_antigen"; do
        local gene_variants=$(awk -F',' -v g="$gene" '$1==g {sum+=$3} END {print sum+0}' "$summary_file")
        local gene_deleterious=$(awk -F',' -v g="$gene" '$1==g {sum+=$4} END {print sum+0}' "$summary_file")
        local gene_pct=0
        
        if [ "$gene_variants" -gt 0 ]; then
            gene_pct=$(echo "scale=1; $gene_deleterious * 100 / $gene_variants" | bc -l)
        fi
        
        cat >> "$report_file" << EOF

$gene Protein:
  Total Variants: $gene_variants
  Deleterious: $gene_deleterious (${gene_pct}%)
  
EOF
    done
    
    # Add most deleterious mutations
    cat >> "$report_file" << EOF

CRITICAL FINDINGS
=================
Most Deleterious Mutations (Score ≤ -10):
EOF

    # Extract highly deleterious mutations
    awk -F',' 'NR>1 && $6 <= -10 {print "  " $1 " (" $2 "): " $4 " deleterious variants, avg score " $6}' "$summary_file" >> "$report_file"
    
    cat >> "$report_file" << EOF

CLINICAL INTERPRETATION
=======================
- High deleterious rate indicates strong functional constraints on HBV proteins
- PreC/Core mutations may affect viral replication and immune recognition
- X protein mutations may impact transcriptional regulation
- Scores ≤ $PROVEAN_THRESHOLD suggest significant functional impact

FILES GENERATED
===============
- Summary data: $summary_file
- Individual results: $FUNCTIONAL_DIR/provean_results/
- This report: $report_file

CITATION
========
Please cite: Choi, Y. & Chan, A.P. (2015). PROVEAN web server: a tool to predict 
the functional effect of amino acid substitutions and indels. Bioinformatics, 
31(16), 2745-2747.
EOF

    log_success "Report generated: $report_file"
}

# Main execution function
main() {
    log_info "Starting PROVEAN functional impact analysis"
    log_info "Configuration:"
    log_info "  PROVEAN threshold: $PROVEAN_THRESHOLD"
    log_info "  Threads: $THREADS"
    log_info "  Output directory: $FUNCTIONAL_DIR"
    log_info "  BLAST database: $BLAST_DB"
    
    # Validate requirements
    validate_requirements
    
    # Process samples
    local failed_analyses
    if ! failed_analyses=$(process_samples); then
        failed_analyses=0
    fi
    
    # Generate report
    generate_report
    
    # Final summary
    log_info "PROVEAN analysis completed"
    if [ "$failed_analyses" -eq 0 ]; then
        log_success "All analyses completed successfully"
        exit 0
    else
        log_warn "$failed_analyses analyses failed"
        exit 1
    fi
}

# Cleanup function
cleanup() {
    log_info "Cleaning up temporary files..."
    # Remove any temporary files if needed
}

trap cleanup EXIT

# Execute main function
main "$@"
