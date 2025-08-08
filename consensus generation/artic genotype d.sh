#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=480000
#SBATCH --job-name=artic_genotype_d
#SBATCH --output=logs/artic_genotype_d.out
#SBATCH --error=logs/artic_genotype_d.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@institution.edu

# HBV MTCT Analysis Pipeline - ARTIC Consensus Generation (Genotype D)
#
# Purpose: Generate consensus sequences using ARTIC pipeline for HBV genotype D samples
# Method: ARTIC minion pipeline with HBV-specific primer schemes
# Input: Concatenated FASTQ files from data preparation step
# Output: Consensus FASTA sequences with quality metrics
#
# Citation: Loman, N.J., et al. (2020). A user-friendly ARTIC pipeline for viral 
#           genome sequencing and real-time surveillance. bioRxiv.
# 
# Author: Your Name
# Date: 2024
# Version: 1.0

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Script directory and configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Load configuration
CONFIG_FILE="$PROJECT_ROOT/config/pipeline_config.yaml"
if [ -f "$CONFIG_FILE" ]; then
    INPUT_DIR=$(grep -A 10 "paths:" "$CONFIG_FILE" | grep "input_dir:" | awk '{print $2}' | tr -d '"')
    OUTPUT_DIR=$(grep -A 10 "paths:" "$CONFIG_FILE" | grep "output_dir:" | awk '{print $2}' | tr -d '"')
    THREADS=$(grep -A 10 "resources:" "$CONFIG_FILE" | grep "threads:" | awk '{print $2}')
    NORMALISE=$(grep -A 10 "consensus:" "$CONFIG_FILE" | grep "normalise_reads:" | awk '{print $2}')
else
    # Default values
    INPUT_DIR="data/raw_fastq"
    OUTPUT_DIR="results"
    THREADS=8
    NORMALISE=200
fi

# Directories
PRIMER_DIR="$PROJECT_ROOT/data/primer_schemes/HBV_D"
CONSENSUS_DIR="$OUTPUT_DIR/consensus"
LOGS_DIR="$OUTPUT_DIR/logs"

# Create output directories
mkdir -p "$CONSENSUS_DIR" "$LOGS_DIR"

# Logging functions
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOGS_DIR/artic_genotype_d.log"
}

log_info() {
    log "INFO: $1"
}

log_warn() {
    log "WARN: $1"
}

log_error() {
    log "ERROR: $1"
}

log_success() {
    log "SUCCESS: $1"
}

# Validate requirements
validate_requirements() {
    log_info "Validating requirements..."
    
    # Check if ARTIC environment is activated
    if [ -z "$CONDA_DEFAULT_ENV" ] || [ "$CONDA_DEFAULT_ENV" != "hbv_analysis" ]; then
        log_error "ARTIC conda environment not activated"
        exit 1
    fi
    
    # Check required tools
    local tools=("artic" "minimap2" "samtools" "bcftools" "medaka")
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            log_error "Required tool not found: $tool"
            exit 1
        fi
    done
    
    # Check primer scheme directory
    if [ ! -d "$PRIMER_DIR" ]; then
        log_error "Primer scheme directory not found: $PRIMER_DIR"
        exit 1
    fi
    
    # Check required primer scheme files
    local required_files=("reference.fasta" "scheme.bed" "primers.tsv")
    for file in "${required_files[@]}"; do
        if [ ! -f "$PRIMER_DIR/$file" ]; then
            log_error "Required primer scheme file not found: $PRIMER_DIR/$file"
            exit 1
        fi
    done
    
    log_success "All requirements validated"
}

# Sample mapping for genotype D
declare -A GENOTYPE_D_SAMPLES=(
    ["barcode01"]="B2"    # Baby 2
    ["barcode02"]="B5"    # Baby 5
    ["barcode03"]="B12"   # Baby 12
    ["barcode04"]="B13"   # Baby 13
    ["barcode05"]="M14"   # Mother 14
    ["barcode06"]="B14"   # Baby 14
    ["barcode07"]="M12"   # Mother 12
    ["barcode08"]="M13"   # Mother 13
    ["barcode83"]="M2"    # Mother 2
    ["barcode85"]="M5"    # Mother 5
)

# Function to run ARTIC for a single sample
run_artic_sample() {
    local barcode="$1"
    local sample="$2"
    local fastq_file="$INPUT_DIR/$barcode/$sample.fastq"
    local sample_output="$CONSENSUS_DIR/$sample"
    
    log_info "Processing sample $sample from $barcode"
    
    # Check if input file exists
    if [ ! -f "$fastq_file" ]; then
        log_error "Input FASTQ file not found: $fastq_file"
        return 1
    fi
    
    # Check file size and read count
    local file_size=$(du -h "$fastq_file" | cut -f1)
    local read_count=$(grep -c "^@" "$fastq_file" || echo "0")
    log_info "Input file: $file_size, $read_count reads"
    
    if [ "$read_count" -lt 1000 ]; then
        log_warn "Low read count for $sample: $read_count reads"
    fi
    
    # Create sample output directory
    mkdir -p "$sample_output"
    
    # Change to output directory (required by ARTIC)
    cd "$sample_output"
    
    log_info "Running ARTIC minion pipeline..."
    
    # Run ARTIC minion pipeline
    if artic minion \
        --normalise "$NORMALISE" \
        --threads "$THREADS" \
        --bed "$PRIMER_DIR/scheme.bed" \
        --ref "$PRIMER_DIR/reference.fasta" \
        --read-file "$fastq_file" \
        "$sample" 2>&1 | tee "$LOGS_DIR/${sample}_artic.log"; then
        
        log_success "ARTIC completed for $sample"
        
        # Validate output files
        validate_artic_output "$sample" "$sample_output"
        
    else
        log_error "ARTIC failed for $sample"
        cd "$PROJECT_ROOT"
        return 1
    fi
    
    # Return to project root
    cd "$PROJECT_ROOT"
    
    return 0
}

# Function to validate ARTIC output
validate_artic_output() {
    local sample="$1"
    local output_dir="$2"
    
    log_info "Validating ARTIC output for $sample..."
    
    # Expected output files
    local expected_files=(
        "$sample.consensus.fasta"
        "$sample.pass.vcf.gz"
        "$sample.fail.vcf.gz" 
        "$sample.primertrimmed.rg.sorted.bam"
        "$sample.coverage_mask.txt"
    )
    
    local missing_files=0
    
    for file in "${expected_files[@]}"; do
        if [ ! -f "$output_dir/$file" ]; then
            log_warn "Missing output file: $file"
            missing_files=$((missing_files + 1))
        fi
    done
    
    # Check consensus sequence quality
    if [ -f "$output_dir/$sample.consensus.fasta" ]; then
        local consensus_length=$(grep -v "^>" "$output_dir/$sample.consensus.fasta" | tr -d '\n' | wc -c)
        local n_count=$(grep -v "^>" "$output_dir/$sample.consensus.fasta" | tr -d '\n' | grep -o "N" | wc -l)
        local coverage_pct=$(echo "scale=2; ($consensus_length - $n_count) * 100 / $consensus_length" | bc -l)
        
        log_info "Consensus statistics:"
        log_info "  Length: $consensus_length bp"
        log_info "  N bases: $n_count"
        log_info "  Coverage: $coverage_pct%"
        
        # Quality assessment
        if (( $(echo "$coverage_pct >= 90" | bc -l) )); then
            log_success "High quality consensus (≥90% coverage)"
        elif (( $(echo "$coverage_pct >= 70" | bc -l) )); then
            log_warn "Moderate quality consensus (70-90% coverage)"
        else
            log_warn "Low quality consensus (<70% coverage)"
        fi
    else
        log_error "Consensus file not generated"
        missing_files=$((missing_files + 1))
    fi
    
    if [ $missing_files -eq 0 ]; then
        log_success "All expected output files present"
    else
        log_warn "$missing_files output files missing"
    fi
}

# Function to generate summary report
generate_summary() {
    log_info "Generating summary report..."
    
    local summary_file="$CONSENSUS_DIR/genotype_d_summary.txt"
    
    cat > "$summary_file" << EOF
HBV Genotype D - ARTIC Consensus Generation Summary
===================================================
Date: $(date)
Pipeline: ARTIC minion
Normalisation: $NORMALISE reads
Threads: $THREADS
Primer scheme: HBV_D

Sample Results:
EOF

    # Process each sample
    local total_samples=0
    local successful_samples=0
    
    for barcode in "${!GENOTYPE_D_SAMPLES[@]}"; do
        local sample="${GENOTYPE_D_SAMPLES[$barcode]}"
        local sample_dir="$CONSENSUS_DIR/$sample"
        local consensus_file="$sample_dir/$sample.consensus.fasta"
        
        total_samples=$((total_samples + 1))
        
        if [ -f "$consensus_file" ]; then
            local length=$(grep -v "^>" "$consensus_file" | tr -d '\n' | wc -c)
            local n_count=$(grep -v "^>" "$consensus_file" | tr -d '\n' | grep -o "N" | wc -l)
            local coverage_pct=$(echo "scale=1; ($length - $n_count) * 100 / $length" | bc -l)
            
            echo "$sample ($barcode): $length bp, ${coverage_pct}% coverage, $n_count Ns" >> "$summary_file"
            successful_samples=$((successful_samples + 1))
        else
            echo "$sample ($barcode): FAILED" >> "$summary_file"
        fi
    done
    
    # Add summary statistics
    cat >> "$summary_file" << EOF

Summary Statistics:
Total samples processed: $total_samples
Successful: $successful_samples
Failed: $((total_samples - successful_samples))
Success rate: $(echo "scale=1; $successful_samples * 100 / $total_samples" | bc -l)%

Output Location: $CONSENSUS_DIR
Log Files: $LOGS_DIR
EOF

    log_info "Summary report: $summary_file"
}

# Main execution
main() {
    log_info "Starting ARTIC consensus generation for HBV genotype D samples"
    log_info "Configuration:"
    log_info "  Input directory: $INPUT_DIR"
    log_info "  Output directory: $OUTPUT_DIR"
    log_info "  Primer scheme: $PRIMER_DIR"
    log_info "  Threads: $THREADS"
    log_info "  Normalisation: $NORMALISE reads"
    
    # Validate requirements
    validate_requirements
    
    # Initialize counters
    local total_samples=${#GENOTYPE_D_SAMPLES[@]}
    local processed=0
    local successful=0
    local failed=0
    
    log_info "Processing $total_samples genotype D samples"
    
    # Process each sample
    for barcode in "${!GENOTYPE_D_SAMPLES[@]}"; do
        sample="${GENOTYPE_D_SAMPLES[$barcode]}"
        processed=$((processed + 1))
        
        log_info "[$processed/$total_samples] Processing $sample"
        
        if run_artic_sample "$barcode" "$sample"; then
            successful=$((successful + 1))
            log_success "Sample $sample completed successfully"
        else
            failed=$((failed + 1))
            log_error "Sample $sample failed"
        fi
    done
    
    # Generate summary
    generate_summary
    
    # Final report
    log_info "ARTIC consensus generation completed"
    log_info "Results:"
    log_info "  Total samples: $total_samples"
    log_info "  Successful: $successful"
    log_info "  Failed: $failed"
    log_info "  Success rate: $(echo "scale=1; $successful * 100 / $total_samples" | bc -l)%"
    
    if [ $failed -eq 0 ]; then
        log_success "All samples processed successfully"
        exit 0
    else
        log_warn "$failed samples failed - check individual logs"
        exit 1
    fi
}

# Trap for cleanup
cleanup() {
    log_info "Cleaning up temporary files..."
    # Add any cleanup commands here
}

trap cleanup EXIT

# Execute main function
main "$@"
