#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=100000
#SBATCH --job-name=fastq_concatenate
#SBATCH --output=logs/fastq_concatenate.out
#SBATCH --error=logs/fastq_concatenate.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@institution.edu

# HBV MTCT Analysis Pipeline - FASTQ File Concatenation
# 
# Purpose: Concatenate FASTQ files from Oxford Nanopore sequencing by barcode
# Input: Raw FASTQ files organized by barcode directories
# Output: Single concatenated FASTQ file per sample
#
# Author: Your Name
# Date: 2024
# Version: 1.0

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Source configuration if available
if [ -f "$PROJECT_ROOT/config/pipeline_config.yaml" ]; then
    # Note: This is a simplified config reader for bash
    BASE_DIR=$(grep -A 10 "paths:" "$PROJECT_ROOT/config/pipeline_config.yaml" | grep "input_dir:" | awk '{print $2}' | tr -d '"')
    BASE_DIR=${BASE_DIR:-"data/raw_fastq"}
else
    BASE_DIR="data/raw_fastq"
fi

# Logging functions
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
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

# Sample mapping array
# Maps barcode directories to sample IDs
declare -A barcode_to_sample=(
    ["barcode01"]="B2"    # Baby 2
    ["barcode02"]="B5"    # Baby 5
    ["barcode03"]="B12"   # Baby 12
    ["barcode04"]="B13"   # Baby 13
    ["barcode05"]="M14"   # Mother 14
    ["barcode06"]="B14"   # Baby 14
    ["barcode07"]="M12"   # Mother 12
    ["barcode08"]="M13"   # Mother 13
    ["barcode09"]="B3"    # Baby 3
    ["barcode10"]="B15"   # Baby 15
    ["barcode11"]="M15"   # Mother 15
    ["barcode83"]="M2"    # Mother 2
    ["barcode84"]="M3"    # Mother 3
    ["barcode85"]="M5"    # Mother 5
)

# Validate input directory
if [ ! -d "$BASE_DIR" ]; then
    log_error "Input directory does not exist: $BASE_DIR"
    exit 1
fi

log_info "Starting FASTQ file concatenation"
log_info "Base directory: $BASE_DIR"
log_info "Processing ${#barcode_to_sample[@]} barcodes"

# Initialize counters
processed=0
successful=0
failed=0

# Process each barcode
for barcode in "${!barcode_to_sample[@]}"; do
    sample_id="${barcode_to_sample[$barcode]}"
    barcode_dir="$BASE_DIR/$barcode"
    output_file="$barcode_dir/$sample_id.fastq"
    
    log_info "Processing $barcode -> $sample_id"
    processed=$((processed + 1))
    
    # Check if barcode directory exists
    if [ ! -d "$barcode_dir" ]; then
        log_error "Directory does not exist: $barcode_dir"
        failed=$((failed + 1))
        continue
    fi
    
    # Check for FASTQ files
    fastq_files=("$barcode_dir"/*.fastq.gz)
    if [ ! -e "${fastq_files[0]}" ]; then
        log_warn "No .fastq.gz files found in $barcode_dir"
        failed=$((failed + 1))
        continue
    fi
    
    # Count input files
    file_count=$(ls "$barcode_dir"/*.fastq.gz 2>/dev/null | wc -l)
    log_info "  Found $file_count FASTQ files to concatenate"
    
    # Remove existing output file if present
    if [ -f "$output_file" ]; then
        log_info "  Removing existing output file: $output_file"
        rm "$output_file"
    fi
    
    # Concatenate files
    log_info "  Concatenating files..."
    if zcat "$barcode_dir"/*.fastq.gz > "$output_file"; then
        # Verify output file was created and is not empty
        if [ -s "$output_file" ]; then
            read_count=$(grep -c "^@" "$output_file" || true)
            file_size=$(du -h "$output_file" | cut -f1)
            log_info "  ✅ Success: $output_file ($read_count reads, $file_size)"
            successful=$((successful + 1))
        else
            log_error "  Output file is empty: $output_file"
            rm -f "$output_file"
            failed=$((failed + 1))
        fi
    else
        log_error "  Failed to concatenate files for $barcode"
        rm -f "$output_file"
        failed=$((failed + 1))
    fi
done

# Summary
log_info "Concatenation completed"
log_info "Summary:"
log_info "  Processed: $processed barcodes"
log_info "  Successful: $successful concatenations"
log_info "  Failed: $failed concatenations"

# Create summary report
summary_file="$BASE_DIR/concatenation_summary.txt"
cat > "$summary_file" << EOF
FASTQ Concatenation Summary
===========================
Date: $(date)
Base Directory: $BASE_DIR
Total Barcodes: $processed
Successful: $successful
Failed: $failed

Sample Details:
EOF

for barcode in "${!barcode_to_sample[@]}"; do
    sample_id="${barcode_to_sample[$barcode]}"
    output_file="$BASE_DIR/$barcode/$sample_id.fastq"
    
    if [ -f "$output_file" ]; then
        read_count=$(grep -c "^@" "$output_file" || echo "0")
        file_size=$(du -h "$output_file" | cut -f1)
        echo "$barcode -> $sample_id: $read_count reads, $file_size" >> "$summary_file"
    else
        echo "$barcode -> $sample_id: FAILED" >> "$summary_file"
    fi
done

log_info "Summary report: $summary_file"

# Exit with appropriate code
if [ $failed -eq 0 ]; then
    log_info "All concatenations completed successfully"
    exit 0
else
    log_error "$failed concatenations failed"
    exit 1
fi
