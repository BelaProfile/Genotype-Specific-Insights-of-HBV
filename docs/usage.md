# Usage Guide

This guide provides detailed instructions for using the HBV MTCT Analysis Pipeline.

## Quick Start

### 1. Prerequisites
- Conda environment `hbv_analysis` activated
- Input FASTQ files in `data/raw_fastq/`
- Reference genomes and databases installed

### 2. Basic Usage
```bash
# Activate environment
conda activate hbv_analysis

# Run complete pipeline
./run_pipeline.sh --input data/raw_fastq --output results --threads 16

# Check results
ls results/
```

## Detailed Usage

### Command Line Options

```bash
./run_pipeline.sh [OPTIONS]

Required:
  -i, --input DIR         Input directory with FASTQ files
  -o, --output DIR        Output directory (default: results)

Optional:
  -c, --config FILE       Configuration file (default: config/pipeline_config.yaml)
  -t, --threads INT       Number of threads (default: 8)
  -m, --memory STRING     Memory limit (default: 32G)
  -M, --module STRING     Run specific module only
  -r, --resume            Resume from previous run
  --test                  Run with test data
  --dry-run               Show commands without executing
  -v, --verbose           Verbose output
  -h, --help              Show help message
```

### Input Data Organization

Your input directory should be organized as follows:
```
data/raw_fastq/
├── barcode01/
│   ├── file1.fastq.gz
│   ├── file2.fastq.gz
│   └── ...
├── barcode02/
│   └── ...
└── barcode83/
    └── ...
```

### Module-Specific Execution

#### Data Preparation
```bash
./run_pipeline.sh -i data/raw_fastq -M data_preparation
```
- Concatenates FASTQ files by barcode
- Creates sample metadata
- Sets up primer schemes

#### Consensus Generation
```bash
./run_pipeline.sh -i data/raw_fastq -M consensus_generation
```
- Runs ARTIC pipeline for both genotypes
- Generates consensus FASTA sequences
- Provides quality metrics

#### Quality Control
```bash
./run_pipeline.sh -i data/raw_fastq -M quality_control
```
- Compares consensus methods
- Generates QC statistics
- Validates sequence quality

#### Phylogenetic Analysis
```bash
./run_pipeline.sh -i data/raw_fastq -M phylogenetic_analysis
```
- Builds combined phylogenetic trees
- Analyzes mother-baby relationships
- Performs evolutionary analysis

#### Gene Extraction
```bash
./run_pipeline.sh -i data/raw_fastq -M gene_extraction
```
- Extracts HBV genes from consensus
- Translates to protein sequences
- Handles circular genome topology

#### Protein Analysis
```bash
./run_pipeline.sh -i data/raw_fastq -M protein_analysis
```
- Cleans and validates ORFs
- Compares mother-baby proteins
- Prepares for structural analysis

#### Functional Prediction
```bash
./run_pipeline.sh -i data/raw_fastq -M functional_prediction
```
- Runs PROVEAN and SIFT analysis
- Predicts mutation functional impact
- Identifies deleterious variants

#### Epitope Mapping
```bash
./run_pipeline.sh -i data/raw_fastq -M epitope_mapping
```
- Maps B-cell and T-cell epitopes
- Analyzes epitope disruption
- Compares with IEDB database

#### Quasispecies Analysis
```bash
./run_pipeline.sh -i data/raw_fastq -M quasispecies_analysis
```
- Calls variants with LoFreq
- Calculates Shannon diversity
- Analyzes transmission bottlenecks

#### Visualization
```bash
./run_pipeline.sh -i data/raw_fastq -M visualization
```
- Generates all plots and figures
- Creates publication-ready images
- Produces interactive dashboard

## Configuration

### Pipeline Configuration File

Edit `config/pipeline_config.yaml`:

```yaml
# Resource allocation
resources:
  threads: 16              # Adjust based on your system
  memory: "64G"           # Adjust based on available RAM

# Analysis parameters
analysis:
  consensus:
    min_coverage: 20
    normalise_reads: 200
  
  functional:
    provean_threshold: -2.5
    sift_threshold: 0.05

# Sample information (customize for your data)
samples:
  genotype_c: ["B3", "B15", "M3", "M15"]
  genotype_d: ["B2", "B5", "B12", "B13", "B14", "M2", "M5", "M12", "M13", "M14"]
```

### HPC/Cluster Usage

For SLURM clusters, modify the SBATCH headers in individual scripts:

```bash
#SBATCH --partition=your_partition
#SBATCH --account=your_account
#SBATCH --qos=your_qos
#SBATCH --time=24:00:00
```

Or use the cluster configuration file:

```yaml
# config/cluster_config.yaml
cluster:
  scheduler: slurm
  partition: defq
  account: your_account
  max_time: "24:00:00"
  max_memory: "64G"
```

## Advanced Usage

### Custom Analysis Workflows

#### Analyzing Specific Samples
```bash
# Edit sample list in config file, then run
./run_pipeline.sh -i data/raw_fastq -c config/custom_config.yaml
```

#### Resume Interrupted Analysis
```bash
./run_pipeline.sh -i data/raw_fastq -o results --resume
```

#### Test Run with Example Data
```bash
./run_pipeline.sh --test
```

#### Dry Run (Show Commands Only)
```bash
./run_pipeline.sh -i data/raw_fastq --dry-run
```

### Custom Functional Analysis

#### PROVEAN Analysis Only
```bash
cd 07_functional_prediction/
./provean_analysis.sh
```

#### Custom Epitope Analysis
```bash
cd 08_epitope_mapping/
python epitope_analysis.py --input ../results/protein_analysis/ --output custom_epitopes/
```

### Data Import/Export

#### Import Reference Genomes
```bash
python scripts/download_references.py --genotypes C,D --output data/reference_genomes/
```

#### Export Results
```bash
python scripts/merge_results.py --input results/ --output combined_results.xlsx
```

#### Backup Important Results
```bash
./scripts/backup_results.sh results/ backup_$(date +%Y%m%d)/
```

## Troubleshooting

### Common Issues

#### Memory Issues
```bash
# Reduce threads and increase memory
./run_pipeline.sh -i data/raw_fastq -t 4 -m 64G
```

#### Missing Dependencies
```bash
# Verify installation
./setup/verify_installation.sh

# Reinstall if needed
./setup/install_dependencies.sh
```

#### Low Quality Consensus
- Check input read quality
- Adjust normalization parameters
- Verify primer scheme compatibility

#### Failed Functional Prediction
- Check protein sequence quality
- Verify BLAST database integrity
- Reduce timeout limits for large proteins

### Log Analysis

Check logs for detailed error information:
```bash
# Main pipeline log
tail -f results/logs/pipeline_*.log

# Module-specific logs
tail -f results/logs/artic_genotype_d.log
tail -f results/logs/provean_analysis.log
```

### Performance Optimization

#### For Large Datasets
```bash
# Use more threads and memory
./run_pipeline.sh -i data/raw_fastq -t 32 -m 128G

# Process in batches
./run_pipeline.sh -i data/batch1/ -o results/batch1/
./run_pipeline.sh -i data/batch2/ -o results/batch2/
```

#### For Limited Resources
```bash
# Reduce resource usage
./run_pipeline.sh -i data/raw_fastq -t 2 -m 16G

# Run modules separately
./run_pipeline.sh -i data/raw_fastq -M consensus_generation
./run_pipeline.sh -i data/raw_fastq -M functional_prediction
```

## Output Interpretation

### Key Result Files

| File | Description |
|------|-------------|
| `results/consensus/*.consensus.fasta` | Consensus genome sequences |
| `results/phylogeny/*.nwk` | Phylogenetic trees |
| `results/functional_analysis/provean_summary.csv` | Functional predictions |
| `results/epitope_mapping/epitope_summary.csv` | Epitope analysis |
| `results/plots/*.png` | All generated visualizations |
| `results/HBV_MTCT_Analysis_Report.html` | Final comprehensive report |

### Quality Metrics

#### Consensus Quality
- **Coverage**: >90% = high quality, 70-90% = moderate, <70% = low
- **N count**: Lower is better (indicates uncertain bases)
- **Length**: Should be ~3200 bp for complete HBV genome

#### Functional Predictions
- **PROVEAN**: ≤-2.5 = deleterious, >-2.5 = neutral
- **SIFT**: ≤0.05 = deleterious, >0.05 = tolerated

#### Epitope Disruption
- **High disruption** (>70%): Major immune impact
- **Moderate disruption** (30-70%): Some immune impact  
- **Low disruption** (<30%): Minimal immune impact

## Integration with Other Tools

### Export to External Analysis

#### GISAID Submission
```bash
# Prepare sequences for GISAID
python scripts/prepare_gisaid_submission.py results/consensus/
```

#### Phylogenetic Analysis (Beast, MrBayes)
```bash
# Convert to NEXUS format
python scripts/convert_alignment.py results/phylogeny/combined_alignment.fasta --format nexus
```

#### Structural Analysis (PyMOL, ChimeraX)
```bash
# Export AlphaFold structures
cp results/protein_analysis/alphafold_structures/*.pdb /path/to/structural/analysis/
```

### Database Integration

#### Update Local Databases
```bash
./setup/setup_databases.sh --update
```

#### Custom Reference Addition
```bash
# Add new reference genome
cp new_reference.fasta data/reference_genomes/
# Update configuration
vim config/pipeline_config.yaml
```

## Reproducibility

### Version Control
```bash
# Record exact versions used
./setup/show_versions.sh > results/versions_used.txt

# Git commit analysis state
git add results/
git commit -m "HBV MTCT analysis results $(date)"
```

### Docker/Singularity Usage
```bash
# Docker execution
docker run -v $(pwd):/workspace hbv-mtct-analysis:latest

# Singularity execution  
singularity exec hbv-mtct-analysis.sif ./run_pipeline.sh
```

### Reproducible Environments
```bash
# Export conda environment
conda env export > environment_exact.yml

# Export pip requirements
pip freeze > requirements_exact.txt
```

## Getting Help

### Built-in Help
```bash
./run_pipeline.sh --help
python scripts/validate_inputs.py --help
Rscript 10_visualization/hbv_visualization.R --help
```

### Documentation
- **Installation**: [docs/INSTALLATION.md](INSTALLATION.md)
- **Troubleshooting**: [docs/TROUBLESHOOTING.md](TROUBLESHOOTING.md)
- **Citations**: [docs/CITATIONS.md](CITATIONS.md)

### Community Support
- **GitHub Issues**: [Report bugs or ask questions](https://github.com/yourusername/hbv-mtct-analysis/issues)
- **Discussions**: [Community discussions](https://github.com/yourusername/hbv-mtct-analysis/discussions)
- **Email**: your.email@institution.edu

---

**Next Steps**: After completing analysis, see the [Results Interpretation Guide](tutorials/interpretation_guide.md) for help understanding your results.
