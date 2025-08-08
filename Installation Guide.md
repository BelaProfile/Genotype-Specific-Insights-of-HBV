# Installation Guide

## System Requirements

### Hardware
- **CPU**: 8+ cores recommended (minimum 4 cores)
- **RAM**: 64GB recommended (minimum 32GB)
- **Storage**: 500GB available space
- **OS**: Linux (Ubuntu 18.04+, CentOS 7+) or macOS 10.15+

### Software Prerequisites
- **Conda/Miniconda**: 4.10+
- **Git**: 2.20+
- **Python**: 3.8+ (installed via conda)
- **R**: 4.0+ (installed via conda)

## Installation Methods

### Method 1: Automated Installation (Recommended)

```bash
# Clone repository
git clone https://github.com/yourusername/hbv-mtct-analysis.git
cd hbv-mtct-analysis

# Make setup script executable
chmod +x setup/install_dependencies.sh

# Run automated installation
./setup/install_dependencies.sh

# Activate environment
conda activate hbv_analysis

# Verify installation
./setup/verify_installation.sh
```

### Method 2: Manual Installation

#### Step 1: Create Conda Environment
```bash
# Create base environment
conda create -n hbv_analysis python=3.8 -y
conda activate hbv_analysis
```

#### Step 2: Install Bioinformatics Tools
```bash
# Core tools
conda install -c bioconda -c conda-forge \
    artic=1.2.1 \
    medaka=1.4.4 \
    minimap2=2.24 \
    samtools=1.15 \
    bcftools=1.15 \
    bedtools=2.30.0 \
    emboss=6.6.0 \
    blast=2.12.0 \
    mafft=7.490 \
    fasttree=2.1.11 \
    lofreq=2.1.5 \
    seqkit=2.3.0 \
    -y

# Additional Python packages
pip install pandas numpy scipy matplotlib seaborn biopython

# R and packages
conda install -c conda-forge r-base=4.2.0 r-essentials -y
```

#### Step 3: Install R Packages
```bash
R -e "
packages <- c('ggplot2', 'dplyr', 'tidyr', 'viridis', 'gridExtra', 
              'pheatmap', 'ggrepel', 'RColorBrewer', 'vegan', 'Biostrings')
install.packages(packages, repos='https://cran.r-project.org/')

# Bioconductor packages
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install('Biostrings')
"
```

#### Step 4: Install PROVEAN (Manual)
```bash
# Download and install PROVEAN
mkdir -p tools/provean
cd tools/provean
wget http://provean.jcvi.org/genome_submit_2/provean-1.1.5.tar.gz
tar -xzf provean-1.1.5.tar.gz
cd provean-1.1.5
./configure --prefix=$CONDA_PREFIX
make && make install
cd ../../../
```

#### Step 5: Install SIFT4G
```bash
conda install -c bioconda sift4g=2.0.0 -y
```

#### Step 6: Setup AlphaFold (Optional)
```bash
# Install AlphaFold dependencies
pip install alphafold-colabfold jax jaxlib

# Note: Full AlphaFold installation requires additional setup
# See: https://github.com/deepmind/alphafold
```

## Database Setup

### BLAST Databases
```bash
# Create database directory
mkdir -p data/databases

# Download Swiss-Prot database
cd data/databases
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# Build BLAST database
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out swissprot -title "Swiss-Prot"
cd ../../
```

### Reference Genomes
```bash
# Download HBV reference genomes
mkdir -p data/reference_genomes
cd data/reference_genomes

# HBV genotype C reference
wget "https://www.ncbi.nlm.nih.gov/nuccore/AB033554.1?report=fasta&format=text" -O HBV_C.reference.fasta

# HBV genotype D reference  
wget "https://www.ncbi.nlm.nih.gov/nuccore/X65259.1?report=fasta&format=text" -O HBV_D.reference.fasta

cd ../../
```

## Verification

### Test Installation
```bash
# Run verification script
./setup/verify_installation.sh

# Expected output:
# ✅ Conda environment: OK
# ✅ Python packages: OK
# ✅ R packages: OK
# ✅ Bioinformatics tools: OK
# ✅ Databases: OK
```

### Quick Test Run
```bash
# Test with provided example data
./run_pipeline.sh --input data/example --output test_results --test

# Should complete in ~10 minutes
```

## Troubleshooting

### Common Issues

#### 1. Conda Environment Issues
```bash
# If conda commands fail
conda clean --all
conda update conda
```

#### 2. R Package Installation Failures
```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install libxml2-dev libcurl4-openssl-dev libssl-dev

# For CentOS/RHEL
sudo yum install libxml2-devel libcurl-devel openssl-devel
```

#### 3. Memory Issues
```bash
# Increase swap space (Ubuntu)
sudo fallocate -l 32G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

#### 4. PROVEAN Installation Issues
```bash
# Install development tools
sudo apt-get install build-essential  # Ubuntu
sudo yum groupinstall "Development Tools"  # CentOS

# Ensure BLAST+ is properly installed
which makeblastdb
which blastp
```

### Getting Help

If you encounter issues:

1. **Check logs**: All installation logs are saved in `setup/logs/`
2. **GitHub Issues**: [Report bugs](https://github.com/yourusername/hbv-mtct-analysis/issues)
3. **Email Support**: your.email@institution.edu

### Version Information

To check installed versions:
```bash
# Show all tool versions
./setup/show_versions.sh

# Expected output includes:
# Python: 3.8.x
# R: 4.2.x  
# ARTIC: 1.2.1
# Medaka: 1.4.4
# etc.
```

## Alternative Installation Methods

### Docker Installation (Experimental)
```bash
# Build Docker image
docker build -t hbv-mtct-analysis .

# Run container
docker run -v $(pwd)/data:/data -v $(pwd)/results:/results hbv-mtct-analysis
```

### Singularity Installation
```bash
# Build Singularity image
singularity build hbv-mtct-analysis.sif docker://yourusername/hbv-mtct-analysis

# Run analysis
singularity exec hbv-mtct-analysis.sif ./run_pipeline.sh
```

## Performance Optimization

### Resource Configuration
Edit `config/pipeline_config.yaml`:
```yaml
resources:
  threads: 16        # Adjust based on available cores
  memory: "64G"      # Adjust based on available RAM
  tmp_dir: "/tmp"    # Use fast storage for temporary files
```

### HPC/Cluster Setup
For SLURM clusters, modify the SBATCH headers in scripts:
```bash
#SBATCH --partition=your_partition
#SBATCH --account=your_account
#SBATCH --qos=your_qos
```

---

**Next Steps**: After successful installation, proceed to the [Usage Guide](USAGE.md).
