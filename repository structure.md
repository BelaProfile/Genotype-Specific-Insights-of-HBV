# Complete HBV MTCT Analysis Repository Structure

This document provides the complete file structure for the GitHub repository with all necessary scripts, documentation, and configuration files.

## Repository Structure

```
hbv-mtct-analysis/
├── README.md                           # Main project documentation
├── LICENSE                             # MIT License
├── .gitignore                          # Git ignore file
├── run_pipeline.sh                     # Main pipeline execution script
├── 01_data_preparation/
│   ├── README.md                       # Module documentation
│   ├── concatenate_fastq.sh           # Concatenate FASTQ files by barcode
│   ├── create_samples_tsv.sh          # Create sample metadata
│   └── setup_primer_schemes.sh        # Setup ARTIC primer schemes
├── 02_consensus_generation/
│   ├── README.md                       # Module documentation
│   ├── artic_genotype_c.sh            # ARTIC pipeline for genotype C
│   ├── artic_genotype_d.sh            # ARTIC pipeline for genotype D
│   └── medaka_consensus.sh            # Alternative Medaka consensus
├── 03_quality_control/
│   ├── README.md                       # Module documentation
│   ├── compare_consensus.R             # Compare ARTIC vs Medaka
│   └── consensus_stats.sh             # Generate QC statistics
├── 04_phylogenetic_analysis/
│   ├── README.md                       # Module documentation
│   ├── combined_tree.sh               # Combined phylogenetic tree
│   ├── pairwise_trees.sh              # Mother-baby pairwise trees
│   └── evolutionary_analysis.sh       # Conservation and divergence analysis
├── 05_gene_extraction/
│   ├── README.md                       # Module documentation
│   ├── extract_genes_translate_proteins.sh  # Extract all HBV genes
│   ├── extract_translate_ref.sh       # Process reference sequences
│   └── extract_polymerase.sh          # Handle circular polymerase
├── 06_protein_analysis/
│   ├── README.md                       # Module documentation
│   ├── clean_orfs.sh                  # Clean and validate ORFs
│   ├── compare_proteins.sh            # Mother-baby protein comparison
│   ├── alphafold_preparation.sh       # Prepare for AlphaFold
│   └── mutation_analysis.py           # Python mutation analysis
├── 07_functional_prediction/
│   ├── README.md                       # Module documentation
│   ├── provean_analysis.sh            # PROVEAN functional prediction
│   ├── sift_analysis.sh               # SIFT functional prediction
│   └── combined_functional_analysis.py # Combine predictions
├── 08_epitope_mapping/
│   ├── README.md                       # Module documentation
│   ├── epitope_analysis.py            # Epitope disruption analysis
│   ├── process_bepipred_results.R     # Process BepiPred output
│   └── iedb_comparison.py             # Compare with IEDB database
├── 09_quasispecies_analysis/
│   ├── README.md                       # Module documentation
│   ├── lofreq_variant_calling.sh      # LoFreq variant calling
│   ├── shannon_diversity.R            # Shannon diversity analysis
│   ├── allele_frequency_analysis.py   # AF comparison analysis
│   └── transmission_bottleneck.R      # Bottleneck analysis
├── 10_visualization/
│   ├── README.md                       # Module documentation
│   ├── hbv_visualization.R             # Main visualization script
│   ├── create_publication_figures.R    # Publication-ready figures
│   ├── plot_phylogeny.R               # Phylogenetic tree plots
│   ├── plot_epitope_maps.R            # Epitope mapping plots
│   └── interactive_dashboard.R        # Interactive Shiny dashboard
├── config/
│   ├── pipeline_config.yaml           # Main pipeline configuration
│   ├── hbv_gene_coordinates.yaml      # HBV gene coordinate definitions
│   ├── sample_metadata.yaml           # Sample information
│   └── cluster_config.yaml            # HPC cluster configuration
├── data/
│   ├── README.md                       # Data directory documentation
│   ├── raw_fastq/                     # Raw sequencing data (user provided)
│   ├── processed/                     # Intermediate processed files
│   ├── reference_genomes/             # HBV reference sequences
│   │   ├── HBV_C.reference.fasta
│   │   ├── HBV_D.reference.fasta
│   │   ├── HBV_C.scheme.bed
│   │   └── HBV_D.scheme.bed
│   ├── databases/                     # External databases
│   │   ├── swissprot/                 # Swiss-Prot BLAST database
│   │   └── iedb_epitopes/             # IEDB epitope data
│   ├── primer_schemes/                # ARTIC primer schemes
│   │   ├── HBV_C/
│   │   └── HBV_D/
│   └── example/                       # Test data
│       └── barcode01/
├── results/
│   ├── README.md                       # Results documentation
│   ├── consensus/                      # Consensus sequences
│   ├── phylogeny/                      # Phylogenetic trees
│   ├── functional_analysis/           # Functional predictions
│   ├── epitope_mapping/               # Epitope analysis results
│   ├── quasispecies/                  # Diversity analysis
│   ├── plots/                         # All generated plots
│   ├── reports/                       # Analysis reports
│   └── logs/                          # Execution logs
├── docs/
│   ├── INSTALLATION.md                # Detailed installation guide
│   ├── USAGE.md                       # Usage documentation
│   ├── CITATIONS.md                   # Software citations
│   ├── TROUBLESHOOTING.md             # Common issues and solutions
│   ├── API.md                         # Pipeline API documentation
│   ├── CHANGELOG.md                   # Version history
│   └── tutorials/
│       ├── quickstart.md              # Quick start tutorial
│       ├── advanced_usage.md          # Advanced usage examples
│       └── interpretation_guide.md    # Results interpretation
├── setup/
│   ├── install_dependencies.sh        # Automated installation script
│   ├── verify_installation.sh         # Installation verification
│   ├── show_versions.sh               # Show tool versions
│   ├── setup_databases.sh             # Download and setup databases
│   ├── create_test_data.sh            # Generate test dataset
│   └── logs/                          # Installation logs
├── utils/
│   ├── __init__.py                    # Python package initialization
│   ├── file_utils.py                  # File handling utilities
│   ├── sequence_utils.py              # Sequence processing utilities
│   ├── plot_utils.py                  # Plotting helper functions
│   ├── config_parser.py               # Configuration file parser
│   ├── logging_utils.py               # Logging utilities
│   ├── validation.py                  # Input validation functions
│   └── cluster_utils.py               # HPC cluster utilities
├── tests/
│   ├── __init__.py                    # Test package initialization
│   ├── test_data_preparation.py       # Test data preparation module
│   ├── test_consensus_generation.py   # Test consensus generation
│   ├── test_phylogenetic_analysis.py  # Test phylogenetic analysis
│   ├── test_functional_prediction.py  # Test functional prediction
│   ├── test_visualization.py          # Test visualization
│   ├── integration_tests.py           # End-to-end integration tests
│   └── test_data/                     # Test datasets
├── workflows/
│   ├── snakemake/                     # Snakemake workflow
│   │   ├── Snakefile                  # Main Snakemake file
│   │   ├── config.yaml                # Snakemake configuration
│   │   └── rules/                     # Individual rule files
│   ├── nextflow/                      # Nextflow workflow
│   │   ├── main.nf                    # Main Nextflow file
│   │   ├── nextflow.config            # Nextflow configuration
│   │   └── modules/                   # Individual modules
│   └── cwl/                           # Common Workflow Language
│       ├── hbv-mtct-workflow.cwl      # Main CWL workflow
│       └── tools/                     # Individual tool definitions
├── docker/
│   ├── Dockerfile                     # Docker container definition
│   ├── docker-compose.yml            # Docker Compose configuration
│   ├── requirements.txt               # Python requirements
│   └── environment.yml               # Conda environment file
├── singularity/
│   ├── hbv-mtct-analysis.def          # Singularity definition file
│   └── build_container.sh            # Container build script
├── scripts/
│   ├── download_references.py         # Download reference genomes
│   ├── validate_inputs.py             # Input validation script
│   ├── merge_results.py               # Merge analysis results
│   ├── generate_report.py             # Generate final report
│   └── backup_results.sh              # Backup important results
├── examples/
│   ├── README.md                       # Examples documentation
│   ├── basic_analysis.sh              # Basic analysis example
│   ├── advanced_analysis.sh           # Advanced analysis example
│   ├── custom_configuration.yaml      # Custom config example
│   └── sample_outputs/                # Example output files
├── .github/
│   ├── workflows/                     # GitHub Actions
│   │   ├── ci.yml                     # Continuous integration
│   │   ├── test.yml                   # Automated testing
│   │   └── release.yml                # Release automation
│   ├── ISSUE_TEMPLATE/                # Issue templates
│   │   ├── bug_report.md              # Bug report template
│   │   ├── feature_request.md         # Feature request template
│   │   └── question.md                # Question template
│   ├── PULL_REQUEST_TEMPLATE.md       # PR template
│   └── CONTRIBUTING.md                # Contribution guidelines
├── CONTRIBUTING.md                    # How to contribute
├── CODE_OF_CONDUCT.md                 # Code of conduct
├── SECURITY.md                        # Security policy
└── MANIFEST.in                        # Python manifest file
```

## Key File Descriptions

### Core Scripts

| File | Purpose | Input | Output |
|------|---------|-------|--------|
| `run_pipeline.sh` | Main pipeline executor | FASTQ files | Complete analysis |
| `concatenate_fastq.sh` | Merge FASTQ by barcode | Raw FASTQ files | Concatenated FASTQ |
| `artic_genotype_c.sh` | ARTIC consensus for genotype C | FASTQ files | Consensus FASTA |
| `provean_analysis.sh` | Functional impact prediction | Protein sequences | PROVEAN scores |
| `epitope_analysis.py` | Epitope mapping analysis | Protein sequences | Epitope maps |

### Configuration Files

| File | Purpose | Format |
|------|---------|--------|
| `pipeline_config.yaml` | Main configuration | YAML |
| `hbv_gene_coordinates.yaml` | Gene coordinates | YAML |
| `sample_metadata.yaml` | Sample information | YAML |
| `cluster_config.yaml` | HPC settings | YAML |

### Documentation Files

| File | Purpose |
|------|---------|
| `README.md` | Main project documentation |
| `INSTALLATION.md` | Installation instructions |
| `USAGE.md` | Usage guide |
| `CITATIONS.md` | Software citations |
| `TROUBLESHOOTING.md` | Common issues |

## Directory Navigation Links

Each directory contains a `README.md` file with:

1. **Module overview** - What the module does
2. **Input requirements** - Required input files
3. **Output description** - Generated output files
4. **Usage examples** - How to run the module
5. **Parameter explanation** - Configurable parameters
6. **Dependencies** - Required software/packages

### Example Module README Structure

```markdown
# Module Name

## Overview
Brief description of what this module does.

## Input Files
- `input1.fastq` - Description
- `input2.fasta` - Description

## Output Files
- `output1.txt` - Description
- `output2.pdf` - Description

## Usage
```bash
./module_script.sh --input data/ --output results/
```

## Parameters
- `--threads` - Number of CPU threads (default: 8)
- `--memory` - Memory allocation (default: 32G)

## Dependencies
- Tool 1 (version)
- Tool 2 (version)
```

## GitHub Repository Features

### Badges
- DOI badge (when published)
- License badge
- Build status
- Version badge
- Download counter

### GitHub Actions
- **Continuous Integration**: Automated testing on push/PR
- **Testing**: Run test suite on multiple OS
- **Release**: Automated version releases
- **Documentation**: Auto-generate docs

### Issue Templates
- Bug reports with system info
- Feature requests with use cases
- Questions with guidelines

### Pull Request Template
- Checklist for contributions
- Testing requirements
- Documentation updates

## Repository Setup Commands

```bash
# Initialize repository
git init
git add .
git commit -m "Initial commit: HBV MTCT Analysis Pipeline v1.0.0"

# Add remote origin
git remote add origin https://github.com/yourusername/hbv-mtct-analysis.git

# Push to GitHub
git branch -M main
git push -u origin main

# Create development branch
git checkout -b develop
git push -u origin develop

# Tag first release
git tag -a v1.0.0 -m "Version 1.0.0: Initial release"
git push origin v1.0.0
```

## Repository Maintenance

### Regular Updates
- Keep dependencies updated
- Update documentation
- Add new features based on user feedback
- Fix bugs and issues

### Version Control
- Use semantic versioning (MAJOR.MINOR.PATCH)
- Tag releases appropriately
- Maintain CHANGELOG.md

### Community Engagement
- Respond to issues promptly
- Review pull requests
- Update documentation based on user feedback
- Provide support through discussions

This complete structure ensures a professional, maintainable, and user-friendly repository that follows GitHub best practices and provides comprehensive documentation for reproducible research.
