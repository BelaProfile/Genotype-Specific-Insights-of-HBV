# Software Citations and Acknowledgments

This pipeline uses numerous open-source tools and databases. Please cite the appropriate papers when using this pipeline.

## Core Bioinformatics Tools

### Sequencing and Assembly

**ARTIC Network Pipeline**
- **Citation**: Loman, N.J., et al. (2020). A user-friendly ARTIC pipeline for viral genome sequencing and real-time surveillance. *bioRxiv*. DOI: 10.1101/2020.04.17.046086
- **GitHub**: https://github.com/artic-network/artic-ncov2019
- **Version used**: 1.2.1
- **License**: MIT

**Medaka**
- **Citation**: Oxford Nanopore Technologies (2021). Medaka: Sequence correction provided by ONT Research
- **GitHub**: https://github.com/nanoporetech/medaka  
- **Version used**: 1.4.4
- **License**: Mozilla Public License 2.0

**Minimap2**
- **Citation**: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094-3100. DOI: 10.1093/bioinformatics/bty191
- **GitHub**: https://github.com/lh3/minimap2
- **Version used**: 2.24
- **License**: MIT

### Sequence Processing

**SAMtools**
- **Citation**: Li, H., et al. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078-2079. DOI: 10.1093/bioinformatics/btp352
- **Website**: http://www.htslib.org/
- **Version used**: 1.15
- **License**: MIT/Expat

**BCFtools**
- **Citation**: Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. DOI: 10.1093/gigascience/giab008
- **Website**: http://www.htslib.org/
- **Version used**: 1.15
- **License**: MIT/Expat

**BEDtools**
- **Citation**: Quinlan, A.R. & Hall, I.M. (2010). BEDtools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, 26(6), 841-842. DOI: 10.1093/bioinformatics/btq033
- **GitHub**: https://github.com/arq5x/bedtools2
- **Version used**: 2.30.0
- **License**: MIT

**EMBOSS**
- **Citation**: Rice, P., et al. (2000). EMBOSS: the European Molecular Biology Open Software Suite. *Trends in Genetics*, 16(6), 276-277. DOI: 10.1016/S0168-9525(00)02024-2
- **Website**: http://emboss.sourceforge.net/
- **Version used**: 6.6.0
- **License**: GPL

### Phylogenetic Analysis

**MAFFT**
- **Citation**: Katoh, K. & Standley, D.M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772-780. DOI: 10.1093/molbev/mst010
- **Website**: https://mafft.cbrc.jp/alignment/software/
- **Version used**: 7.490
- **License**: BSD

**FastTree**
- **Citation**: Price, M.N., et al. (2010). FastTree 2 – approximately maximum-likelihood trees for large alignments. *PLoS ONE*, 5(3), e9490. DOI: 10.1371/journal.pone.0009490
- **Website**: http://www.microbesonline.org/fasttree/
- **Version used**: 2.1.11
- **License**: GPL

### Variant Calling

**LoFreq**
- **Citation**: Wilm, A., et al. (2012). LoFreq: a sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. *Nucleic Acids Research*, 40(22), 11189-11201. DOI: 10.1093/nar/gks918
- **Website**: https://csb5.github.io/lofreq/
- **Version used**: 2.1.5
- **License**: MIT

### Functional Prediction

**PROVEAN**
- **Citation**: Choi, Y. & Chan, A.P. (2015). PROVEAN web server: a tool to predict the functional effect of amino acid substitutions and indels. *Bioinformatics*, 31(16), 2745-2747. DOI: 10.1093/bioinformatics/btv195
- **Website**: http://provean.jcvi.org/
- **Version used**: 1.1.5
- **License**: Academic License

**SIFT4G**
- **Citation**: Vaser, R., et al. (2016). SIFT missense predictions for genomes. *Nature Protocols*, 11(1), 1-9. DOI: 10.1038/nprot.2015.123
- **GitHub**: https://github.com/rvaser/sift4g
- **Version used**: 2.0.0
- **License**: MIT

### Epitope Prediction

**BepiPred 2.0**
- **Citation**: Jespersen, M.C., et al. (2017). BepiPred-2.0: improving sequence-based B-cell epitope prediction using conformational epitopes. *Nucleic Acids Research*, 45(W1), W24-W29. DOI: 10.1093/nar/gkx346
- **Website**: https://services.healthtech.dtu.dk/service.php?BepiPred-2.0
- **Version used**: 2.0
- **License**: Academic License

### Structural Prediction

**AlphaFold**
- **Citation**: Jumper, J., et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596(7873), 583-589. DOI: 10.1038/s41586-021-03819-2
- **Website**: https://alphafold.ebi.ac.uk/
- **Version used**: 2.0
- **License**: Apache 2.0

### Sequence Analysis

**SeqKit**
- **Citation**: Shen, W., et al. (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLoS ONE*, 11(10), e0163962. DOI: 10.1371/journal.pone.0163962
- **GitHub**: https://github.com/shenwei356/seqkit
- **Version used**: 2.3.0
- **License**: MIT

**BLAST+**
- **Citation**: Camacho, C., et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. DOI: 10.1186/1471-2105-10-421
- **Website**: https://blast.ncbi.nlm.nih.gov/
- **Version used**: 2.12.0
- **License**: Public Domain

## Programming Languages and Libraries

### Python (3.8+)
- **Citation**: Van Rossum, G. & Drake, F.L. (2009). Python 3 Reference Manual. Scotts Valley, CA: CreateSpace.
- **Website**: https://www.python.org/
- **License**: Python Software Foundation License

**Python Libraries:**

**Pandas**
- **Citation**: McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 51-56.
- **Version used**: Latest stable
- **License**: BSD

**NumPy**
- **Citation**: Harris, C.R., et al. (2020). Array programming with NumPy. *Nature*, 585(7825), 357-362. DOI: 10.1038/s41586-020-2649-2
- **Version used**: Latest stable
- **License**: BSD

**SciPy**
- **Citation**: Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17(3), 261-272. DOI: 10.1038/s41592-019-0686-2
- **Version used**: Latest stable
- **License**: BSD

**Matplotlib**
- **Citation**: Hunter, J.D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90-95. DOI: 10.1109/MCSE.2007.55
- **Version used**: Latest stable
- **License**: PSF-based

**Seaborn**
- **Citation**: Waskom, M.L. (2021). seaborn: statistical data visualization. *Journal of Open Source Software*, 6(60), 3021. DOI: 10.21105/joss.03021
- **Version used**: Latest stable
- **License**: BSD

**Biopython**
- **Citation**: Cock, P.J., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423. DOI: 10.1093/bioinformatics/btp163
- **Version used**: Latest stable
- **License**: Biopython License Agreement

### R (4.2+)
- **Citation**: R Core Team (2023). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.
- **Website**: https://www.R-project.org/
- **License**: GPL-2 | GPL-3

**R Libraries:**

**ggplot2**
- **Citation**: Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.
- **Version used**: Latest stable
- **License**: MIT

**dplyr**
- **Citation**: Wickham, H., et al. (2023). dplyr: A Grammar of Data Manipulation.
- **Version used**: Latest stable
- **License**: MIT

**tidyr**
- **Citation**: Wickham, H., et al. (2023). tidyr: Tidy Messy Data.
- **Version used**: Latest stable
- **License**: MIT

**viridis**
- **Citation**: Garnier, S. (2021). viridis: Colorblind-Friendly Color Maps for R.
- **Version used**: Latest stable
- **License**: MIT

**pheatmap**
- **Citation**: Kolde, R. (2019). pheatmap: Pretty Heatmaps.
- **Version used**: Latest stable
- **License**: GPL-2

**Biostrings (Bioconductor)**
- **Citation**: Pagès, H., et al. (2023). Biostrings: Efficient manipulation of biological strings.
- **Version used**: Latest stable
- **License**: Artistic-2.0

**vegan**
- **Citation**: Oksanen, J., et al. (2022). vegan: Community Ecology Package.
- **Version used**: Latest stable
- **License**: GPL-2

## Databases

**Swiss-Prot (UniProt)**
- **Citation**: The UniProt Consortium (2023). UniProt: the universal protein knowledgebase in 2023. *Nucleic Acids Research*, 51(D1), D523-D531. DOI: 10.1093/nar/gkac1052
- **Website**: https://www.uniprot.org/
- **License**: CC BY 4.0

**IEDB (Immune Epitope Database)**
- **Citation**: Vita, R., et al. (2019). The Immune Epitope Database (IEDB): 2018 update. *Nucleic Acids Research*, 47(D1), D339-D343. DOI: 10.1093/nar/gky1006
- **Website**: https://www.iedb.org/
- **License**: Public Domain

**NCBI Reference Sequences**
- **Citation**: O'Leary, N.A., et al. (2016). Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Research*, 44(D1), D733-D745. DOI: 10.1093/nar/gkv1189
- **Website**: https://www.ncbi.nlm.nih.gov/refseq/
- **License**: Public Domain

## Infrastructure and Environment Management

**Conda**
- **Citation**: Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.
- **Website**: https://docs.conda.io/
- **License**: BSD

**SLURM Workload Manager**
- **Citation**: Yoo, A.B., et al. (2003). SLURM: Simple Linux Utility for Resource Management. *Job Scheduling Strategies for Parallel Processing*, 44-60. DOI: 10.1007/10968987_3
- **Website**: https://slurm.schedmd.com/
- **License**: GPL

## Data Visualization Tools

**iTOL (Interactive Tree of Life)**
- **Citation**: Letunic, I. & Bork, P. (2021). Interactive Tree Of Life (iTOL) v5: an online tool for phylogenetic tree display and annotation. *Nucleic Acids Research*, 49(W1), W293-W296. DOI: 10.1093/nar/gkab301
- **Website**: https://itol.embl.de/
- **License**: Academic License

## Acknowledgments

We acknowledge the developers and maintainers of all the above tools and databases. This work would not be possible without their contributions to the open-source bioinformatics community.

## Version Information

For exact version information used in this analysis, run:
```bash
./setup/show_versions.sh
```

## License Compliance

This pipeline is distributed under the MIT License. All included tools and databases are used in accordance with their respective licenses. For commercial use, please ensure compliance with all component licenses, particularly:

- PROVEAN (Academic License)
- BepiPred 2.0 (Academic License)  
- iTOL (Academic License)

## How to Cite This Pipeline

If you use this pipeline in your research, please cite this repository along with the relevant tools used in your specific analysis:

```bibtex
@software{hbv_mtct_pipeline,
  author = {Your Name},
  title = {HBV Mother-to-Child Transmission Analysis Pipeline},
  url = {https://github.com/yourusername/hbv-mtct-analysis},
  version = {1.0.0},
  year = {2024}
}
```

---

*Last updated: December 2024*
