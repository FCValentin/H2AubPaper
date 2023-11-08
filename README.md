## Sperm derived H2AK119ub1 is required for embryonic development in Xenopus Laevis

__Valentin Francois--Campion, Florian Berger, Mami Oikawa, Maissa Goumeidane, Romain Guibeaux, and Jérôme Jullien__

H3K4me3, H3K27me3 and H2AK119Ub1 have been analyzed during this protocol with ChIP-Sequencing. Transcriptomic of Xenopus laevis embryo development was analysed by RNA-Sequencing

## Getting Started

This repository reports all scripts used for bio-informatics analysis. Data are deposited on ENA (GSEXXXXX). 

### Prerequisites and genome version
Genome version and annotation

* Xenopus Laevis 9.2 (from Xenbase)
* 601DNA sequence (100% ubiquitylated semi-synthetic nucleosome) : tagctcaattggtcgtagcaagctctagcaccgcttaaacgcacgtacgcgctgtcccccgcgttttaaccgccaaggggattactccctagtctccagg
* CpG islands (cpgIslandExt.txt.gz) and repeats elements (rmsk.txt.gz) annotation file from UCSC database https://hgdownload.soe.ucsc.edu/goldenPath/xenLae2/database/)

### Tools used for analysis (conda environment)

* Chip.yml => ChIP-Seq conda environment

* rnaseq_align.yml => RNA-Seq conda environment

### Scripts used for analysis

* https://gitlab.univ-nantes.fr/E114424Z/rnaseq_align => RNA-Seq alignment snakemake

* https://gitlab.univ-nantes.fr/E114424Z/BulkRNAseq => RNA-Seq R scripts for transcriptomics analysis

* ProcessingChIPdata.R => All R scripts for permutation, TF motif heatmap, P-Value of particles enrichment and genomic caracterisation

* ChIP_Sequencing.sh => All ChIP-Seq command line for data processing

* H2A_enrichment.sh => Scripts for H2A (particles) enrichment detection

## Authors

**Valentin FRANCOIS--CAMPION** 
