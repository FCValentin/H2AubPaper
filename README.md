# Sperm derived H2AK119ub1 is required for embryonic development in Xenopus Laevis

__Valentin Francois--Campion, Florian Berger, Mami Oikawa, Maissa Goumeidane, Romain Guibeaux, and Jérôme Jullien__

H3K4me3, H3K27me3 and H2AK119ub1 have been analyzed during this protocol with ChIP-Sequencing. Transcriptomic of Xenopus laevis embryo development was analysed by RNA-Sequencing

## Getting Started

This repository reports all scripts used for bio-informatics analysis. Data are deposited on ENA (PRJEB56442). 

### Prerequisites and genome version
Genome version and annotation

* Xenopus laevis 9.2 from [Xenbase](https://download.xenbase.org/xenbase/Genomics/JGI/Xenla9.2/)
* 601DNA sequence (100% ubiquitylated semi-synthetic nucleosome) : tagctcaattggtcgtagcaagctctagcaccgcttaaacgcacgtacgcgctgtcccccgcgttttaaccgccaaggggattactccctagtctccagg
* CpG islands (cpgIslandExt.txt.gz) and repeats elements (rmsk.txt.gz) annotation file from [UCSC database](https://hgdownload.soe.ucsc.edu/goldenPath/xenLae2/database/)

### Tools used for analysis (conda environment)

* ChIP-Seq conda environment => [Chip.yml](https://github.com/FCValentin/H2AubPaper/blob/main/Chip.yml)
* RNA-Seq conda environment => [rnaseq_align.yml](https://github.com/FCValentin/H2AubPaper/blob/main/rnaseq_align.yml)
* R libraries used are shown in [ProcessingChIPdata.R](https://github.com/FCValentin/H2AubPaper/blob/main/ProcessingChIPdata.R)

### Scripts used for analysis

* RNA-Seq alignment [snakemake](https://gitlab.univ-nantes.fr/E114424Z/rnaseq_align)
* RNA-Seq R scripts for [gene expression analysis](https://gitlab.univ-nantes.fr/E114424Z/BulkRNAseq)
* All R scripts for permutation, TF motif heatmap, P-Value of particles enrichment and genomic caracterisation => [ProcessingChIPdata.R](https://github.com/FCValentin/H2AubPaper/blob/main/ProcessingChIPdata.R) 
* All ChIP-Seq command line for data processing => [ChIP_Sequencing.sh](https://github.com/FCValentin/H2AubPaper/blob/main/ChIP_Sequencing.sh)
* Scripts for H2A (particles) enrichment detection => [H2A_enrichment.sh](https://github.com/FCValentin/H2AubPaper/blob/main/H2A_enrichment.sh)

## Authors

**Valentin FRANCOIS--CAMPION** 
