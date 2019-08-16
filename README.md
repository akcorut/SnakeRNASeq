# SnakemakeRNASeq
A Snakemake workflow to process paired-end RNA-Seq data.

**Steps:**

The workflow consists following steps:

- Quality control of the raw and trimmed data (FastQC, MultiQC)
- Adapter trimming w/ trim_galore
- Alignment/mapping to the reference genome (hisat2, STAR)
- Quality control with RSeQC

Potential future additions:

- Transcript quantification (stringtie)
- De-novo assembly (trinity, rnaSPAdes)
- Alignment-free transcript quantification (kallisto/salmon)
- Machine learning based gene expression level quality control (GeneQC)
- Differential gene expression analysis (deseq2/ballgown/sleuth) 

![dag](https://user-images.githubusercontent.com/42179487/63191116-c0aeed80-c035-11e9-8392-45c40d385c74.png)



