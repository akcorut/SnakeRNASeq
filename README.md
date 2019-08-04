# SnakemakeRNASeq
A Snakemake workflow to process paired-end RNA-Seq data.

**Steps:**
The workflow consists following steps/rules with potential additions and extractions:

- Quality control of the raw data (FastQC, MultiQC)
- Adapter trimming (cutadapt, trim_galore)
- Alignment/mapping to the reference genome (hisat2, STAR)
- Quality control with RSeQC
- Transcript quantification

Potential future additions:

- De-novo assembly (trinity)
- Alignment-free transcript quantification (kallisto/salmon)
- Machine learning based gene expression level quality control
- Differential gene expression analysis (deseq2/ballgown/sleuth) 

![dag](https://user-images.githubusercontent.com/42179487/61956450-20f8c500-af8b-11e9-847b-faa70ad33772.png)



