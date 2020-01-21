# SnakemakeRNASeq
A Snakemake workflow to process paired-end RNA-Seq data.

**Steps:**

The workflow consists following steps:

- Quality control of the raw and trimmed data (FastQC, MultiQC)
- Adapter trimming w/ trim_galore (Optional)
- Alignment/mapping to the reference genome (hisat2, STAR)
- Quality control with RSeQC
- Transcript quantification (StringTie, featureCounts)
- Alignment-free transcript quantification (kallisto/salmon)

Potential future additions:
- Machine learning based gene expression level quality control (GeneQC)
- Differential gene expression analysis (deseq2/ballgown/sleuth) 

![dag](https://user-images.githubusercontent.com/42179487/65396555-4cbede00-dd75-11e9-8eaa-4c58991f3dcc.png)



