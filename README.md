# SnakemakeRNASeq
A Snakemake workflow to process paired-end RNA-Seq data.

## Steps:

**The workflow consists following steps:**

- Quality control of the raw and/or trimmed data (FastQC, MultiQC)
- Adapter trimming w/ trim_galore (Optional)
- Contamination check and decontamination (Optional)
- Alignment to the reference genome (hisat2, STAR)
- Quality control with RSeQC, QualiMap
- Transcript/gene quantification (StringTie, featureCounts, RSEM)
- Alignment-free transcript quantification (kallisto/salmon)

**Future additions:**
- Differential gene expression analysis (deseq2)

![dag](https://user-images.githubusercontent.com/42179487/74106710-ec1f1a80-4b36-11ea-94b2-d52f4bee5574.png)





