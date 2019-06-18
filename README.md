# SnakemakeRNASeq
A Snakemake workflow to process paired-end RNA-Seq data.

**Steps:**

- Quality control check on raw data (FastQC)
- Summarise the quality control results (MultiQC)
- Adapter trimming with cutadapt
- Post quality control check after trimming
- Summarise the post trimming QC results 
- Reference genome indexing
- Mapping to reference genome using hisat2

![dag](https://user-images.githubusercontent.com/42179487/59709094-f66d5c80-91d3-11e9-873d-cc0e003f4585.png)



