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

![dag](https://user-images.githubusercontent.com/42179487/60461360-53293800-9c14-11e9-8161-ff5a0973b786.png)



