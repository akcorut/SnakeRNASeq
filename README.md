# SnakemakeRNASeq
A Snakemake workflow to process paired-end RNA-Seq data.

**Steps:**

- Quality control check on raw data (FastQC)
- Summarise the quality control results (MultiQC)
- Adapter trimming with cutadapt / trim_galore
- Post quality control check after trimming
- Summarise the post trimming QC results 
- Reference genome indexing
- Mapping to reference genome using hisat2 and/or STAR
- RSeQC evaluation
- Assemble and quantify expressed genes and transcripts with Stringtie

![dag](https://user-images.githubusercontent.com/42179487/61956450-20f8c500-af8b-11e9-847b-faa70ad33772.png)



