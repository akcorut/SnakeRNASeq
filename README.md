# SnakemakeRNASeq
A Snakemake workflow to process paired-end RNA-Seq data.

**Steps:**

- Quality control check on raw data (FastQC)
- Summarise the quality control results (MultiQC)
- Adapter trimming with cutadapt
- Post quality control check after trimming
- Summarise the post trimming QC results 
- Reference genome indexing
- Mapping to reference genome using hisat2 and/or STAR
- RSeQC evaluation
- Assemble and quantify expressed genes and transcripts with Stringtie

![dag](https://user-images.githubusercontent.com/42179487/60604416-a2db4100-9d85-11e9-9c35-aa84c0da85fa.png)



