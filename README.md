# SnakemakeRNASeq
A Snakemake workflow to process paired-end RNA-Seq data.

**Steps:**

- Quality control check on raw data (FastQC)
- Summarise the quality control results (MultiQC)
- Adapter trimming with cutadapt
- Post quality control check after trimming
- Summarise the post trimming QC results 
- De novo assembly with trinity (Protocol:[_Haas et al.,2013_](https://www.nature.com/articles/nprot.2013.084.pdf))

![dag](https://user-images.githubusercontent.com/42179487/59511884-6baff900-8e85-11e9-8865-599dae245db1.png)


