# SnakemakeRNASeq
A Snakemake workflow to process paired-end RNA-Seq data.

**Steps:**

The workflow consists following steps:

- Quality control of the raw and trimmed data (FastQC, MultiQC)
- Adapter trimming w/ trim_galore (Optional)
- Alignment/mapping to the reference genome (hisat2, STAR)
- Quality control with RSeQC, QualiMap
- Checking against rRNA contamination
- Transcript quantification (StringTie, featureCounts. RSEM)
- Alignment-free transcript quantification (kallisto/salmon)

Future additions:
- Differential gene expression analysis (deseq2)

![dag](https://user-images.githubusercontent.com/42179487/72843534-b938e380-3c68-11ea-8319-5aa78f392fba.png)





