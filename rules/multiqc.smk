rule multiqc:
    input:
        expand("results/FastQC/{smp}_R1_fastqc.html", smp=sample_id),
        expand("results/FastQC/{smp}_R2_fastqc.html", smp=sample_id) 
    output:
        "results/MultiQC/fastq_multiqc.html"
    log:
        "results/MultiQC/logs/multiqc.log"
    wrapper:
        "0.35.0/bio/multiqc"

rule multiqc_after:
    input:
        expand("results/FastQCCut/{smp}_cutadapt_R1_fastqc.html", smp=sample_id),
        expand("results/FastQCCut/{smp}_cutadapt_R2_fastqc.html", smp=sample_id)
    output:
        "results/MultiQCCut/fastq_cutadapt_multiqc.html"
    log:
        "results/MultiQCCut/logs/multiqc_cutadapt.log"
    wrapper:
        "0.35.0/bio/multiqc"