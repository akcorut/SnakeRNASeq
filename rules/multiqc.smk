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

rule multiqc_trim_galore:
    input:
        expand("results/FastQCTrim/{smp}_R1_val_1_fastqc.html", smp=sample_id),
        expand("results/FastQCTrim/{smp}_R2_val_2_fastqc.html", smp=sample_id)
    output:
        "results/MultiQCTrim/multiqc_report_trim_galore.html"
    log:
        "results/MultiQCTrim/logs/multiqc.log"
    wrapper:
        "0.31.1/bio/multiqc"

rule multiqc_feature:
    input:
        "results/feature/featureCounts_results.txt.summary"
    output:
        "results/feature/feature_multiqc.html"
    log:
        "results/feature/logs/multiqc.log"
    wrapper:
        "0.35.0/bio/multiqc"

rule multiqc_kallisto:
    input:
        expand("results/kallisto/logs/kallisto_quant_{smp}.log", smp=sample_id),
    output:
        "results/kallisto/kallisto_multiqc.html"
    log:
        "results/kallisto/logs/multiqc.log"
    wrapper:
        "0.35.0/bio/multiqc"