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
        "0.39.0/bio/multiqc"

rule multiqc_kallisto:
    input:
        expand("results/kallisto/logs/kallisto_quant_{smp}.log", smp=sample_id),
    output:
        "results/kallisto/kallisto_multiqc.html"
    log:
        "results/kallisto/logs/multiqc.log"
    wrapper:
        "0.35.0/bio/multiqc"

rule multiqc_salmon:
    input:
        expand("results/salmon/quant/{smp}_salmon_quant", smp=sample_id)
    output:
        "results/salmon/salmon_multiqc.html"
    log:
        "results/salmon/logs/multiqc.log"
    wrapper:
        "0.35.0/bio/multiqc"

rule multiqc_rsem:
    input:
        expand("results/rsem/genes/{smp}.stat", smp=sample_id)
    output:
        "results/rsem/rsem_multiqc.html"
    log:
        "results/rsem/logs/multiqc.log"
    wrapper:
        "0.47.0/bio/multiqc"

rule multiqc_qualimap:
    input:
        expand("results/qualimap/{smp}", smp=sample_id)
    output:
        "results/qualimap/qualimap_multiqc.html"
    log:
        "results/qualimap/logs/multiqc.log"
    wrapper:
        "0.47.0/bio/multiqc"