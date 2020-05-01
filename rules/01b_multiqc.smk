rule multiqc:
    input:
        expand("results/01_qc/01a_fqc/fqc_init/{smp}_R1_fastqc.html", smp=sample_id),
        expand("results/01_qc/01a_fqc/fqc_init/{smp}_R2_fastqc.html", smp=sample_id) 
    output:
        "results/01b_multiqc/multiqc_init/fastq_multiqc.html"
    log:
        "results/01b_multiqc/multiqc_init/logs/multiqc.log"
    wrapper:
        "0.35.0/bio/multiqc"

rule multiqc_trim:
    input:
        expand("results/01_qc/01a_fqc/fqc_trim/{smp}_R1_val_1_fastqc.html", smp=sample_id),
        expand("results/01_qc/01a_fqc/fqc_trim/{smp}_R2_val_2_fastqc.html", smp=sample_id)
    output:
        "results/01b_multiqc/multiqc_trim/multiqc_report_trim_galore.html"
    log:
        "results/01b_multiqc/multiqc_trim/logs/multiqc.log"
    wrapper:
        "0.31.1/bio/multiqc"

rule multiqc_clean:
    input:
        expand("results/01_qc/01a_fqc/fqc_clean/{smp}_R1_clean_fastqc.html", smp=sample_id),
        expand("results/01_qc/01a_fqc/fqc_clean/{smp}_R2_clean_fastqc.html", smp=sample_id)
    output:
        "results/01_qc/01b_multiqc/multiqc_clean/multiqc_report_decont.html"
    log:
        "results/01_qc/01b_multiqc/multiqc_clean/logs/multiqc.log"
    wrapper:
        "0.47.0/bio/multiqc"
