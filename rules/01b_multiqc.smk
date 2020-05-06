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
