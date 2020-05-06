rule fastqc:
    input:
        GetFastq
    output:
        html_fwd="results/01_qc/01a_fqc/fqc_init/{smp}_R1_fastqc.html",
        zip_fwd="results/01_qc/01a_fqc/fqc_init/{smp}_R1_fastqc.zip",
        html_rev="results/01_qc/01a_fqc/fqc_init/{smp}_R2_fastqc.html",
        zip_rev="results/01_qc/01a_fqc/fqc_init/{smp}_R2_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    threads:20
    log:
        log_fwd="results/01_qc/01a_fqc/fqc_init/logs/{smp}_R1.fastqc.log",
        log_rev="results/01_qc/01a_fqc/fqc_init/logs/{smp}_R2.fastqc.log"
    shell:
        """
        fastqc --outdir /results/01_qc/01a_fqc/fqc_init -t {threads} -f fastq {input}
        """
