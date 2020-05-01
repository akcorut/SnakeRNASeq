rule fastqc:
    input:
        fwd= INPUT_DIR + "/{smp}_R1.fastq.gz",
        rev= INPUT_DIR + "/{smp}_R2.fastq.gz"
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
        fastqc --outdir /results/01_qc/01a_fqc/fqc_init -t {threads} -f fastq {input.fwd} {input.rev}
        """

######################################### After Trimming ###########################################

rule fastqc_trim:
    input:
        r1="results/02_trim/{smp}_R1_val_1.fq.gz",
        r2="results/02_trim/{smp}_R2_val_2.fq.gz"
    output:
        "results/01_qc/01a_fqc/fqc_trim/{smp}_R1_val_1_fastqc.html",
        "results/01_qc/01a_fqc/fqc_trim/{smp}_R1_val_1_fastqc.zip",
        "results/01_qc/01a_fqc/fqc_trim/{smp}_R2_val_2_fastqc.html",
        "results/01_qc/01a_fqc/fqc_trim/{smp}_R2_val_2_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    threads:20
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} -q -f fastq -o results/01_qc/01a_fqc/fqc_trim/
        """

######################################### After Decontamination ####################################

rule fastqc_clean:
    input:
        r1="results/decontamination/03b_rrna_cleaned/{smp}_R1_clean.fq",
        r2="results/decontamination/03b_rrna_cleaned/{smp}_R2_clean.fq",
    output:
        "results/01_qc/01a_fqc/fqc_clean/{smp}_R1_clean_fastqc.html",
        "results/01_qc/01a_fqc/fqc_clean/{smp}_R1_clean_fastqc.zip",
        "results/01_qc/01a_fqc/fqc_clean/{smp}_R2_clean_fastqc.html",
        "results/01_qc/01a_fqc/fqc_clean/{smp}_R2_clean_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    threads:20
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} -q -f fastq -o results/01_qc/01a_fqc/fqc_clean/
        """