rule trim_galore_pe:
    input:
        GetFastq
    output:
        trim1="results/02_trim/{smp}_R1_val_1.fq.gz",
        report1="results/02_trim/{smp}_R1.fastq.gz_trimming_report.txt",
        trim2="results/02_trim/{smp}_R2_val_2.fq.gz",
        report2="results/02_trim/{smp}_R2.fastq.gz_trimming_report.txt"
    params:
        extra="--illumina -q 20"
    log:
        "results/02_trim/logs/{smp}.log"
    priority:4
    wrapper:
        "0.35.2/bio/trim_galore/pe"

rule fastqc_trim:
    input:
        r1=rules.trim_galore_pe.output.trim1,
        r2=rules.trim_galore_pe.output.trim2
    output:
        html1="results/01_qc/01a_fqc/fqc_trim/{smp}_R1_val_1_fastqc.html",
        zip1="results/01_qc/01a_fqc/fqc_trim/{smp}_R1_val_1_fastqc.zip",
        html2="results/01_qc/01a_fqc/fqc_trim/{smp}_R2_val_2_fastqc.html",
        zip2="results/01_qc/01a_fqc/fqc_trim/{smp}_R2_val_2_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    priority:3
    threads:20
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} -q -f fastq -o results/01_qc/01a_fqc/fqc_trim/
        """

rule multiqc_trim:
    input:
        expand("results/01_qc/01a_fqc/fqc_trim/{smp}_R1_val_1_fastqc.html", smp=sample_id),
        expand("results/01_qc/01a_fqc/fqc_trim/{smp}_R2_val_2_fastqc.html", smp=sample_id)
    output:
        "results/01b_multiqc/multiqc_trim/multiqc_report_trim_galore.html"
    log:
        "results/01b_multiqc/multiqc_trim/logs/multiqc.log"
    priority:2
    wrapper:
        "0.31.1/bio/multiqc"

rule trim_galore_multiqc:
    input:
        expand("results/02_trim/{smp}_R1.fastq.gz_trimming_report.txt", smp=sample_id),
        expand("results/02_trim/{smp}_R2.fastq.gz_trimming_report.txt", smp=sample_id)
    output:
        "results/02_trim/trim_galore_multiqc_report.html"
    log:
        "results/02_trim/logs/multiqc.log"
    priority:1
    wrapper:
        "0.31.1/bio/multiqc"