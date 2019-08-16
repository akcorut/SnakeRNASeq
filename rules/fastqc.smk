rule fastqc:
    input:
        fwd= raw_data + "/{smp}_R1.fastq.gz",
        rev= raw_data + "/{smp}_R2.fastq.gz"
    output:
        html_fwd="results/FastQC/{smp}_R1_fastqc.html",
        zip_fwd="results/FastQC/{smp}_R1_fastqc.zip",
        html_rev="results/FastQC/{smp}_R2_fastqc.html",
        zip_rev="results/FastQC/{smp}_R2_fastqc.zip"
    threads:20
    log:
        log_fwd="results/FastQC/logs/{smp}_R1.fastqc.log",
        log_rev="results/FastQC/logs/{smp}_R2.fastqc.log"
    shell:
        """
        fastqc --outdir /results/FastQC -t {threads} -f fastq {input.fwd} {input.rev}
        """

rule fastqc_trim_galore:
    input:
        r1="results/trimmed/{smp}_R1_val_1.fq.gz",
        r2="results/trimmed/{smp}_R2_val_2.fq.gz"
    output:
        "results/FastQCTrim/{smp}_R1_val_1_fastqc.html",
        "results/FastQCTrim/{smp}_R1_val_1_fastqc.zip",
        "results/FastQCTrim/{smp}_R2_val_2_fastqc.html",
        "results/FastQCTrim/{smp}_R2_val_2_fastqc.zip"
    threads:40
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} -q -f fastq -o results/FastQCTrim/
        """