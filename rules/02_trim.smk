rule trim_galore_pe:
    input:
        [INPUT_DIR + "/{smp}_R1.fastq.gz", INPUT_DIR + "/{smp}_R2.fastq.gz"]
    output:
        "results/02_trim/{smp}_R1_val_1.fq.gz",
        "results/02_trim/{smp}_R1.fastq.gz_trimming_report.txt",
        "results/02_trim/{smp}_R2_val_2.fq.gz",
        "results/02_trim/{smp}_R2.fastq.gz_trimming_report.txt"
    params:
        extra="--illumina -q 20"
    log:
        "results/02_trim/logs/{smp}.log"
    wrapper:
        "0.35.2/bio/trim_galore/pe"

rule trim_galore_multiqc:
    input:
        expand("results/02_trim/{smp}_R1.fastq.gz_trimming_report.txt", smp=sample_id),
        expand("results/02_trim/{smp}_R2.fastq.gz_trimming_report.txt", smp=sample_id)
    output:
        "results/02_trim/trim_galore_multiqc_report.html"
    log:
        "results/02_trim/logs/multiqc.log"
    wrapper:
        "0.31.1/bio/multiqc"