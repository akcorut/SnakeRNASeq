rule trim_galore_pe:
    input:
        [raw_data + "/{smp}_R1.fastq.gz", raw_data + "/{smp}_R2.fastq.gz"]
    output:
        "/scratch/ac32082/PeanutRnaSeq/results/trimmed/{smp}_R1_val_1.fq.gz",
        "/scratch/ac32082/PeanutRnaSeq/results/trimmed/{smp}_R1.fastq.gz_trimming_report.txt",
        "/scratch/ac32082/PeanutRnaSeq/results/trimmed/{smp}_R2_val_2.fq.gz",
        "/scratch/ac32082/PeanutRnaSeq/results/trimmed/{smp}_R2.fastq.gz_trimming_report.txt"
    params:
        extra="--illumina -q 20"
    log:
        "/scratch/ac32082/PeanutRnaSeq/results/trimmed/logs/{smp}.log"
    wrapper:
        "0.35.2/bio/trim_galore/pe"

rule trim_galore_multiqc:
    input:
        expand("/scratch/ac32082/PeanutRnaSeq/results/trimmed/{smp}_R1.fastq.gz_trimming_report.txt", smp=sample_id),
        expand("/scratch/ac32082/PeanutRnaSeq/results/trimmed/{smp}_R2.fastq.gz_trimming_report.txt", smp=sample_id)
    output:
        "results/trimmed/trim_galore_multiqc_report.html"
    log:
        "results/trimmed/logs/multiqc.log"
    wrapper:
        "0.31.1/bio/multiqc"