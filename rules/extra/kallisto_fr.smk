rule kallisto_quant_fr:
    input:
        fastq = ["results/trimmed/{smp}_R1_val_1.fq.gz", "results/trimmed/{smp}_R2_val_2.fq.gz"],
        index = "results/kallisto/index/transcripts.idx"
    output:
        directory("results/kallisto/quant_fr/{smp}_fr")
    params:
        extra = "-b 100 --fr-stranded --bias --gtf work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf"
    log:
        "results/kallisto/quant_fr/logs/kallisto_quant_{smp}.log"
    threads: 20
    wrapper:
        "0.47.0/bio/kallisto/quant"