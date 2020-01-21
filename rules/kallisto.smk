rule kallisto_index:
    input:
        tcp = TRANSCRIPTS
    output:
        "results/kallisto/index/transcripts.idx"
    log:
        "results/kallisto/logs/index.log"
    priority:50
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input} 2> {log}"

rule kallisto_quant:
    input:
        fastq = ["results/trimmed/{smp}_R1_val_1.fq.gz", "results/trimmed/{smp}_R2_val_2.fq.gz"],
        index = "results/kallisto/index/transcripts.idx"
    output:
        directory("results/kallisto/quant/quant_results_{smp}")
    params:
<<<<<<< HEAD
        extra = "-b 100 --bias --fusion --gtf work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf"
=======
        extra = "-b 100"
>>>>>>> 58c7e0000fcb1f754282d482021c355d4887289a
    log:
        "results/kallisto/logs/kallisto_quant_{smp}.log"
    threads: 16
    wrapper:
        "0.36.0/bio/kallisto/quant"