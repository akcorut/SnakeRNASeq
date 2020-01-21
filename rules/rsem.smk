rule rsem_genome:
    input:
        gtf= "/work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf",
        fasta = REFERENCE
    output:
        rsemindex = config["rsem"]["rsemindex"] + ".n2g.idx.fa"
    threads: 8
    params:
        ref_name = config["rsem"]["rsemindex"]
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-prepare-reference \
        -p {threads} \
        --gtf {input.gtf} {input.fasta} {params.ref_name}
        """

rule rsem:
    input:
        bam= "results/star/{smp}/Aligned.toTranscriptome.out.bam"
    output:
        genes = "results/rsem/genes/{smp}.genes.results"
    params:
        genomedir = config["rsem"]["rsemindex"],
        prefix = config["rsem"]["genes"] + "/" + "{smp}"
    threads: 16
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-calculate-expression \
        --paired-end \
        --no-bam-output \
        --quiet \
        --no-qualities \
        -p {threads} \
        --forward-prob 0.5 \
        --seed-length 25 \
        --fragment-length-mean -1.0 \
        --bam {input.bam} {params.genomedir} {params.prefix}
        """