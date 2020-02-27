rule rsem_genome:
    input:
        gff= config["ref"]["annotation"],
        fasta = REFERENCE
    output:
        rsemindex = config["rsem"]["rsemindex"] + ".n2g.idx.fa"
    threads: 18
    params:
        ref_name = config["rsem"]["rsemindex"],
        pattern = "mRNA"
    conda:
        "../envs/rsem.yaml"
    priority:50
    shell:
        """
        rsem-prepare-reference \
        -p {threads} \
        --gff3 {input.gff} --gff3-RNA-patterns {params.pattern} {input.fasta} {params.ref_name}
        """

rule rsem_meristem:
    input:
        bam_meristem= "results/star/{sample_meristem}/Aligned.toTranscriptome.out.bam"
    output:
        genes_meristem = "results/rsem/genes/{sample_meristem}.genes.results"
    params:
        genomedir_meristem = config["rsem"]["rsemindex"],
        prefix_meristem = config["rsem"]["genes"] + "/" + "{sample_meristem}"
    threads: 18
    conda:
        "../envs/rsem.yaml"
    priority: -1
    shell:
        """
        rsem-calculate-expression \
        --no-bam-output \
        --no-qualities \
        -p {threads} \
        --strandedness forward \
        --alignments --paired-end {input.bam_meristem} {params.genomedir_meristem} {params.prefix_meristem}
        """

rule rsem_leaf:
    input:
        bam_leaf= "results/star/{sample_leaf}/Aligned.toTranscriptome.out.bam"
    output:
        genes_leaf = "results/rsem/genes/{sample_leaf}.genes.results"
    params:
        genomedir_leaf = config["rsem"]["rsemindex"],
        prefix_leaf = config["rsem"]["genes"] + "/" + "{sample_leaf}"
    threads: 18
    conda:
        "../envs/rsem.yaml"
    priority: -10
    shell:
        """
        rsem-calculate-expression \
        --no-bam-output \
        --no-qualities \
        -p {threads} \
        --strandedness reverse \
        --alignments --paired-end {input.bam_leaf} {params.genomedir_leaf} {params.prefix_leaf}
        """

ruleorder: rsem_leaf > rsem_meristem