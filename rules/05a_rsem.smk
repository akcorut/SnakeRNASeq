rule rsem_genome:
    input:
        gff= config["ref"]["annotation"],
        fasta = REFERENCE
    output:
        rsemindex = config["rsem"]["rsemindex"] + ".n2g.idx.fa"
    threads: 24
    params:
        ref_name = lambda wildcards, output: output[0][:-11]
    conda:
        "../envs/rsem.yaml"
    priority:50
    shell:
        """
        rsem-prepare-reference \
        -p {threads} \
        --gff3 {input.gff} {input.fasta} {params.ref_name}
        """

rule rsem_calculate:
    input:
        bam= "results/04_alignment/04a_alignment_results/star/pass2/{smp}/Aligned.toTranscriptome.out.bam"
    output:
        genes = "results/05_quantification/05a_rsem/genes/{smp}.genes.results",
    params:
        genomedir= config["rsem"]["rsemindex"],
        prefix = lambda wildcards, output: output[0][:-14]
    threads: 24
    conda:
        "../envs/rsem.yaml"
    priority: -10
    shell:
        """
        rsem-calculate-expression \
        --no-qualities \
        -p {threads} \
        --strandedness reverse \
        --alignments --paired-end {input.bam} {params.genomedir} {params.prefix}
        """

rule multiqc_rsem:
    input:
       expand("results/05_quantification/05a_rsem/genes/{smp}.stat", smp=sample_id)
    output:
        "results/05_quantification/05a_rsem/rsem_multiqc.html"
    log:
        "results/05_quantification/05a_rsem/logs/multiqc.log"
    wrapper:
        "0.49.0/bio/multiqc"