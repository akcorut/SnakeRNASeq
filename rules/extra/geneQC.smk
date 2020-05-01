rule salmon_mappings_to_bam:
    input:
        sam= "results/salmon/mappings/{smp}_salmon_mappings"
    output:
        bam= temp("results/salmon/mappings/{smp}_salmon.bam")
    conda:
        "../envs/samtools.yaml"
    threads: 24
    shell:
        """
        samtools view -@ {threads} -bS {input.sam} > {output.bam}
        """

rule GeneQC:
    input:
        bam="results/salmon/mappings/{smp}_salmon.bam",
        gff= ANNOTATION,
        ref= REFERENCE
    output:
        "results/GeneQC/{smp}_salmon_out.txt"
    params: "results/GeneQC/{smp}"
    conda:
        "../envs/geneqc.yaml"
    priority:-1
    shell:
        "python scripts/GeneQC_Python/GeneQC2.py 1 {input.ref} {input.gff} {input.bam} {params}"
        
