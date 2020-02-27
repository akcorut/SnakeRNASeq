rule featureCounts:
    input:
        bam= expand("results/star/{smp}/Aligned.out.bam", smp=sample_id),
        gtf= "/work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf",
        fasta = REFERENCE
    output:
        "results/feature/featureCounts_results.txt"
    conda:
        "../envs/subread.yaml"
    threads:20
    shell:
        """
        featureCounts -T {threads} -p -t exon -g gene_id -M --fraction -J -G {input.fasta} -a {input.gtf} -o {output} {input.bam} 
        """