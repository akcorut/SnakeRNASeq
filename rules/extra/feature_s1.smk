rule featureCounts_s1:
    input:
        bam= expand("results/star/{smp}/Aligned.out.bam", smp=sample_id),
        gtf= "/work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf",
        fasta = REFERENCE
    output:
        "results/feature_s1/featureCounts_s1_results.txt"
    conda:
        "../envs/subread.yaml"
    threads:20
    shell:
        """
        featureCounts -T {threads} -s 1 -p -t exon -g gene_id -M --fraction -J -G {input.fasta} -a {input.gtf} -o {output} {input.bam} 
        """