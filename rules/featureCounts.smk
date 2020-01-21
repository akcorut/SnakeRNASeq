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
<<<<<<< HEAD
        featureCounts -T {threads} -p -t exon -g gene_id -M --fraction -J -G {input.fasta} -a {input.gtf} -o {output} {input.bam} 
=======
        featureCounts -T {threads} -p -t exon -g gene_id -M --fraction -J -G {input.fasta} -a {input.gtf} -o {output} {input.bam}
>>>>>>> 58c7e0000fcb1f754282d482021c355d4887289a
        """