rule sort_star_bam:
    input:
        bam= "results/04_alignment/04a_alignment_results/star/pass2/{smp}/Aligned.sortedByCoord.out.bam"
    output:
        temp("results/feature/sorted_bam/{smp}_sorted.bam")
    params:
        prfx= lambda wildcards, output: output[0][:-4]
    priority:50
    conda:
        "../envs/samtools.yaml"
    priority:50
    threads:20
    shell:
        '''
        samtools sort -@ {threads} -n -o {output} -T {params.prfx} {input.bam}
        '''

rule featureCounts:
    input:
        bam= expand("results/feature/sorted_bam/{smp}_sorted.bam", smp=sample_id),
        gtf= rules.gff3_to_gtf.output.gtf,
        fasta = REFERENCE
    output:
        "results/feature/featureCounts_results.txt"
    log:
        "results/feature/featureCount.log"
    conda:
        "../envs/subread.yaml"
    threads:24
    shell:
        """
        featureCounts -T {threads} -p -s 2 -t exon -g gene_id -J -G {input.fasta} -a {input.gtf} -o {output} {input.bam} 2> {log} 
        """

rule multiqc_feature:
    input:
        "results/feature/featureCounts_results.txt.summary"
    output:
        "results/feature/feature_multiqc.html"
    log:
        "results/feature/logs/multiqc.log"
    wrapper:
        "0.49.0/bio/multiqc"