rule htseq_count:
    input:
        rules.star_pass2.output.sorted_bam
    output:
        "results/05_quantification/05b_htseq/{smp}_htseq.cnt"
    conda:
        "../envs/htseq.yaml"
    params:
        gtf= rules.gff3_to_gtf.output.gtf
    log:
        "results/05_quantification/05b_htseq/{smp}_htseq_count.log"
    threads: 24
    priority:2
    shell:
        """
        htseq-count -m intersection-nonempty --stranded=reverse --idattr gene_id -r pos -f bam {input} {params.gtf} > {output} 2> {log}
        """

rule htseq_multiqc:
    input:
        expand("results/05_quantification/05b_htseq/{smp}_htseq.cnt", smp=sample_id)
    output:
        "results/05_quantification/05b_htseq/htseq_multiqc.html"
    log:
        "results/05_quantification/05b_htseq/logs/multiqc.log"
    params:
        prefix = lambda wildcards, output: output[0][:-18]
    conda:
        "../envs/multiqc.yaml"
    priority:1
    shell:
        """
        multiqc -m htseq {params.prefix} --filename {output} 2> {log}
        """
