rule sort_bam:
    input:
        bam= rules.star_pass2.output.sorted_bam
    output:
        temp("results/04_alignment/04b_alignment_qc/qualimap/SortedBam/{smp}_sorted.bam")
    params:
        prfx=lambda wildcards, output: output[0][:-4]
    priority:50
    conda:
        "../envs/samtools.yaml"
    priority:50
    threads:20
    shell:
        '''
        samtools sort -@ {threads} -n -o {output} -T {params.prfx} {input.bam}
        '''

rule qualimap:
    input: 
        sorted_bam= rules.sort_bam.output
    output: 
        "results/04_alignment/04b_alignment_qc/qualimap/{smp}/qualimapReport.html"
    params:
        outdir=lambda wildcards, output: output[0][:-20],
        gtf= rules.gff3_to_gtf.output.gtf
    conda:
        "../envs/qualimap.yaml"
    shell:
        '''
        qualimap rnaseq -bam {input.sorted_bam} -gtf {params.gtf} --outdir {params.outdir} --sorted --paired
        '''

rule multiqc_qualimap:
    input:
        expand("results/04_alignment/04b_alignment_qc/qualimap{smp}/qualimapReport.html", smp=sample_id)
    output:
        "results/04_alignment/04b_alignment_qc/qualimap/qualimap_multiqc.html"
    priority:-1
    log:
        "results/04_alignment/04b_alignment_qc/qualimap/logs/multiqc.log"
    wrapper:
        "0.49.0/bio/multiqc"

ruleorder: sort_bam > qualimap > multiqc_qualimap