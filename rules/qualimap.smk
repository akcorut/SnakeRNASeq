rule sort_bam:
    input:
        bam= "results/star/{smp}/Aligned.sortedByCoord.out.bam"
    output:
        temp("results/SortedBam/{smp}_sorted.bam")
    params:
        prfx="/work/jawlab/kivanc/tmp/{smp}_sorted"
    priority:50
    conda:
        "../envs/samtools.yaml"
    priority:50
    threads:8
    shell:
        '''
        samtools sort -@ {threads} -n -o {output} -T {params.prfx} {input.bam}
        '''

rule qualimap:
    input: 
        sorted_bam= "results/SortedBam/{smp}_sorted.bam"
    output: 
        "results/qualimap/{smp}/qualimapReport.html"
    params:
        outdir="results/qualimap/{smp}",
        gtf= rules.gff3_to_gtf.output.gtf
    conda:
        "../envs/qualimap.yaml"
    shell:
        '''
        qualimap rnaseq -bam {input.sorted_bam} -gtf {params.gtf} --outdir {params.outdir} --sorted --paired
        '''