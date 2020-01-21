rule gtf2bed:
    input:
        rules.gff3_to_gtf.output.gtf
    output:
        bed="results/rseqc/tifrunner_annotation.bed",
        db=temp("results/rseqc/tifrunner_annotation.db")
    priority: 50
    log:
        "results/rseqc/logs/gtf2bed.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"

rule rseqc_junction_annotation:
    input:
        bam="results/star/{smp}/Aligned.out.bam",
        bed=rules.gtf2bed.output.bed
    output:
        "results/rseqc/{smp}.junctionanno.junction.bed"
    priority: 1
    log:
        "results/rseqc/logs/rseqc_junction_annotation/{smp}.log"
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix="results/rseqc/{smp}.junctionanno"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"

rule rseqc_stat:
    input:
        "results/star/{smp}/Aligned.out.bam"
    output:
        "results/rseqc/{smp}.stats.txt"
    priority: 1
    log:
        "results/rseqc/logs/rseqc_stat/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

        
rule rseqc_infer:
    input:
        bam="results/star/{smp}/Aligned.out.bam",
        bed=rules.gtf2bed.output.bed
    output:
        "results/rseqc/{smp}.infer_experiment.txt"
    priority: 1
    log:
        "results/rseqc/logs/rseqc_infer/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

rule rseqc_readdis:
    input:
        bam="results/star/{smp}/Aligned.out.bam",
        bed=rules.gtf2bed.output.bed
    output:
        "results/rseqc/{smp}.readdistribution.txt"
    priority: 1
    log:
        "results/rseqc/logs/rseqc_readdis/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"
        
rule rseqc_readgc:
    input:
        "results/star/{smp}/Aligned.out.bam"
    output:
        "results/rseqc/{smp}.readgc.GC_plot.pdf"
    priority: 1
    log:
        "results/rseqc/logs/rseqc_readgc/{smp}.log"
    params:
        prefix="results/rseqc/{smp}.readgc"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"

#rule rseqc_readqual:
#    input:
#        "results/star/{smp}/Aligned.out.bam"
#    output:
#        "results/rseqc/{smp}.readqual.heatmap.pdf",
#        "results/rseqc/{smp}.readqual.boxplot.pdf"
#    priority: 1
#    log:
#        "results/rseqc/logs/rseqc_readqual/{smp}.log"
#    params:
#        prefix="results/rseqc/{smp}.readqual"
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "read_quality.py -i {input} -o {params.prefix} > {log} 2>&1"

#rule rseqc_infer_hisat2:
#    input:
#        bam="results/hisat2/{smp}.trimmed.bam",
#        bed=rules.gtf2bed.output.bed
#    output:
#        "results/rseqc/hisat2/{smp}.infer_experiment.txt"
#    priority: 1
#    log:
#        "results/rseqc/hisat2/logs/rseqc_infer/{smp}.log"
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

rule multiqc_rseqc:
    input:
        expand("results/star/{smp}/Aligned.out.bam", smp=sample_id),
        expand("results/rseqc/{smp}.junctionanno.junction.bed", smp=sample_id),
        expand("results/rseqc/{smp}.infer_experiment.txt", smp=sample_id),
        expand("results/rseqc/{smp}.stats.txt", smp=sample_id),
        expand("results/rseqc/{smp}.readdistribution.txt", smp=sample_id),
        expand("results/rseqc/{smp}.readgc.GC_plot.pdf", smp=sample_id),
        expand("results/rseqc/logs/rseqc_junction_annotation/{smp}.log", smp=sample_id),
        #expand("results/rseqc/{smp}.readqual.heatmap.pdf", smp=sample_id),
        #expand("results/rseqc/{smp}.readqual.boxplot.pdf", smp=sample_id),
        expand("results/star/{smp}/Log.final.out", smp=sample_id)
    output:
        "results/rseqc/multiqc_report.html"
    log:
        "results/rseqc/logs/multiqc.log"
    wrapper:
        "0.31.1/bio/multiqc"