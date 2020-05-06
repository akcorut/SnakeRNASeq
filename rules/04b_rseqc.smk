rule gtf2bed:
    input:
        rules.gff3_to_gtf.output.gtf
    output:
        bed="results/04_alignment/04b_alignment_qc/rseqc/tifrunner_annotation.bed",
        db=temp("results/04_alignment/04b_alignment_qc/rseqc/tifrunner_annotation.db")
    priority:7
    log:
        "results/04_alignment/04b_alignment_qc/rseqc/logs/gtf2bed.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"

rule rseqc_junction_annotation:
    input:
        bam=rules.star_pass2.output.bam,
        bed=rules.gtf2bed.output.bed
    output:
        "results/04_alignment/04b_alignment_qc/rseqc/{smp}.junctionanno.junction.bed"
    priority:6
    log:
        "results/04_alignment/04b_alignment_qc/rseqc/logs/rseqc_junction_annotation/{smp}.log"
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda wildcards, output: output[0][:-13]
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"

rule rseqc_stat:
    input:
        rules.star_pass2.output.bam
    output:
        "results/04_alignment/04b_alignment_qc/rseqc/{smp}.stats.txt"
    priority:5
    log:
        "results/04_alignment/04b_alignment_qc/rseqc/logs/rseqc_stat/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

        
rule rseqc_infer:
    input:
        bam=rules.star_pass2.output.bam,
        bed=rules.gtf2bed.output.bed
    output:
        "results/04_alignment/04b_alignment_qc/rseqc/{smp}.infer_experiment.txt"
    priority:4
    log:
        "results/04_alignment/04b_alignment_qc/rseqc/logs/rseqc_infer/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

rule rseqc_readdis:
    input:
        bam=rules.star_pass2.output.bam,
        bed=rules.gtf2bed.output.bed
    output:
        "results/04_alignment/04b_alignment_qc/rseqc/{smp}.readdistribution.txt"
    priority:3
    log:
        "results/04_alignment/04b_alignment_qc/rseqc/logs/rseqc_readdis/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"
        
rule rseqc_readgc:
    input:
        rules.star_pass2.output.bam
    output:
        "results/04_alignment/04b_alignment_qc/rseqc/{smp}.readgc.GC.xls"
    priority:2
    log:
        "results/04_alignment/04b_alignment_qc/rseqc/logs/rseqc_readgc/{smp}.log"
    params:
        prefix=lambda wildcards, output: output[0][:-7]
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"

#rule rseqc_readqual:
#    input:
#        "results/04_alignment/04a_alignment_results/star/{smp}/Aligned.out.bam"
#    output:
#        "results/04_alignment/04b_alignment_qc/rseqc/{smp}.readqual.heatmap.pdf",
#        "results/04_alignment/04b_alignment_qc/rseqc/{smp}.readqual.boxplot.pdf"
#    priority: 1
#    log:
#        "results/04_alignment/04b_alignment_qc/rseqc/logs/rseqc_readqual/{smp}.log"
#    params:
#        prefix="results/04_alignment/04b_alignment_qc/rseqc/{smp}.readqual"
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "read_quality.py -i {input} -o {params.prefix} > {log} 2>&1"


rule multiqc_rseqc:
    input:
        expand("results/04_alignment/04a_alignment_results/star/pass2/{smp}/Aligned.out.bam", smp=sample_id),
        expand("results/04_alignment/04b_alignment_qc/rseqc/{smp}.junctionanno.junction.bed", smp=sample_id),
        expand("results/04_alignment/04b_alignment_qc/rseqc/{smp}.infer_experiment.txt", smp=sample_id),
        expand("results/04_alignment/04b_alignment_qc/rseqc/{smp}.stats.txt", smp=sample_id),
        expand("results/04_alignment/04b_alignment_qc/rseqc/{smp}.readdistribution.txt", smp=sample_id),
        expand("results/04_alignment/04b_alignment_qc/rseqc/{smp}.readgc.GC.xls", smp=sample_id),
        expand("results/04_alignment/04b_alignment_qc/rseqc/logs/rseqc_junction_annotation/{smp}.log", smp=sample_id),
        #expand("results/04_alignment/04b_alignment_qc/rseqc/{smp}.readqual.heatmap.pdf", smp=sample_id),
        #expand("results/04_alignment/04b_alignment_qc/rseqc/{smp}.readqual.boxplot.pdf", smp=sample_id),
        expand("results/04_alignment/04a_alignment_results/star/pass2/{smp}/Log.final.out", smp=sample_id)
    output:
        "results/04_alignment/04b_alignment_qc/rseqc/rseqc_multiqc_report.html"
    priority:1
    log:
        "results/04_alignment/04b_alignment_qc/rseqc/logs/multiqc.log"
    wrapper:
        "0.49.0/bio/multiqc"