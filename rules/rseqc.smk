rule gff2bed:
    input:
        anno=ANNOTATION
    output:
        bed="results/rseqc/tifrunner_annotation.bed",
    priority: 2
    log:
        "results/rseqc/logs/gtf2bed.log"
    conda:
        "../envs/bedops.yaml"
    script:
        "gff2bed < {input.anno} > {output.bed}"

rule rseqc_junction_annotation:
    input:
        bam="results/star/{smp}/Aligned.out.bam",
        bed=rules.gff2bed.output.bed
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

        
rule rseqc_junction_saturation:
    input:
        bam="results/star/{smp}/Aligned.out.bam",
        bed=rules.gff2bed.output.bed
    output:
        "results/rseqc/{smp}.junctionsat.junctionSaturation_plot.pdf"
    priority: 1
    log:
        "results/rseqc/logs/rseqc_junction_saturation/{smp}.log"
    params:
        extra=r"-q 255", 
        prefix="results/rseqc/{smp}.junctionsat"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


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
        bed=rules.gff2bed.output.bed
    output:
        "results/rseqc/{smp}.infer_experiment.txt"
    priority: 1
    log:
        "results/rseqc/logs/rseqc_infer/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

        
rule rseqc_innerdis:
    input:
        bam="results/star/{smp}/Aligned.out.bam",
        bed=rules.gff2bed.output.bed
    output:
        "results/rseqc/{smp}.inner_distance_freq.inner_distance.txt"
    priority: 1
    log:
        "results/rseqc/logs/rseqc_innerdis/{smp}.log"
    params:
        prefix="results/rseqc/{smp}.inner_distance_freq"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="results/star/{smp}/Aligned.out.bam",
        bed=rules.gff2bed.output.bed
    output:
        "results/rseqc/{smp}.readdistribution.txt"
    priority: 1
    log:
        "results/rseqc/logs/rseqc_readdis/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "results/star/{smp}/Aligned.out.bam"
    output:
        "results/rseqc/{smp}.readdup.DupRate_plot.pdf"
    priority: 1
    log:
        "results/rseqc/logs/rseqc_readdup/{smp}.log"
    params:
        prefix="results/rseqc/{smp}.readdup"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"

        
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
        

rule multiqc:
    input:
        expand("results/star/{smp}/Aligned.out.bam", smp=sample_id),
        expand("results/rseqc/{smp}.junctionanno.junction.bed", smp=sample_id),
        expand("results/rseqc/{smp}.junctionsat.junctionSaturation_plot.pdf", smp=sample_id),
        expand("results/rseqc/{smp}.infer_experiment.txt", smp=sample_id),
        expand("results/rseqc/{smp}.stats.txt", smp=sample_id),
        expand("results/rseqc/{smp}.inner_distance_freq.inner_distance.txt", smp=sample_id),
        expand("results/rseqc/{smp}.readdistribution.txt", smp=sample_id),
        expand("results/rseqc/{smp}.readdup.DupRate_plot.pdf", smp=sample_id),
        expand("results/rseqc/{smp}.readgc.GC_plot.pdf", smp=sample_id),
        expand("results/rseqc/logs/rseqc_junction_annotation/{smp}.log", smp=sample_id)
    output:
        "results/rseqc/multiqc_report.html"
    log:
        "results/rseqc/logs/multiqc.log"
    wrapper:
"0.31.1/bio/multiqc"