rule rseqc_junction_annotation_h2:
    input:
        bam="results/hisat2/{smp}.cutadapt.bam",
        bed=rules.gtf2bed.output.bed
    output:
        "results/rseqc/hisat2/{smp}.junctionanno.junction.bed"
    priority: 1
    log:
        "results/rseqc/hisat2/logs/rseqc_junction_annotation/{smp}.log"
    params:
        extra=r"-q 40",  # hisat2 seems to be uses 40
        prefix="results/rseqc/hisat2/{smp}.junctionanno"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"

        
rule rseqc_junction_saturation_h2:
    input:
        bam="results/hisat2/{smp}.cutadapt.bam",
        bed=rules.gtf2bed.output.bed
    output:
        "results/rseqc/hisat2/{smp}.junctionsat.junctionSaturation_plot.pdf"
    priority: 1
    log:
        "results/rseqc/hisat2/logs/rseqc_junction_saturation/{smp}.log"
    params:
        extra=r"-q 40", 
        prefix="results/rseqc/hisat2/{smp}.junctionsat"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat_h2:
    input:
        "results/hisat2/{smp}.cutadapt.bam"
    output:
        "results/rseqc/hisat2/{smp}.stats.txt"
    priority: 1
    log:
        "results/rseqc/hisat2/logs/rseqc_stat/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

        
rule rseqc_infer_h2:
    input:
        bam="results/hisat2/{smp}.cutadapt.bam",
        bed=rules.gtf2bed.output.bed
    output:
        "results/rseqc/hisat2/{smp}.infer_experiment.txt"
    priority: 1
    log:
        "results/rseqc/hisat2/logs/rseqc_infer/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

        
rule rseqc_innerdis_h2:
    input:
        bam="results/hisat2/{smp}.cutadapt.bam",
        bed=rules.gtf2bed.output.bed
    output:
        "results/rseqc/hisat2/{smp}.inner_distance_freq.inner_distance.txt"
    priority: 1
    log:
        "results/rseqc/hisat2/logs/rseqc_innerdis/{smp}.log"
    params:
        prefix="results/rseqc/hisat2/{smp}.inner_distance_freq"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis_h2:
    input:
        bam="results/hisat2/{smp}.cutadapt.bam",
        bed=rules.gtf2bed.output.bed
    output:
        "results/rseqc/hisat2/{smp}.readdistribution.txt"
    priority: 1
    log:
        "results/rseqc/hisat2/logs/rseqc_readdis/{smp}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup_h2:
    input:
        "results/hisat2/{smp}.cutadapt.bam"
    output:
        "results/rseqc/hisat2/{smp}.readdup.DupRate_plot.pdf"
    priority: 1
    log:
        "results/rseqc/hisat2/logs/rseqc_readdup/{smp}.log"
    params:
        prefix="results/rseqc/hisat2/{smp}.readdup"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"

        
rule rseqc_readgc_h2:
    input:
        "results/hisat2/{smp}.cutadapt.bam"
    output:
        "results/rseqc/hisat2/{smp}.readgc.GC_plot.pdf"
    priority: 1
    log:
        "results/rseqc/hisat2/logs/rseqc_readgc/{smp}.log"
    params:
        prefix="results/rseqc/hisat2/{smp}.readgc"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"
        

rule multiqc_rseqc_h2:
    input:
        expand("results/hisat2/{smp}.cutadapt.bam", smp=sample_id),
        expand("results/rseqc/hisat2/{smp}.junctionanno.junction.bed", smp=sample_id),
        expand("results/rseqc/hisat2/{smp}.junctionsat.junctionSaturation_plot.pdf", smp=sample_id),
        expand("results/rseqc/hisat2/{smp}.infer_experiment.txt", smp=sample_id),
        expand("results/rseqc/hisat2/{smp}.stats.txt", smp=sample_id),
        expand("results/rseqc/hisat2/{smp}.inner_distance_freq.inner_distance.txt", smp=sample_id),
        expand("results/rseqc/hisat2/{smp}.readdistribution.txt", smp=sample_id),
        expand("results/rseqc/hisat2/{smp}.readdup.DupRate_plot.pdf", smp=sample_id),
        expand("results/rseqc/hisat2/{smp}.readgc.GC_plot.pdf", smp=sample_id),
        expand("results/rseqc/hisat2/logs/rseqc_junction_annotation/{smp}.log", smp=sample_id)
    output:
        "results/rseqc/hisat2/multiqc_rseqc_h2_report.html"
    log:
        "results/rseqc/hisat2/logs/multiqc.log"
    wrapper:
        "0.31.1/bio/multiqc"