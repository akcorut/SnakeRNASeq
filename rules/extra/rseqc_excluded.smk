rule rseqc_innerdis:
    input:
        bam="results/star/{smp}/Aligned.out.bam",
        bed=rules.gtf2bed.output.bed
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

rule rseqc_junction_saturation:
    input:
        bam="results/star/{smp}/Aligned.out.bam",
        bed=rules.gtf2bed.output.bed
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

