rule kallisto_index:
    input:
        tcp= TRANSCRIPTS
    output:
        "results/06_alignment_free/06b_kallisto/index/transcripts.idx"
    log:
        "results/06_alignment_free/06b_kallisto/logs/index.log"
    priority:3
    conda:
        "../envs/kallisto.yaml"
    threads:20
    shell:
        "kallisto index -i {output} {input} 2> {log}"

rule kallisto_quant:
    input:
        r1=GetClean(0),
        r2=GetClean(1),
        idx= rules.kallisto_index.output,
        gtf= rules.gff3_to_gtf.output.gtf
    output:
        directory("results/06_alignment_free/06b_kallisto/quant/{smp}")
    log:
        "results/06_alignment_free/06b_kallisto/quant/{smp}.log"
    params:
        extra="--rf-stranded -b 100 --bias"
    conda:
        "../envs/kallisto.yaml"
    threads:20
    priority:2
    shell:
        "kallisto quant --threads {threads} -i {input.idx} -o {output} --gtf {input.gtf} "
        "{params.extra} {input.r1} {input.r2} 2> {log}"

rule multiqc_kallisto:
    input:
        expand("results/06_alignment_free/06b_kallisto/quant/{smp}.log", smp=sample_id)
    output:
        "results/06_alignment_free/06b_kallisto/kallisto_multiqc.html"
    log:
        "results/06_alignment_free/06b_kallisto/logs/multiqc.log"
    priority:1
    wrapper:
        "0.49.0/bio/multiqc"