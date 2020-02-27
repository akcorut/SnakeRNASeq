rule kallisto_index:
    input:
        tcp= TRANSCRIPTS
    output:
        "results/kallisto/index/transcripts.idx"
    log:
        "results/kallisto/logs/index.log"
    priority:50
    conda:
        "../envs/kallisto.yaml"
    threads:20
    shell:
        "kallisto index -i {output} {input} 2> {log}"

rule kallisto_quant:
    input:
        r1=GetClean(0),
        r2=GetClean(1),
        idx= "results/kallisto/index/transcripts.idx",
        gtf= rules.gff3_to_gtf.output.gtf
    output:
        directory("results/kallisto/quant/{smp}")
    log:
        "results/kallisto/quant/{smp}.log"
    params:
        extra="--rf-stranded -b 100 --bias"
    conda:
        "../envs/kallisto.yaml"
    threads:20
    shell:
        "kallisto quant --threads {threads} -i {input.idx} -o {output} --gtf {input.gtf} "
        "{params.extra} {input.r1} {input.r2} 2> {log}"