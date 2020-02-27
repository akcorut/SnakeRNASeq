rule kallisto_index_V2:
    input:
        tcp = config["ref"]["transcript2"]
    output:
        "results/kallistoV2/index/transcripts.idx"
    log:
        "results/kallistoV2/logs/index.log"
    priority:50
    conda:
        "../envs/kallisto.yaml"
    threads:20
    shell:
        "kallisto index -i {output} {input} 2> {log}"

rule kallisto_quant_V2:
    input:
        r1="results/bbsplit/{smp}_R1_clean.fq",
        r2="results/bbsplit/{smp}_R2_clean.fq",
        idx= "results/kallistoV2/index/transcripts.idx",
        gtf= config["ref"]["annotation2"]
    output:
        directory("results/kallistoV2/quant/{smp}")
    log:
        "results/kallistoV2/quant/{smp}.log"
    params:
        extra="--rf-stranded -b 100 --bias"
    conda:
        "../envs/kallisto.yaml"
    threads:20
    shell:
        "kallisto quant --threads {threads} -i {input.idx} -o {output} --gtf {input.gtf} "
        "{params.extra} {input.r1} {input.r2} 2> {log}"