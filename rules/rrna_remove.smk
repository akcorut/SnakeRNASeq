rule bbmap:
    input:
        r1="results/trimmed/{smp}_R1_val_1.fq.gz",
        r2="results/trimmed/{smp}_R2_val_2.fq.gz",
        rrna="/work/jawlab/kivanc/PeanutRnaSeq/rRNA/rrna.fa"
    output:
        out1="results/bbsplit/{smp}_R1_clean.fq",
        out2="results/bbsplit/{smp}_R2_clean.fq",
        stats="results/bbsplit/stats/{smp}_bbsplit_stats.txt"
    params:
        out_rrna="results/bbsplit/rrna/{smp}",
        index="/work/jawlab/kivanc/PeanutRnaSeq/rRNA/bb_index"
    log:
        "results/bbsplit/logs/{smp}_decontamination.log"
    conda:
        "../envs/bbmap.yaml"
    threads:16
    shell:
        """
        bbsplit.sh -Xmx120g threads={threads} in1={input.r1} in2={input.r2} ref_rrna={input.rrna} path={params.index} basename={params.out_rrna}/out_%.fq outu1={output.out1} outu2={output.out2} refstats={output.stats} 2> {log}
        """