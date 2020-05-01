rule bbmap:
    input:
        r1="results/02_trim/{smp}_R1_val_1.fq.gz",
        r2="results/02_trim/{smp}_R2_val_2.fq.gz",
        rrna=config["files"]["rRNA"]
    output:
        out1="results/03_decontamination/03b_rrna_cleaned/{smp}_R1_clean.fq",
        out2="results/03_decontamination/03b_rrna_cleaned/{smp}_R2_clean.fq",
        stats="results/03_decontamination/03b_rrna_cleaned/stats/{smp}_bbsplit_stats.txt"
    params:
        out_rrna="results/03_decontamination/03b_rrna_cleaned/rrna/{smp}",
        index="/work/jawlab/kivanc/PeanutRnaSeq/rRNA/bb_index"
    log:
        "results/03_decontamination/03b_rrna_cleaned/logs/{smp}_decontamination.log"
    conda:
        "../envs/bbmap.yaml"
    threads:16
    shell:
        """
        bbsplit.sh -Xmx120g threads={threads} in1={input.r1} in2={input.r2} ref_rrna={input.rrna} path={params.index} basename={params.out_rrna}/out_%.fq outu1={output.out1} outu2={output.out2} refstats={output.stats} 2> {log}
        """