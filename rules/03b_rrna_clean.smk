rule bbmap:
    input:
        r1="results/02_trim/{smp}_R1_val_1.fq.gz",
        r2="results/02_trim/{smp}_R2_val_2.fq.gz",
        rrna=config["rRNA"]
    output:
        out1="results/03_decontamination/03b_rrna_cleaned/{smp}_R1_clean.fq",
        out2="results/03_decontamination/03b_rrna_cleaned/{smp}_R2_clean.fq",
        stats="results/03_decontamination/03b_rrna_cleaned/stats/{smp}_bbsplit_stats.txt"
    params:
        out_rrna="results/03_decontamination/03b_rrna_cleaned/rrna/{smp}",
        index= config["bb_index"]
    log:
        "results/03_decontamination/03b_rrna_cleaned/logs/{smp}_decontamination.log"
    conda:
        "../envs/bbmap.yaml"
    priority:3
    threads:16
    shell:
        """
        bbsplit.sh -Xmx120g threads={threads} in1={input.r1} in2={input.r2} ref_rrna={input.rrna} path={params.index} basename={params.out_rrna}/out_%.fq outu1={output.out1} outu2={output.out2} refstats={output.stats} 2> {log}
        """

rule fastqc_clean:
    input:
        r1=rules.bbmap.output.out1,
        r2=rules.bbmap.output.out2
    output:
        "results/01_qc/01a_fqc/fqc_clean/{smp}_R1_clean_fastqc.html",
        "results/01_qc/01a_fqc/fqc_clean/{smp}_R1_clean_fastqc.zip",
        "results/01_qc/01a_fqc/fqc_clean/{smp}_R2_clean_fastqc.html",
        "results/01_qc/01a_fqc/fqc_clean/{smp}_R2_clean_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    priority:2
    threads:20
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} -q -f fastq -o results/01_qc/01a_fqc/fqc_clean/
        """

rule multiqc_clean:
    input:
        expand("results/01_qc/01a_fqc/fqc_clean/{smp}_R1_clean_fastqc.html", smp=sample_id),
        expand("results/01_qc/01a_fqc/fqc_clean/{smp}_R2_clean_fastqc.html", smp=sample_id)
    output:
        "results/01_qc/01b_multiqc/multiqc_clean/multiqc_report_decont.html"
    priority:1
    log:
        "results/01_qc/01b_multiqc/multiqc_clean/logs/multiqc.log"
    wrapper:
        "0.47.0/bio/multiqc"