rule bwa_mem:
    input:
        r1="results/02_trim/{smp}_R1_val_1.fq.gz",
        r2="results/02_trim/{smp}_R2_val_2.fq.gz",
        rrna=config["rRNA"]
    output:
        temp("results/03_decontamination/03a_rrna_check/{smp}_rrna.sam")
    threads:16
    conda:
        "../envs/bwa.yaml"
    priority:5
    shell:
        '''
        bwa mem -t {threads} {input.rrna} {input.r1} {input.r2} > {output}
        '''

rule sam_to_bam:
    input:
        sam="results/03_decontamination/03a_rrna_check/{smp}_rrna.sam"
    output:
        temp("results/03_decontamination/03a_rrna_check/bam/{smp}_rrna.bam")
    threads:16
    conda:
        "../envs/samtools.yaml"
    priority:4
    shell:
        '''
        samtools view -@ {threads} -bS -o {output} {input.sam}
        '''

rule flagstat:
    input:
        bam="results/03_decontamination/03a_rrna_check/bam/{smp}_rrna.bam"
    output:
        "results/03_decontamination/03a_rrna_check/flagstat/{smp}_rrna.out"
    threads:16
    conda:
        "../envs/samtools.yaml"
    priority:3
    shell:
        '''
        samtools flagstat -@ {threads} {input.bam} > {output}
        '''

rule stats:
    input:
        bam="results/03_decontamination/03a_rrna_check/bam/{smp}_rrna.bam"
    output:
        "results/03_decontamination/03a_rrna_check/stats/{smp}_rrna_stats.out"
    threads:16
    priority:2
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        samtools stats -@ {threads} {input.bam} > {output}
        '''

rule multiqc_rrna:
    input:
        expand("results/03_decontamination/03a_rrna_check/flagstat/{smp}_rrna.out", smp=sample_id),
        expand("results/03_decontamination/03a_rrna_check/bam/{smp}_rrna.bam", smp=sample_id),
        expand("results/03_decontamination/03a_rrna_check/stats/{smp}_rrna_stats.out", smp=sample_id)
    output:
        "results/03_decontamination/03a_rrna_check/rrna_multiqc_report.html"
    log:
        "results/03_decontamination/03a_rrna_check/logs/multiqc.log"
    priority:1
    wrapper:
        "0.48.0/bio/multiqc"