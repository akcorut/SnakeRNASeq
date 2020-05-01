rule bwa_mem:
    input:
        r1="results/02_trim/{smp}_R1_val_1.fq.gz",
        r2="results/02_trim/{smp}_R2_val_2.fq.gz",
        rrna=config["files"]["rRNA"]
    output:
        temp("results/03_decontamination/03a_rrna_check/{smp}_rrna.sam")
    threads:16
    conda:
        "../envs/bwa.yaml"
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
    shell:
        '''
        samtools view -@ {threads} -bS -o {output.bam} {input.sam}
        '''

rule flagstat:
    input:
        bam="results/03_decontamination/03a_rrna_check/bam/{smp}_rrna.bam"
    output:
        "results/03_decontamination/03a_rrna_check/flagstat/{smp}_rrna.out"
    threads:16
    conda:
        "../envs/samtools.yaml"
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
    priority:-10
    wrapper:
        "0.48.0/bio/multiqc"


############################ After Decontamination #####################################


rule bwa_mem_decontamination:
    input:
        r1="results/03_decontamination/03b_rrna_cleaned/{smp}_R1_clean.fq",
        r2="results/03_decontamination/03b_rrna_cleaned/{smp}_R2_clean.fq",
        rrna=config["files"]["rRNA"]
    output:
        temp("results/03_decontamination/03c_rrna_cleaned_check/{smp}_clean_rrna.sam")
    threads:16
    priority:50
    conda:
        "../envs/bwa.yaml"
    shell:
        '''
        bwa mem -t {threads} {input.rrna} {input.r1} {input.r2} > {output}
        '''

rule sam_to_bam_decontamination:
    input:
        sam="results/03_decontamination/03c_rrna_cleaned_check/{smp}_clean_rrna.sam"
    output:
        temp("results/03_decontamination/03c_rrna_cleaned_check/bam/{smp}_clean_rrna.bam")
    threads:16
    priority:10
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        samtools view -@ {threads} -bS -o {output.bam} {input.sam}
        '''

rule stats_decontamination:
    input:
        bam="results/03_decontamination/03c_rrna_cleaned_check/bam/{smp}_clean_rrna.bam"
    output:
        "results/03_decontamination/03c_rrna_cleaned_check/stats/{smp}_clean_rrna_stats.out"
    threads:16
    priority:1
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        samtools stats -@ {threads} {input.bam} > {output}
        '''

rule multiqc_rrna_decontamination:
    input:
        expand("results/03_decontamination/03c_rrna_cleaned_check/bam/{smp}_clean_rrna.bam", smp=sample_id),
        expand("results/03_decontamination/03c_rrna_cleaned_check/stats/{smp}_clean_rrna_stats.out", smp=sample_id)
    output:
        "results/03_decontamination/03c_rrna_cleaned_check/rrna_clean_multiqc_report.html"
    log:
        "results/03_decontamination/03c_rrna_cleaned_check/logs/multiqc.log"
    priority:-10
    wrapper:
        "0.48.0/bio/multiqc"