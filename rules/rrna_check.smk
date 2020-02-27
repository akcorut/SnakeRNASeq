rule bwa_mem:
    input:
        r1="results/trimmed/{smp}_R1_val_1.fq.gz",
        r2="results/trimmed/{smp}_R2_val_2.fq.gz",
        rrna="/work/jawlab/kivanc/PeanutRnaSeq/rRNA/rrna.fa"
    output:
        temp("results/rrnaCheck/{smp}_rrna.sam")
    threads:16
    conda:
        "../envs/bwa.yaml"
    shell:
        '''
        bwa mem -t {threads} {input.rrna} {input.r1} {input.r2} > {output}
        '''

rule sam_to_bam:
    input:
        sam="results/rrnaCheck/{smp}_rrna.sam"
    output:
        temp("results/rrnaCheck/bam/{smp}_rrna.bam")
    threads:16
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        samtools view -@ {threads} -bS -o {output.bam} {input.sam}
        '''

rule flagstat:
    input:
        bam="results/rrnaCheck/bam/{smp}_rrna.bam"
    output:
        "results/rrnaCheck/flagstat/{smp}_rrna.out"
    threads:16
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        samtools flagstat -@ {threads} {input.bam} > {output}
        '''

rule stats:
    input:
        bam="results/rrnaCheck/bam/{smp}_rrna.bam"
    output:
        "results/rrnaCheck/stats/{smp}_rrna_stats.out"
    threads:16
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        samtools stats -@ {threads} {input.bam} > {output}
        '''

rule multiqc_rrna:
    input:
        expand("results/rrnaCheck/flagstat/{smp}_rrna.out", smp=sample_id),
        expand("results/rrnaCheck/bam/{smp}_rrna.bam", smp=sample_id),
        expand("results/rrnaCheck/stats/{smp}_rrna_stats.out", smp=sample_id)
    output:
        "results/rrnaCheck/rrna_multiqc_report.html"
    log:
        "results/rrnaCheck/logs/multiqc.log"
    priority:-10
    wrapper:
        "0.48.0/bio/multiqc"


############################ After Decontamination #####################################

rule bwa_mem_decontamination:
    input:
        r1="results/bbsplit/{smp}_R1_clean.fq",
        r2="results/bbsplit/{smp}_R2_clean.fq",
        rrna="/work/jawlab/kivanc/PeanutRnaSeq/rRNA/rrna.fa"
    output:
        temp("results/rrna_decont/{smp}_clean_rrna.sam")
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
        sam="results/rrna_decont/{smp}_clean_rrna.sam"
    output:
        temp("results/rrna_decont/bam/{smp}_clean_rrna.bam")
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
        bam="results/rrna_decont/bam/{smp}_clean_rrna.bam"
    output:
        "results/rrna_decont/stats/{smp}_clean_rrna_stats.out"
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
        expand("results/rrna_decont/bam/{smp}_clean_rrna.bam", smp=sample_id),
        expand("results/rrna_decont/stats/{smp}_clean_rrna_stats.out", smp=sample_id)
    output:
        "results/rrna_decont/rrna_clean_multiqc_report.html"
    log:
        "results/rrna_decont/logs/multiqc.log"
    priority:-10
    wrapper:
        "0.48.0/bio/multiqc"