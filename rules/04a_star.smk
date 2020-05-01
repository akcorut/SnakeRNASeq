rule gff3_to_gtf:
    input:
        anno = ANNOTATION
    output:
        gtf = "/work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf"
    priority:100
    shell:
        "gffread {input.anno} -T -o {output.gtf}"

rule star_index:
    input:
        fasta = REFERENCE
    output:
        directory(config["ref"]["star_index"])
    threads:15
    priority:50
    params:
        extra = "",
        gtf = rules.gff3_to_gtf.output.gtf
    log:
        config["ref"]["star_index"] + "/log/star_index_.log"
    wrapper:
        "0.49.0/bio/star/index"

rule star_pass1:
    input:
        fq1=GetClean(0),
        fq2=GetClean(1)
    output:
        "results/04_alignment/04a_alignment_results/star/pass1/{smp}/Aligned.out.bam",
        "results/04_alignment/04a_alignment_results/star/pass1/{smp}/SJ.out.tab"
    log:
        "results/04_alignment/04a_alignment_results/star/pass1/logs/{smp}.log"
    priority:10
    params:
        # path to STAR reference genome index
        index=config["ref"]["star_index"],
        extra="--outSAMtype BAM Unsorted --alignIntronMax 10000 --sjdbGTFfile {}".format(
              rules.gff3_to_gtf.output.gtf)
    threads:20
    wrapper:
        "0.49.0/bio/star/align"

rule get_junctions:
    input:
        expand("results/04_alignment/04a_alignment_results/star/pass1/{smp}/SJ.out.tab", smp=sample_id)
    output:
        sj="results/04_alignment/04a_alignment_results/star/junctions/SJ.filtered.tab"
    priority:1
    shell:
        """
        cat {input} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output}
        """

rule star_pass2:
    input:
        fq1=GetClean(0),
        fq2=GetClean(1),
    output:
        "results/04_alignment/04a_alignment_results/star/pass2/{smp}/Aligned.out.bam",
        "results/04_alignment/04a_alignment_results/star/pass2/{smp}/Aligned.toTranscriptome.out.bam",
        "results/04_alignment/04a_alignment_results/star/pass2/{smp}/Aligned.sortedByCoord.out.bam",
        "results/04_alignment/04a_alignment_results/star/pass2/{smp}/ReadsPerGene.out.tab"
    log:
        "results/04_alignment/04a_alignment_results/star/pass2/logs/{smp}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["star_index"],
        extra="--outSAMunmapped Within --outSAMtype BAM SortedByCoordinate Unsorted --quantMode GeneCounts TranscriptomeSAM --alignIntronMax 10000 --sjdbFileChrStartEnd {} --sjdbGTFfile {}".format(
              rules.get_junctions.output.sj, rules.gff3_to_gtf.output.gtf)
    priority:-1
    threads:20
    wrapper:
        "0.49.0/bio/star/align"