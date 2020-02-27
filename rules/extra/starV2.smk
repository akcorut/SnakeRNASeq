rule gff3_to_gtf:
    input:
        anno = ANNOTATION
    output:
        gtf = "/work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf"
    shell:
        "gffread {input.anno} -T -o {output.gtf}"

rule star_index_V2:
    input:
        fasta = config["ref"]["reference2"],
    output:
        directory("/work/jawlab/kivanc/PeanutRnaSeq/StarIndexV2")
    threads:15
    priority:50
    params:
        extra = "",
        gtf = config["ref"]["annotation2"]
    log:
        "/work/jawlab/kivanc/PeanutRnaSeq/StarIndexV2/log/star_index_V2.log"
    wrapper:
        "0.49.0/bio/star/index"

rule starV2_pass1:
    input:
        fq1=GetClean(0),
        fq2=GetClean(1)
    output:
        "results/starV2/pass1/{smp}/Aligned.out.bam"
    log:
        "results/starV2/pass1/logs/{smp}.log"
    priority:1
    params:
        # path to STAR reference genome index
        index="/work/jawlab/kivanc/PeanutRnaSeq/StarIndexV2",
        extra="--outSAMtype BAM Unsorted --sjdbGTFfile {}".format(
              config["ref"]["annotation2"])
    threads:30
    wrapper:
        "0.49.0/bio/star/align"

rule get_junctions:
    input:
        expand("results/starV2/pass1/{smp}/SJ.out.tab", smp=sample_id)
    output:
        sj="results/starV2/junctions/SJ.filtered.tab"
    shell:
        """
        cat {input} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output}
        """

rule starV2_pass2:
    input:
        fq1=GetClean(0),
        fq2=GetClean(1),
        sj="results/starV2/junctions/SJ.filtered.tab"
    output:
        "results/starV2/pass2/{smp}/Aligned.out.bam",
        "results/starV2/pass2/{smp}/Aligned.toTranscriptome.out.bam",
        "results/starV2/pass2/{smp}/Aligned.sortedByCoord.out.bam",
        "results/starV2/pass2/{smp}/ReadsPerGene.out.tab"
    log:
        "results/starV2/pass2/logs/{smp}.log"
    params:
        # path to STAR reference genome index
        index="/work/jawlab/kivanc/PeanutRnaSeq/StarIndexV2",
        extra="--outSAMunmapped Within --outSAMtype BAM SortedByCoordinate Unsorted --quantMode GeneCounts TranscriptomeSAM --alignIntronMax 10000 --sjdbFileChrStartEnd {} --sjdbGTFfile {}".format(
              rules.get_junctions.output.sj, config["ref"]["annotation2"])
    threads:30
    wrapper:
        "0.49.0/bio/star/align"
    
