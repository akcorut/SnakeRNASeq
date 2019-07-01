rule star_index:
    input:
        fasta = REFERENCE
    output:
        directory("/work/jawlab/kivanc/PeanutRnaSeq/reference/star_index")
    threads:
        40
    params:
        extra = "",
        gtf = rules.gff3_to_gtf.output.gtf
    log:
        "/work/jawlab/kivanc/PeanutRnaSeq/reference/star_index/logs/star_index.log"
    wrapper:
        "0.35.1/bio/star/index"

rule star_align:
    input:
        fq1= raw_data + "/{smp}_R1.fastq.gz",
        fq2= raw_data + "/{smp}_R2.fastq.gz"
    output:
        "results/star/{smp}/Aligned.out.bam",
        "results/star/{smp}/ReadsPerGene.out.tab",
        "results/star/{smp}/Aligned.toTranscriptome.out.bam",
        "results/star/{smp}/Aligned.sortedByCoord.out.bam"

    log:
        "results/star/logs/{smp}.log"
    params:
        index= "/work/jawlab/kivanc/PeanutRnaSeq/reference/star_index",
        extra= "--quantMode GeneCounts TranscriptomeSAM --outSAMtype BAM SortedByCoordinate Unsorted --sjdbOverhang 100 --sjdbGTFfile {}".format(rules.gff3_to_gtf.output.gtf)
    threads: 40
    wrapper:
        "0.35.1/bio/star/align"
