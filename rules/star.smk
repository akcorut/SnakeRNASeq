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
        fq1= trim_data + "/{smp}_cutadapt_R1.fastq.gz",
        fq2= trim_data + "/{smp}_cutadapt_R2.fastq.gz
    output:
        "results/star/{smp}/Aligned.out.sam",
        "results/star/{smp}/ReadsPerGene.out.tab"
    log:
        "results/star/logs/{smp}.log"
    params:
        index= rules.star_index.output.directory,
        extra= "--quantMode GeneCounts --sjdbGTFfile {}".format(rules.gff3_to_gtf.output.gtf)
    threads: 40
    wrapper:
        "0.35.1/bio/star/align"
