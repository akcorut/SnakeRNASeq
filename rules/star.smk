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

def GetReads(num):
    raw1="/work/jawlab/kivanc/PeanutRnaSeq/data/fastq/{smp}_R1.fastq.gz"
    raw2="/work/jawlab/kivanc/PeanutRnaSeq/data/fastq/{smp}_R2.fastq.gz"
    trim1="results/trimmed/{smp}_R1_val_1.fq.gz"
    trim2="results/trimmed/{smp}_R2_val_2.fq.gz"
    if config["trimming"]["skip"]:
        sample_list = [raw1,raw2]
        return sample_list[num];
    else:
        sample_list = [trim1, trim2]
        return sample_list[num];

rule star_align:
    input:
        fq1=GetReads(0),
        fq2=GetReads(1)
    output:
        "results/star/{smp}/Aligned.out.bam",
        "results/star/{smp}/ReadsPerGene.out.tab",
        "results/star/{smp}/Aligned.toTranscriptome.out.bam"
    log:
        "results/star/logs/{smp}.log"
    params:
        index= "/work/jawlab/kivanc/PeanutRnaSeq/reference/star_index",
        extra= "--twopassMode Basic --outSAMunmapped Within --limitOutSJcollapsed 1000000 --limitSjdbInsertNsj 1000000 --outFilterMultimapNmax 100 --outFilterMismatchNmax 33 --outFilterMismatchNoverLmax 0.3 --seedSearchStartLmax 12 --alignSJoverhangMin 15 --alignEndsType Local --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0.3 --winAnchorMultimapNmax 50 --alignSJDBoverhangMin 3 --quantMode GeneCounts TranscriptomeSAM --outSAMtype BAM Unsorted --sjdbOverhang 100 --sjdbGTFfile {}".format(rules.gff3_to_gtf.output.gtf)
    threads: 16
    wrapper:
        "0.35.1/bio/star/align"
