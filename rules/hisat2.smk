rule extract_splice_sites:
    """Get splice sites from ref. annotation file"""
    input:
        gtf = "/work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf"
    output:
        splice_sites = SS_DIR + "/tifrunner_splice_sites.tsv"
    shell:
        "hisat2_extract_splice_sites.py {input.gtf} > {output.splice_sites}"

rule extract_exons:
    """Get exon sites from ref. annotation file"""
    input:
        gtf = "/work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf"
    output:
        exons = EXON_DIR + "/tifrunner_exons.tsv"
    shell:
        "hisat2_extract_exons.py {input.gtf} > {output.exons}"

rule hisat2_index:
    input:
        fasta = REFERENCE,
        splice_sites = SS_DIR + "/tifrunner_splice_sites.tsv",
        exons = EXON_DIR + "/tifrunner_exons.tsv"
    output:
        protected(expand(
            INDEX_DIR + "/tifrunner.{extension}.ht2",
            extension="1 2 3 4 5 6 7 8".split()))
    params:
        prefix = INDEX_DIR + "/tifrunner"
    log:
        INDEX_DIR + "/logs/hisat2_index.log"
    threads: 40
    shell:
        "hisat2-build "
            "-p {threads} "
            "--ss {input.splice_sites} "
            "--exon {input.exons} "
            "{input.fasta} "
            "{params.prefix} "
        "2> {log}"

rule hisat2:
    input:
        index= expand(
            INDEX_DIR + "/tifrunner.{extension}.ht2",
            extension="1 2 3 4 5 6 7 8".split()),
        r1= trim_data + "/{smp}_R1_val_1.fq.gz",
        r2= trim_data + "/{smp}_R2_val_2.fq.gz"
    output:
        r1 = temp("results/hisat2/{smp}.trimmed.sam")
    params:
        index= INDEX_DIR + "/tifrunner"
    threads: 32
    shell: 
        """
        hisat2 -p {threads} -x {params.index} --dta --rna-strandness RF -1 {input.r1} -2 {input.r2} -S {output.r1}
        """

rule sam2bam:
    input:
        r1 = "results/hisat2/{smp}.trimmed.sam"
    output:
        r1 = protected("results/hisat2/{smp}.trimmed.bam")
    threads:
        32
    shell:
        """
        samtools sort -@ {threads} -o {output.r1} {input.r1}
        """