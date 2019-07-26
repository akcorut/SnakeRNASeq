rule hisat2:
    input:
        index= expand(
            INDEX_DIR + "/tifrunner.{extension}.ht2",
            extension="1 2 3 4 5 6 7 8".split()),
        r1= trim_data + "/{smp}_cutadapt_R1.fastq.gz",
        r2= trim_data + "/{smp}_cutadapt_R2.fastq.gz"
    output:
        r1 = "results/hisat2/{smp}.cutadapt.sam"
    params:
        index= INDEX_DIR + "/tifrunner"
    threads: 40
    shell: 
        """
        hisat2 -p {threads} -x {params.index} --dta --rna-strandness RF -1 {input.r1} -2 {input.r2} -S {output.r1}
        """

rule sam2bam:
    input:
        r1 = "results/hisat2/{smp}.cutadapt.sam"
    output:
        r1 = "results/hisat2/{smp}.cutadapt.bam"
    threads:
        32
    shell:
        """
        samtools sort -@ {threads} -o {output.r1} {input.r1}
        """