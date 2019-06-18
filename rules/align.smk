rule hisat2:
    input:
        index_prefix= INDEX_DIR + "/tifrunner",
        index= expand(
            INDEX_DIR + "/tifrunner.{extension}.ht2",
            extension="1 2 3 4 5 6 7 8".split()),
        r1= trim_data + "/{smp}_cutadapt_R1.fastq.gz",
        r2= trim_data + "/{smp}_cutadapt_R2.fastq.gz"
    output:
        r1 = "results/hisat2/{smp}.cutadapt.sam"
    threads: 40
    shell: 
        """
        hisat2  -p {threads} -x {input.index_prefix} --dta  --rna-strandness RF -1 {input.r1}  -2 {input.r2} -S {output.r1}
        """