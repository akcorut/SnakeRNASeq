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