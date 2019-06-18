rule cutadapt:
    input:
        r1= raw_data + "/{smp}_R1.fastq.gz",
        r2= raw_data + "/{smp}_R2.fastq.gz"
    output:
        r1="/work/jawlab/kivanc/PeanutRnaSeq/data/cutadapt/{smp}_cutadapt_R1.fastq.gz",
        r2="/work/jawlab/kivanc/PeanutRnaSeq/data/cutadapt/{smp}_cutadapt_R2.fastq.gz",
    log:
        "/work/jawlab/kivanc/PeanutRnaSeq/data/cutadapt/logs/{smp}.cutadapt.log"
    threads:40
    shell:
        """
        cutadapt -j {threads} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o {output.r1} -p {output.r2} {input.r1} {input.r2} &> {log}
        """