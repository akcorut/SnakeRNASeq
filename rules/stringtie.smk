rule stringtie:
    input:
        r1 = "results/hisat2/{smp}.trimmed.bam"
    output:
        r1 = "results/stringtie/{smp}/transcript.gtf",
        r2 = "results/stringtie/{smp}/gene_abundances.tsv",
        r3 = "results/stringtie/{smp}/cov_ref.gtf"
    conda:
        "../envs/stringtie.yaml"
    params:
        gff = ANNOTATION 
    threads: 16
    shell:
        """
        stringtie -p {threads} -G {params.gff} --rf -e -B -o {output.r1} -A {output.r2} -C {output.r3} --rf {input.r1}
        """

rule gtf_merge:
    input:
        tcp = expand("results/stringtie/{smp}/transcript.gtf", smp=sample_id)
    output:
        "results/stringtie/merge_transcripts.gtf"
    params:
        gtf= rules.gff3_to_gtf.output.gtf,
    threads: 16
    shell: 
        """
        stringtie -p {threads} --merge  -G {params.gtf} -o {output} {input.tcp}
        """

rule gff_compare:
    input:
        "results/stringtie/merge_transcripts.gtf"
    output:
        "results/stringtie/gffcompare/stringtie.stats"
    conda:
        "../envs/gffcompare.yaml"
    params:
        prefix= "results/stringtie/gffcompare/stringtie",
        gtf= rules.gff3_to_gtf.output.gtf
    threads: 16
    shell:
        """
        gffcompare -V -G -r {params.gtf} -o {params.prefix} {input}
        """ 
