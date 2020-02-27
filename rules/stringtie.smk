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
