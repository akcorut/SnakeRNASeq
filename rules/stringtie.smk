rule stringtie:
    input:
        r1 = "results/hisat2/{smp}.cutadapt.bam"
    output:
        r1 = "results/stringtie/{smp}/transcript.gtf",
        r2 = "results/stringtie/{smp}/gene_abundances.tsv",
        r3 = "results/stringtie/{smp}/cov_ref.gtf"
    params:
        gtf= rules.gff3_to_gtf.output.gtf

    shell: """
    stringtie -G {params.gtf} --rf -e -B -o {output.r1} -A {output.r2} -C {output.r3} --rf {input.r1}
"""