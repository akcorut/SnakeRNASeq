rule salmon_decoy:
    input:
        ref= REFERENCE,
        gtf= rules.gff3_to_gtf.output.gtf,
        tcp= TRANSCRIPTOME
    output:
        directory("results/salmon/decoy")
    conda:
        "../envs/salmon.yaml"
    threads:8
    shell:
        """
        bash scripts/SalmonTools/scripts/generateDecoyTranscriptome.sh -j {threads} -g {input.ref} -t {input.tcp} -a {input.gtf} -o {output}
        """

rule salmon_index:
    input:
        fasta= "results/salmon/decoy/gentrome.fa",
        decoy= "results/salmon/decoy/decoys.txt"
    output:
        directory("results/salmon/index")
    log:
        "results/salmon/logs/index.log"
    conda:
        "../envs/salmon.yaml"
    priority:50
    threads:20
    shell:
        """
        salmon index -p {threads} -t {input.fasta} -i {output} --decoys {input.decoy} -k 31 &> {log}
        """

rule salmon_quant:
    input:
        r1="results/trimmed/{smp}_R1_val_1.fq.gz",
        r2="results/trimmed/{smp}_R2_val_2.fq.gz",
        index = "results/salmon/index"
    output:
        directory("results/salmon/quant/{smp}_salmon_quant")
    log:
		"results/salmon/logs/{smp}.salmon.log"
    conda:
        "../envs/salmon.yaml"
    threads:20
    shell:
        """
        salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -o {output} --validateMappings --gcBias --seqBias --posBias --writeUnmappedNames -p {threads} --numBootstraps 100
        """
    