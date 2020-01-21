rule salmon_index_k21:
    input:
        gentrome= "results/salmon_k21/decoy/gentrome_gnm1.fa",
        decoy= "results/salmon_k21/decoy/decoys_gnm1.txt"
    output:
        directory("results/salmon_k21/index_k21")
    priority:1
    log:
        "results/salmon_k21/logs/index.log"
    conda:
        "../envs/salmon.yaml"
    priority:50
    threads:20
    shell:
        """
        salmon index -t {input.gentrome} -d {input.decoy} -p {threads} -i {output} &> {log} -k 21
        """

rule salmon_quant_k21:
    input:
        r1="results/trimmed/{smp}_R1_val_1.fq.gz",
        r2="results/trimmed/{smp}_R2_val_2.fq.gz",
        index = "results/salmon_k21/index_k21"
    output:
        directory("results/salmon_k21/quant/{smp}_salmon_quant")
    log:
		"results/salmon_k21/logs/{smp}.salmon.log"
    conda:
        "../envs/salmon.yaml"
    threads:20
    shell:
        """
        salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -o {output} --validateMappings --gcBias --writeUnmappedNames -p {threads} --numBootstraps 100
        """
    