rule salmon_new_index:
    input:
        gentrome= "results/salmon_new/decoy/gentrome.fa",
        decoy= "results/salmon_new/decoy/decoys.txt"
    output:
        directory("results/salmon_new/index")
    priority:1
    log:
        "results/salmon_new/logs/index.log"
    conda:
        "../envs/salmon.yaml"
    priority:50
    threads:20
    shell:
        """
        salmon index -t {input.gentrome} -d {input.decoy} -p {threads} -i {output} &> {log}
        """

rule salmon_new_quant:
    input:
        r1="results/trimmed/{smp}_R1_val_1.fq.gz",
        r2="results/trimmed/{smp}_R2_val_2.fq.gz",
        index = "results/salmon_new/index"
    output:
        directory("results/salmon_new/quant/{smp}_salmon_quant")
    log:
		"results/salmon_new/logs/{smp}.salmon.log"
    conda:
        "../envs/salmon.yaml"
    threads:20
    shell:
        """
        salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -o {output} --validateMappings --gcBias --writeUnmappedNames -p {threads} --numBootstraps 100
        """
    