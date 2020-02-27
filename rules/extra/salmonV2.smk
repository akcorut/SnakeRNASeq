rule salmon_meta_V2:
    input:
        ref= config["ref"]["reference2"],
        tcp= config["ref"]["transcript2"]
    output:
        gent= "results/salmonV2/decoy/gentromeV2.fa",
        decoy= "results/salmonV2/decoy/decoysV2.txt",
        bak= "results/salmonV2/decoy/decoysV2.txt.bak"
    priority:50
    conda:
        "../envs/salmon.yaml"
    threads:8
    shell:
        """
        grep "^>" {input.ref} | cut -d " " -f 1 > {output.decoy}
        sed -i.bak -e 's/>//g' {output.decoy}
        cat {input.tcp} {input.ref} > {output.gent}
        """

rule salmon_index_V2:
    input:
        gent= "results/salmonV2/decoy/gentromeV2.fa",
        decoy= "results/salmonV2/decoy/decoysV2.txt",
    output:
        directory("results/salmonV2/index")
    priority:1
    log:
        "results/salmonV2/logs/index.log"
    conda:
        "../envs/salmon.yaml"
    priority:50
    threads:20
    shell:
        """
        salmon index -p {threads} -t {input.gent} -d {input.decoy} -i {output} &> {log}
        """

if config["salmon_mode"]["mapping_mode"]:
    rule salmon_quant_mapping_V2:
        input:
            r1="results/bbsplit/{smp}_R1_clean.fq",
            r2="results/bbsplit/{smp}_R2_clean.fq",
            index = "results/salmonV2/index"
        output:
            directory("results/salmonV2/quant/{smp}_salmon_quant"),
            mappings="results/salmonV2/mappings/{smp}_salmon_mappings"
        log:
    		"results/salmonV2/logs/{smp}.salmon.log"
        conda:
            "../envs/salmon.yaml"
        priority:-1
        threads:20
        shell:
            """
            salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -o {output} --validateMappings --gcBias --seqBias --writeUnmappedNames --writeMappings={output.mappings} -p {threads} --numBootstraps 100
            """