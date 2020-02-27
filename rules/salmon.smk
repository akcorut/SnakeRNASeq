rule salmon_meta:
    input:
        ref= REFERENCE,
        tcp= TRANSCRIPTS
    output:
        gent= "results/salmon/decoy/gentrome.fa",
        decoy= "results/salmon/decoy/decoys.txt",
        bak="results/salmon/decoy/decoys.txt.bak"
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

rule salmon_index:
    input:
        gent= "results/salmon/decoy/gentrome.fa",
        decoy= "results/salmon/decoy/decoys.txt",
    output:
        directory("results/salmon/index")
    priority:1
    log:
        "results/salmon/logs/index.log"
    conda:
        "../envs/salmon.yaml"
    priority:50
    threads:24
    shell:
        """
        salmon index -p {threads} -t {input.gent} -d {input.decoy} -i {output} &> {log}
        """

if config["salmon_mode"]["mapping_mode"]:
    rule salmon_quant_mapping:
        input:
            r1=GetClean(0),
            r2=GetClean(1),
            index = "results/salmon/index"
        output:
            directory("results/salmon/quant/{smp}"),
            mappings="results/salmon/mappings/{smp}_salmon_mappings"
        log:
    		"results/salmon/logs/{smp}.salmon.log"
        conda:
            "../envs/salmon.yaml"
        priority:-1
        threads:24
        shell:
            """
            salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -o {output} --validateMappings --gcBias --seqBias --writeUnmappedNames --writeMappings={output.mappings} -p {threads} --numBootstraps 100
            """

if config["salmon_mode"]["alignment_mode"]:
    rule make_transcript:
        input:
            ref= REFERENCE,
            gtf= rules.gff3_to_gtf.output.gtf
        output:
            "results/salmon_align/transcript/arahy_transcripts.fa"
        conda:
            "../envs/gffread.yaml"
        shell:
            '''
            gffread -w {output} -g {input.ref} {input.gtf}
            '''

    rule salmon_quant_alignment:
        input:
            bam="results/star/{smp}/Aligned.toTranscriptome.out.bam",
            tcp = "results/salmon_align/transcript/arahy_transcripts.fa",
            gtf= rules.gff3_to_gtf.output.gtf,
        output:
            directory("results/salmon_align/quant/{smp}_salmon_quant_align")
        priority:-1
        log:
    		"results/salmon_align/logs/{smp}.salmon_align.log"
        conda:
            "../envs/salmon.yaml"
        threads:24
        shell:
            '''
            salmon quant -t {input.tcp} -l A -a {input.bam} -o {output} --gcBias --seqBias --writeUnmappedNames -p {threads} -g {input.gtf} --numBootstraps 100
            '''