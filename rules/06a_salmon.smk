rule salmon_meta:
    input:
        ref= REFERENCE,
        tcp= TRANSCRIPTS
    output:
        gent= "results/06_alignment_free/06a_salmon/decoy/gentrome.fa",
        decoy= "results/06_alignment_free/06a_salmon/decoy/decoys.txt",
        bak="results/06_alignment_free/06a_salmon/decoy/decoys.txt.bak"
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
        gent= "results/06_alignment_free/06a_salmon/decoy/gentrome.fa",
        decoy= "results/06_alignment_free/06a_salmon/decoy/decoys.txt",
    output:
        directory("results/06_alignment_free/06a_salmon/index")
    priority:1
    log:
        "results/06_alignment_free/06a_salmon/logs/index.log"
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
            index = "results/06_alignment_free/06a_salmon/index"
        output:
            directory("results/06_alignment_free/06a_salmon/quant/{smp}"),
            mappings="results/06_alignment_free/06a_salmon/mappings/{smp}_salmon_mappings"
        log:
    		"results/06_alignment_free/06a_salmon/logs/{smp}.salmon.log"
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
            bam="results/04_alignment/04a_alignment_results/star/pass2/{smp}/Aligned.toTranscriptome.out.bam",
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

rule multiqc_salmon:
    input:
        expand("results/06_alignment_free/06a_salmon/quant/{smp}", smp=sample_id)
    output:
        "results/06_alignment_free/06a_salmon/salmon_multiqc.html"
    log:
        "results/06_alignment_free/06a_salmon/logs/multiqc.log"
    wrapper:
        "0.49.0/bio/multiqc"

if config["salmon_mode"]["alignment_mode"]:
    rule multiqc_salmon_align:
        input:
            expand("results/salmon_align/quant/{smp}_salmon_quant_align", smp=sample_id)
        output:
            "results/salmon_align/salmon_align_multiqc.html"
        log:
            "results/salmon_align/logs/multiqc.log"
        wrapper:
            "0.47.0/bio/multiqc"