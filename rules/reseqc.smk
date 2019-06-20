rule gtf2bed:
    input:
        anno=ANNOTATION
    output:
        bed="results/rseqc/tifrunner_annotation.bed",
        db=temp("results/rseqc/annotation.db")
    log:
        "results/rseqc/logs/gtf2bed.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"

