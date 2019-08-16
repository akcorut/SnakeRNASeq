rule trinity:
    input:
        left=expand("results/trimmed/{smp}_R1_val_1.fq.gz", smp=sample_id),
        right=expand("results/trimmed/{smp}_R2_val_2.fq.gz", smp=sample_id)
    output:
        "results/trinity/Trinity.fasta"
    log:
        'results/trinity/logs/trinity.log'
    params:
        extra="--SS_lib_type RF",
        max_memory="900G"
    threads: 24
    wrapper:
        "0.36.0/bio/trinity"