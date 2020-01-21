"""
Author: Kivanc Corut
A Snakemake workflow to process paired end RNA-Seq data.
"""
import glob
import os
import pandas as pd
from snakemake.io import expand
from snakemake.utils import R
from snakemake.io import glob_wildcards
import re
from os.path import join

configfile: "config.yaml"

raw_data = "/work/jawlab/kivanc/PeanutRnaSeq/data/fastq"
trim_data = "/scratch/ac32082/PeanutRnaSeq/results/trimmed"

REFERENCE = config["ref"]["reference"]
ANNOTATION = config["ref"]["annotation"]
INDEX_DIR = config["ref"]["index"]
EXON_DIR = config["ref"]["exon"]
SS_DIR = config["ref"]["ss"]
TRANSCRIPTOME = config["ref"]["tcp"]
TRANSCRIPTS = config["ref"]["transcript"]

sample_id, run_id = glob_wildcards(raw_data + "/{smp}_R{run}.fastq.gz")

rule all:
    input:
        expand("results/FastQC/{smp}_R1_fastqc.html", smp=sample_id),
        expand("results/FastQC/{smp}_R2_fastqc.html", smp=sample_id),
        "results/MultiQC/fastq_multiqc.html",
        "/work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf",
        SS_DIR + "/tifrunner_splice_sites.tsv",
        EXON_DIR + "/tifrunner_exons.tsv",
        expand(
            INDEX_DIR + "/tifrunner.{extension}.ht2",
            extension="1 2 3 4 5 6 7 8".split()),
        "/work/jawlab/kivanc/PeanutRnaSeq/reference/star_index",
        expand("results/star/{smp}/Aligned.out.bam", smp=sample_id),
        expand("results/star/{smp}/ReadsPerGene.out.tab", smp=sample_id),
        expand("results/star/{smp}/Aligned.toTranscriptome.out.bam", smp=sample_id),
        expand("results/star/{smp}/Aligned.sortedByCoord.out.bam", smp=sample_id),
        expand("results/trimmed/{smp}_R1_val_1.fq.gz", smp=sample_id),
        expand("results/trimmed/{smp}_R2_val_2.fq.gz", smp=sample_id),
        "results/trimmed/trim_galore_multiqc_report.html",
        "results/MultiQCTrim/multiqc_report_trim_galore.html",
        expand("results/hisat2/{smp}.trimmed.bam", smp=sample_id),
        "results/rseqc/multiqc_report.html",
        "results/feature/featureCounts_results.txt",
        expand("results/stringtie/{smp}/transcript.gtf", smp=sample_id),
        expand("results/stringtie/{smp}/gene_abundances.tsv", smp=sample_id),
        expand("results/stringtie/{smp}/cov_ref.gtf", smp=sample_id),
        "results/stringtie/merge_transcripts.gtf",
        "results/stringtie/gffcompare/stringtie.stats",
        "results/feature/feature_multiqc.html",
        expand("results/kallisto/quant/quant_results_{smp}", smp=sample_id),
        "results/kallisto/kallisto_multiqc.html",
        expand("results/salmon/quant/{smp}_salmon_quant", smp=sample_id),
        "results/salmon/salmon_multiqc.html",
        config["rsem"]["rsemindex"] + ".n2g.idx.fa",
        expand("results/rsem/genes/{smp}.genes.results", smp=sample_id),
        "results/rsem/rsem_multiqc.html",
        expand("results/qualimap/{smp}/qualimapReport.html", smp=sample_id),
        "results/qualimap/qualimap_multiqc.html",
        expand("results/rrnaCheck/bam/{smp}_rrna.bam", smp=sample_id),
        "results/rrnaCheck/rrna_multiqc_report.html"

include: "rules/fastqc.smk"
include: "rules/multiqc.smk"
include: "rules/index.smk"
include: "rules/hisat2.smk"
include: "rules/star.smk"
include: "rules/trim_galore.smk"
include: "rules/rseqc.smk"
include: "rules/featureCounts.smk"
include: "rules/stringtie.smk"
include: "rules/kallisto.smk"
include: "rules/salmon.smk"
include: "rules/rsem.smk"
include: "rules/qualimap.smk"
include: "rules/rrna_check.smk"
