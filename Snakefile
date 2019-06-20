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
trim_data = "/work/jawlab/kivanc/PeanutRnaSeq/data/cutadapt"

REFERENCE = config["ref"]["reference"]
ANNOTATION = config["ref"]["annotation"]
INDEX_DIR = config["ref"]["index"]
EXON_DIR = config["ref"]["exon"]
SS_DIR = config["ref"]["ss"]


sample_id, run_id = glob_wildcards(raw_data + "/{smp}_R{run}.fastq.gz")

rule all:
    input:
        expand("results/FastQC/{smp}_R1_fastqc.html", smp=sample_id),
        expand("results/FastQC/{smp}_R2_fastqc.html", smp=sample_id),
        "results/MultiQC/fastq_multiqc.html",
        expand("/work/jawlab/kivanc/PeanutRnaSeq/data/cutadapt/{smp}_cutadapt_R1.fastq.gz", smp=sample_id),
        expand("/work/jawlab/kivanc/PeanutRnaSeq/data/cutadapt/{smp}_cutadapt_R2.fastq.gz", smp=sample_id),
        expand("results/FastQCCut/{smp}_cutadapt_R1_fastqc.html", smp=sample_id),
        expand("results/FastQCCut/{smp}_cutadapt_R2_fastqc.html", smp=sample_id),
        "results/MultiQCCut/fastq_cutadapt_multiqc.html",
        "/work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf",
        SS_DIR + "/tifrunner_splice_sites.tsv",
        EXON_DIR + "/tifrunner_exons.tsv",
        expand(
            INDEX_DIR + "/tifrunner.{extension}.ht2",
            extension="1 2 3 4 5 6 7 8".split()),
        expand("results/hisat2/{smp}.cutadapt.sam", smp=sample_id),
        expand("results/hisat2/{smp}.cutadapt.bam", smp=sample_id


include: "rules/fastqc.smk"
include: "rules/multiqc.smk"
include: "rules/cutadapt.smk"
include: "rules/index.smk"
include: "rules/align.smk"