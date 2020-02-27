"""
Author: Kivanc Corut
A Snakemake workflow to process paired end RNA-Seq data.
"""
import glob
import os
import json
import pandas as pd
from snakemake.io import expand
from snakemake.utils import R
from snakemake.io import glob_wildcards
import re
from os.path import join, basename, dirname

configfile: "config.yaml"

raw_data = "/work/jawlab/kivanc/PeanutRnaSeq/data/fastq"
trim_data = "/scratch/ac32082/PeanutRnaSeq/results/trimmed"

REFERENCE = config["ref"]["reference"]
ANNOTATION = config["ref"]["annotation"]
INDEX_DIR = config["ref"]["index"]
EXON_DIR = config["ref"]["exon"]
SS_DIR = config["ref"]["ss"]
TRANSCRIPTS = config["ref"]["transcript"]

sample_id, run_id = glob_wildcards(raw_data + "/{smp}_R{run}.fastq.gz")

meristem_samples = pd.read_csv(config["samples"]["meristem"]).set_index("sample_meristem", drop=False)
leaf_samples = pd.read_csv(config["samples"]["leaf"]).set_index("sample_leaf", drop=False)

ALL_TARGET = []

if not config["trimming"]["skip"]:
    TRIMMED_R1= expand("results/trimmed/{smp}_R1_val_1.fq.gz", smp=sample_id)
    TRIMMED_R2= expand("results/trimmed/{smp}_R2_val_2.fq.gz", smp=sample_id)
    TRIMMED_MQC= expand("results/trimmed/trim_galore_multiqc_report.html")
    ALL_TARGET.extend(TRIMMED_R1)
    ALL_TARGET.extend(TRIMMED_R2)
    ALL_TARGET.extend(TRIMMED_MQC)

if not config["decontamination"]["skip"]:
    CHECK_RRNA= expand("results/rrnaCheck/rrna_multiqc_report.html")
    CLEAN_R1= expand("results/bbsplit/{smp}_R1_clean.fq", smp=sample_id)
    CLEAN_R2= expand("results/bbsplit/{smp}_R2_clean.fq", smp=sample_id)
    CHECK_CLEAN= expand("results/rrna_decont/rrna_clean_multiqc_report.html")
    CLEAN_MQC= expand("results/MultiQCDecont/multiqc_report_decont.html")
    ALL_TARGET.extend(CHECK_RRNA)
    ALL_TARGET.extend(CLEAN_R1)
    ALL_TARGET.extend(CLEAN_R2)
    ALL_TARGET.extend(CHECK_CLEAN)
    ALL_TARGET.extend(CLEAN_MQC)

if config["alignment"]["activate"]:
    
    if config["star"]["activate"]:
        GTF= expand("/work/jawlab/kivanc/PeanutRnaSeq/reference/tifrunner_gene_models.gtf")
        STAR_INDEX= expand("/work/jawlab/kivanc/PeanutRnaSeq/StarIndex")
        STAR_PASS1= expand("results/star/pass1/{smp}/Aligned.out.bam", smp=sample_id)
        STAR_PASS1_SJ= expand("results/star/pass1/{smp}/SJ.out.tab", smp=sample_id)
        SPLICE_JUNC= expand("results/star/junctions/SJ.filtered.tab")
        STAR_BAM= expand("results/star/pass2/{smp}/Aligned.out.bam", smp=sample_id)
        STAR_TCP_BAM= expand("results/star/pass2/{smp}/Aligned.toTranscriptome.out.bam", smp=sample_id)
        STAR_SORTED_BAM= expand("results/star/pass2/{smp}/Aligned.sortedByCoord.out.bam", smp=sample_id)
        STAR_GENE_COUNT= expand("results/star/pass2/{smp}/ReadsPerGene.out.tab", smp=sample_id)
        STAR_RSEQC= expand("results/rseqc/multiqc_report.html")
        QUALIMAP_REPORT= expand("results/qualimap/{smp}/qualimapReport.html", smp=sample_id)
        QUALIMAP_MQC= expand("results/qualimap/qualimap_multiqc.html")
        FCOUNTS_RESULTS= expand("results/feature/featureCounts_results.txt")
        FCOUNTS_MQC= expand("results/feature/feature_multiqc.html")
        #config["rsem"]["rsemindex"] + ".n2g.idx.fa"
        #expand("results/rsem/genes/{smp}.genes.results", smp=sample_id)
        #"results/rsem/rsem_multiqc.html"
        ALL_TARGET.extend(GTF)
        ALL_TARGET.extend(STAR_INDEX)
        ALL_TARGET.extend(STAR_PASS1)
        ALL_TARGET.extend(SPLICE_JUNC)
        ALL_TARGET.extend(STAR_BAM)
        ALL_TARGET.extend(STAR_TCP_BAM)
        ALL_TARGET.extend(STAR_SORTED_BAM)
        ALL_TARGET.extend(STAR_GENE_COUNT)
        ALL_TARGET.extend(STAR_RSEQC)
        ALL_TARGET.extend(QUALIMAP_REPORT)
        ALL_TARGET.extend(QUALIMAP_MQC)
        ALL_TARGET.extend(FCOUNTS_RESULTS)
        ALL_TARGET.extend(FCOUNTS_MQC)
    
    if config["hisat2"]["activate"]:
        SPLICE_SITES= expand("/work/jawlab/kivanc/PeanutRnaSeq/reference/splice")
        EXON= expand("/work/jawlab/kivanc/PeanutRnaSeq/reference/exon")
        HISAT_INDEX= expand(
            INDEX_DIR + "/tifrunner.{extension}.ht2",
            extension="1 2 3 4 5 6 7 8".split())
        HISAT_BAM= expand("results/hisat2/{smp}.bam", smp=sample_id)
        STRINGTIE_GTF= expand("results/stringtie/{smp}/transcript.gtf", smp=sample_id)
        STRINGTIE_ABUN= expand("results/stringtie/{smp}/gene_abundances.tsv", smp=sample_id)
        STRINGTIE_COV= expand("results/stringtie/{smp}/cov_ref.gtf", smp=sample_id)
        ALL_TARGET.extend(SPLICE_SITES)
        ALL_TARGET.extend(EXON)
        ALL_TARGET.extend(HISAT_INDEX)
        ALL_TARGET.extend(HISAT_BAM)
        ALL_TARGET.extend(STRINGTIE_GTF)
        ALL_TARGET.extend(STRINGTIE_ABUN)
        ALL_TARGET.extend(STRINGTIE_COV)
        

if config["alignment_free"]["activate"]:

    if config["salmon"]["activate"]:

        if config["salmon_mode"]["mapping_mode"]: 
            SALMON_QUASI= expand("results/salmon/quant/{smp}", smp=sample_id)
            SALMON__QUASI_MQC= expand("results/salmon/salmon_multiqc.html")
            ALL_TARGET.extend(SALMON_QUASI)
            ALL_TARGET.extend(SALMON__QUASI_MQC)

        if config["salmon_mode"]["alignment_mode"]:
            SALMON_ALIGN= expand("results/salmon_align/quant/{smp}_salmon_quant_align", smp=sample_id)
            SALMON_ALIGN_MQC=  expand("results/salmon_align/salmon_align_multiqc.html")
            ALL_TARGET.extend(SALMON_ALIGN)
            ALL_TARGET.extend(SALMON_ALIGN_MQC)

    if config["kallisto"]["activate"]:
        KALLISTO_QUANT= expand("results/kallisto/quant/{smp}", smp=sample_id)
        KALLISTO_MQC= expand("results/kallisto/kallisto_multiqc.html")
        ALL_TARGET.extend(KALLISTO_QUANT)
        ALL_TARGET.extend(KALLISTO_MQC)

INIT_QC_R1= expand("results/FastQC/{smp}_R1_fastqc.html", smp=sample_id)
INIT_QC_R2= expand("results/FastQC/{smp}_R2_fastqc.html", smp=sample_id)
INIT_MQC= expand("results/MultiQC/fastq_multiqc.html")
ALL_TARGET.extend(INIT_QC_R1)
ALL_TARGET.extend(INIT_QC_R2)
ALL_TARGET.extend(INIT_MQC)

rule all:
    input: ALL_TARGET

include: "rules/common.smk"
include: "rules/fastqc.smk"
include: "rules/multiqc.smk"
include: "rules/trim_galore.smk"
include: "rules/rrna_check.smk"
include: "rules/rrna_remove.smk"
include: "rules/hisat2.smk"
include: "rules/stringtie.smk"
include: "rules/star.smk"
#include: "rules/starV2.smk"
include: "rules/rseqc.smk"
include: "rules/qualimap.smk"
include: "rules/featureCounts.smk"
#include: "rules/rsem.smk"
include: "rules/kallisto.smk"
include: "rules/salmon.smk"
#include: "rules/salmonV2.smk"
#include: "rules/kallistoV2.smk"
