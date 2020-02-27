def GetClean(num):
    clean1="results/bbsplit/{smp}_R1_clean.fq"
    clean2="results/bbsplit/{smp}_R2_clean.fq"
    trim1="results/trimmed/{smp}_R1_val_1.fq.gz"
    trim2="results/trimmed/{smp}_R2_val_2.fq.gz"
    if config["decontamination"]["skip"]:
        sample_list = [trim1,trim2]
        return sample_list[num];
    else:
        sample_list = [clean1, clean2]
        return sample_list[num]

def GetReads(num):
    raw1="/work/jawlab/kivanc/PeanutRnaSeq/data/fastq/{smp}_R1.fastq.gz"
    raw2="/work/jawlab/kivanc/PeanutRnaSeq/data/fastq/{smp}_R2.fastq.gz"
    trim1="results/trimmed/{smp}_R1_val_1.fq.gz"
    trim2="results/trimmed/{smp}_R2_val_2.fq.gz"
    if config["trimming"]["skip"]:
        sample_list = [raw1,raw2]
        return sample_list[num];
    else:
        sample_list = [trim1, trim2]
        return sample_list[num];