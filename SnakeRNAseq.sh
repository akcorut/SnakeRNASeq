#PBS -S /bin/bash
#PBS -q wallace_q
#PBS -N SnakemakeRnaSeq
#PBS -l walltime=240:00:00
#PBS -l nodes=1:ppn=24
#PBS -l mem=240gb

#PBS -M ac32082@uga.edu
#PBS -m ae
#PBS -j oe

cd /scratch/ac32082/PeanutRnaSeq

module load Anaconda3/5.0.1
source activate snakemake.515

export LC_ALL=en_SG.utf8
export LANG=en_SG.utf8

snakemake --use-conda --cores 24 --latency-wait 120 -s Snakefile