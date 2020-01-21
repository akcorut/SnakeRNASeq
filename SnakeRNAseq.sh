#PBS -S /bin/bash
#PBS -q wallace_q 
#PBS -N SnakemakeRnaSeq
<<<<<<< HEAD
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=40
#PBS -l mem=210gb
=======
#PBS -l walltime=480:00:00
#PBS -l nodes=1:ppn=40
#PBS -l mem=230gb
>>>>>>> 58c7e0000fcb1f754282d482021c355d4887289a

#PBS -M ac32082@uga.edu
#PBS -m ae
#PBS -j oe

cd /scratch/ac32082/PeanutRnaSeq

module load Anaconda3/5.0.1
source activate snakemake-rnaseq

export LC_ALL=en_SG.utf8
export LANG=en_SG.utf8

snakemake --use-conda --cores 40 -s Snakefile