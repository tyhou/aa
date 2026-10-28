#!/bin/bash
#SBATCH --job-name=trimmomatic_cut1
#SBATCH -p sonmi
#SBATCH -n 32
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=compute-0-3
##SBATCH -mem=200G
##SBATCH -t 72:00:00

ulimit -s unlimited
ulimit -l unlimited

cd $SLURM_SUBMIT_DIR

for i in $(cat SampleID_cut1)
do
	java -jar /share/home/hty/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
	rawdata/${i}_R1.fq.gz rawdata/${i}_R2.fq.gz \
	cleandata/${i}_trim_pair_R1.fq.gz cleandata/unpair/${i}_trim_unpair_R1.fq.gz \
	cleandata/${i}_trim_pair_R2.fq.gz cleandata/unpair/${i}_trim_unpair_R2.fq.gz \
	-threads 32 \
	-phred33 \
	ILLUMINACLIP:/share/home/hty/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
	LEADING:3 \
	TRAILING:3 \
	SLIDINGWINDOW:4:15 \
	MINLEN:60
echo $i
done
exit
