#!/bin/bash
#SBATCH --job-name=mapping_ref
#SBATCH -p sonmi
#SBATCH -n 32
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=compute-0-0
##SBATCH -mem=200G
##SBATCH -t 72:00:00

ulimit -s unlimited
ulimit -l unlimited

module load bwa-mem2
cd $SLURM_SUBMIT_DIR

for id in $(cat refID)
do
bwa-mem2 mem -t 32 /biodata/hty/reference/index_BWA_Cpar_IId_11730_Harbin_A20G1_TGS/11730 \
cleandata/${id}_trim_pair_R1.fq.gz cleandata/${id}_trim_pair_R2.fq.gz | samtools sort -O BAM -@ 32 -o BWA_ref_HLJ/${id}_sort.bam > BWA_ref_HLJ/${id}.log
echo "$id in bwa-mem2"
done
echo "bwa-mem2 is done"

for id in $(cat refID)
do
sambamba markdup -r -t 32 BWA_ref_HLJ/${id}_sort.bam markdup_ref_HLJ/${id}_markdup.bam
samtools flagstat markdup_ref_HLJ/${id}_markdup.bam -@ 32 > markdup_ref_HLJ/${id}_stat.log
echo "$id in sambamba"
done
echo "sambamba is done"
exit
