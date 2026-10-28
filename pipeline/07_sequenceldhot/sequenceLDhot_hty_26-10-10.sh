#!/bin/bash
#SBATCH --job-name=sequenceLDhot
#SBATCH -p sonmi
#SBATCH -n 56
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=compute-0-1
##SBATCH -mem=200G
##SBATCH -t 72:00:00

ulimit -s unlimited
ulimit -l unlimited

BASE_DIR=/share/home/hty/work/IId-cross-Single-Oocyst-26-06-03/phase-sequenceLD-dip

mkdir -p output_W20_Q75_26-10-10


cd $SLURM_SUBMIT_DIR
for i in {1..8};
do
./sequenceLDhot in1 \
	-P ${BASE_DIR}/input_sequenceLDhot/chr${i}.phase.o \
	-R ${BASE_DIR}/input_sequenceLDhot/chr${i}.phase.o_recom \
	-W 20 -Q 75 output_W20_Q75_26-10-10/chr${i}_sequenceLDhot

done
echo "sequenceLDhot finished"
exit

