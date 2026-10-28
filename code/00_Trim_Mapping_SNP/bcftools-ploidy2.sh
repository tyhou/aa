#!/bin/bash
#SBATCH --job-name=bcftools-ploidy2
#SBATCH -p sonmi
#SBATCH -n 56
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=compute-0-0
##SBATCH -mem=200G
##SBATCH -t 72:00:00

ulimit -s unlimited
ulimit -l unlimited

cd $SLURM_SUBMIT_DIR
##二倍体call
bcftools mpileup -Ou -a FORMAT/AD,FORMAT/DP -f /biodata/hty/reference/Cpar/IIdA20G1-Harbin-11730-TGS-24-10-01/11730_spades_r.fasta \
	-b Sample_path -d 500 --threads 56 | bcftools call -mv -O z \
	-o ../vcf_dip/all_dip.vcf.gz --threads 56

echo "bcftools ploidy2 finished"
exit

