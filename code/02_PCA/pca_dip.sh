#!/bin/bash
#SBATCH --job-name=plink-pca
#SBATCH -p sonmi
#SBATCH -n 16
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=compute-0-1
##SBATCH -mem=200G
##SBATCH -t 72:00:00

ulimit -s unlimited
ulimit -l unlimited

cd $SLURM_SUBMIT_DIR

for i in F1SOSS;
do
plink --allow-extra-chr --vcf ${i}_dip_Marker_1065.vcf --noweb --make-bed --out ${i}_dip_Marker_1065
plink --allow-extra-chr -bfile ${i}_dip_Marker_1065 --pca 20 --out ${i}_dip_Marker_1065
done
exit
