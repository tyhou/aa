#!/bin/bash
#SBATCH --job-name=coverage_depth
#SBATCH -p sonmi
#SBATCH -n 32
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=compute-0-2
##SBATCH -mem=200G
##SBATCH -t 72:00:00

ulimit -s unlimited
ulimit -l unlimited

module load bwa-mem2
cd $SLURM_SUBMIT_DIR
mkdir coverage
for id in $(cat SampleID)
do
samtools mpileup markdup_ref_HLJ/${id}_markdup.bam |perl -alne '{$pos{$F[0]}++;$depth{$F[0]}+=$F[3]} END{print "$_\t$pos{$_}\t$depth{$_}" foreach sort keys %pos}' 1> coverage/coverage_depth_${id}.txt 2>coverage/coverage_depth_${id}.log
sed -i "s/$/\t${id}/" coverage/coverage_depth_${id}.txt
echo "$id in samtools mpileup"
done
cat coverage/coverage_depth_*.txt | grep -v 'mpileup' > all_coverage_depth.txt
samtools view -H markdup_ref_HLJ/11730_markdup.bam |awk '{if($0~"SQ")print$2,$3}'|sed 's/SN://g'|sed 's/ LN:/\t/g' > ChrLengthIIdA20G1.txt
awk 'NR==FNR{var[$1]=$2;next}NR>FNR{if($1 in var){print $0 "\t" var[$1]} else {print $0 "\t" "None"}}' ChrLengthIIdA20G1.txt all_coverage_depth.txt | awk '{print $0"\t"$2/$5*100"\t"$3/$5}' |sed "1i Chr\tcov_length\tcov_base\tSample\tChr_length\tcov\tdepth" > all_coverage_depth_final.txt
echo "coverage and depth were calculated"
exit
