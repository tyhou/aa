for i in F1SOSS SOLine;
do
vcftools --vcf ${i}_dip_Marker_1065.vcf --plink --out ${i}_dip_Marker_1065
plink --vcf ${i}_dip_Marker_1065.vcf --make-bed --allow-extra-chr --out ${i}_dip_Marker_1065
plink -bfile ${i}_dip_Marker_1065 --het --allow-extra-chr --out het/sample_${i}
plink --file ${i}_dip_Marker_1065 --hardy --out het/site_${i} &> /dev/zero
done