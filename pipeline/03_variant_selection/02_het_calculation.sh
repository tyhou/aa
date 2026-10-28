#!/bin/bash

# 创建输出目录
mkdir -p het

# 样本名称（根据您的文件命名）
samples=("F1SOSS" "SOLine")
# SNP 集名称
snp_sets=("Marker_1065" "overlap_1391")

for sample in "${samples[@]}"; do
    for snp_set in "${snp_sets[@]}"; do
        vcf="${sample}_dip_${snp_set}.vcf"
        
        if [[ ! -f "$vcf" ]]; then
            echo "Warning: $vcf not found, skipping..."
            continue
        fi
        
        prefix="${sample}_dip_${snp_set}"
        echo "Processing $vcf ..."
        
        # 1. 生成 Plink 文本格式
        vcftools --vcf "$vcf" --plink --out "$prefix"
        
        # 2. 生成 Plink 二进制格式
        plink --vcf "$vcf" --make-bed --allow-extra-chr --out "$prefix"
        
        # 3. 计算个体杂合度（产生 .het 文件）
        plink -bfile "$prefix" --het --allow-extra-chr --out "het/sample_${prefix}"
        
        # 4. 计算位点 Hardy-Weinberg 平衡
        plink --file "$prefix" --hardy --allow-extra-chr --out "het/site_${prefix}" &> /dev/null
        
        echo "Finished $vcf"
    done
done

# 创建汇总文件（含表头）
echo -e "Sample\tFID\tIID\tO(HOM)\tE(HOM)\tN(NM)\tF\tHet_Rate" > het/het_with_rate_summary.txt

# 处理每个 .het 文件
for het_file in het/sample_*.het; do
    # 提取基础名称（如 sample_F1SOSS_dip_Marker_1065）
    base=$(basename "$het_file" .het)
    # 去掉前缀 "sample_" 得到 "F1SOSS_dip_Marker_1065"
    sample_id=${base#sample_}
    
    # 输出包含杂合率的新文件，文件名：het/F1SOSS_dip_Marker_1065_het_with_rate.het
    out_file="het/${sample_id}_het_with_rate.het"
    
    # 添加表头
    echo -e "FID\tIID\tO(HOM)\tE(HOM)\tN(NM)\tF\tHet_Rate" > "$out_file"
    
    # 处理数据行（跳过表头）
    awk 'NR>1 {
        total=$5; hom=$3; 
        het_rate = (total - hom) / total;
        print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"het_rate;
    }' "$het_file" >> "$out_file"
    
    # 追加到汇总文件（同样添加样品标识）
    awk -v sample="$sample_id" 'NR>1 {
        total=$5; hom=$3; het_rate=(total-hom)/total;
        print sample"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"het_rate;
    }' "$het_file" >> het/het_with_rate_summary.txt
    
    echo "Generated: $out_file"
done

echo ""
echo "All jobs done. Results in ./het/"
echo "  - sample_*.het : raw PLINK heterozygosity output"
echo "  - site_*.hwe   : Hardy-Weinberg equilibrium results"
echo "  - *_het_with_rate.het : individual .het files with added Het_Rate column"
echo "  - het_with_rate_summary.txt : combined summary for all samples and SNP sets"
echo ""
echo "=== Preview of summary ==="
head -5 het/het_with_rate_summary.txt
