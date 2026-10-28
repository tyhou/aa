#!/bin/bash

# 切换到脚本所在目录
cd "$(dirname "$(readlink -f "$0")")"

echo "Working directory: $(pwd)"

# 处理样本
for i in F1SOSS SOLine F1-bulk-SOSS; do
    echo "=== Processing sample: $i ==="
    
    # 提取样本特定的 VCF
    echo "  Extracting sample $i from VCF..."
    /share/home/hty/software/bcftools-1.10/bcftools view -S ${i}.txt all_hap.id.vcf.gz -Ov > ${i}_hap.id.vcf
    
    # 压缩并索引（-f 强制覆盖）
    echo "  Compressing and indexing..."
    bgzip -f ${i}_hap.id.vcf
    tabix -p vcf ${i}_hap.id.vcf.gz
    
    # 提取 1065 个 SNP
    echo "  Extracting 1065 SNPs (Marker_1065.txt)..."
    /share/home/hty/software/bcftools-1.10/bcftools view -R Marker_1065.txt ${i}_hap.id.vcf.gz > ${i}_hap_Marker_1065.vcf
    
    # 提取 1391 个重叠 SNP
    echo "  Extracting 1391 SNPs (overlap_1391.txt)..."
    /share/home/hty/software/bcftools-1.10/bcftools view -R overlap_1391.txt ${i}_hap.id.vcf.gz > ${i}_hap_overlap_1391.vcf
    
    echo "  Finished sample: $i"
done

echo "=== All samples processed ==="
echo "Output files:"
ls -lh *_hap_Marker_1065.vcf *_hap_overlap_1391.vcf 2>/dev/null

# 统计 SNP 数量
echo ""
echo "=== SNP Count Statistics ==="
for f in *_hap_Marker_1065.vcf *_hap_overlap_1391.vcf; do
    if [ -f "$f" ] && [ -s "$f" ]; then
        count=$(grep -v "^#" "$f" | wc -l)
        echo "  $f: $count SNPs"
    else
        echo "  $f: EMPTY or NOT FOUND"
    fi
done
