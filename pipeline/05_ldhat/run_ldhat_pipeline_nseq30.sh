#!/bin/bash
set -e

# 工具路径
BCF="/share/home/hty/software/bcftools-1.10/bcftools"
LDHAT="/share/home/hty/software/LDhat-master"

echo "=== Step 0: 提取15个样本 ==="
$BCF view -S 15_SO_line.txt ../SOLine_dip_overlap_1391.vcf -Ov > 15_SOLine_dip_overlap_1391.vcf

echo "=== Step 1: 去除重复位置（保留第一次出现） ==="
# 备份原始文件（可选）
cp 15_SOLine_dip_overlap_1391.vcf 15_SOLine_dip_overlap_1391.bak.vcf

# 利用 awk 按染色体:位置去重，保留头部
$BCF view -Ov 15_SOLine_dip_overlap_1391.vcf | \
awk '!/^#/ {key=$1":"$2; if (!seen[key]++) print; next} /^#/ {print}' \
> 15_SOLine_dip_overlap_1391.uniq.vcf

mv 15_SOLine_dip_overlap_1391.uniq.vcf 15_SOLine_dip_overlap_1391.vcf

echo "=== Step 2: vcftools 转换为 LDhat 格式（8条染色体） ==="
for i in {1..8}; do
    echo "  转换染色体 ${i} ..."
    vcftools --vcf 15_SOLine_dip_overlap_1391.vcf --out chr${i} --ldhat-geno --chr ${i}
done

echo "=== Step 3: 准备似然率表 ==="
cp -r $LDHAT/lk_files/ ./
cp lk_files/lk_n50_t0.1.gz ./
gunzip -f lk_n50_t0.1.gz
$LDHAT/lkgen -lk lk_n50_t0.1 -nseq 30

echo "=== Step 4: 运行 interval ==="
for i in {1..8}; do
    echo "  interval 染色体 ${i} ..."
    $LDHAT/interval -seq chr${i}.ldhat.sites -loc chr${i}.ldhat.locs \
        -lk new_lk.txt -its 1000000 -bpen 5 -samp 2000 -prefix chr${i}.Rho.
done

echo "=== Step 5: 运行 stat ==="
for i in {1..8}; do
    echo "  stat 染色体 ${i} ..."
    $LDHAT/stat -input chr${i}.Rho.rates.txt -loc chr${i}.ldhat.locs -prefix chr${i}.stat.
done

echo "=== Step 6: 合并结果 ==="
for i in {1..8}; do
    sed "s/^/Chr${i}\t/" chr${i}.stat.res.txt | sed '1d' > chr${i}.stat.final.txt
    echo >> chr${i}.stat.final.txt
done
cat chr{1..8}.stat.final.txt | sed '1iChr\tLoci\tMean_rho\tMedian\tL95\tU95' > Rho.stat.final.txt

echo "=== 全部完成！结果在 Rho.stat.final.txt ==="