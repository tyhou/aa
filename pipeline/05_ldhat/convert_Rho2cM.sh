#!/bin/bash
# 功能：将 LDhat stat 结果转换为遗传距离（cM）
# 输入：result_nseq30 目录下的 chr*.stat.res.txt 文件
# 输出：chr*.genetic_map.txt 和 genetic_lengths.txt

set -e

# ================== 配置 ==================
Ne=39.604
IN_DIR="result_nseq30"
OUT_DIR="."
summary_file="$OUT_DIR/genetic_lengths.txt"
# ==========================================

if [ ! -d "$IN_DIR" ]; then
    echo "错误：目录 $IN_DIR 不存在！"
    exit 1
fi

echo -e "Chr\tTotal_bp\tAvg_rho\tTotal_cM" > "$summary_file"

for chr in {1..8}; do
    infile="$IN_DIR/chr${chr}.stat.res.txt"
    outfile="$OUT_DIR/chr${chr}.genetic_map.txt"

    if [ ! -f "$infile" ]; then
        echo "警告：$infile 不存在，跳过染色体 $chr"
        continue
    fi

    echo "处理染色体 $chr ..."

    awk -v Ne="$Ne" -v chr="$chr" '
    BEGIN {
        print "Loci\tcum_cM" > "'$outfile'"
        prev_pos = 0
        prev_rho = 0
        cum_cM = 0
        line_count = 0
        sum_rho = 0
        rho_count = 0
        min_pos = 0
        max_pos = 0
        first = 1
    }
    NR == 1 { next }                           # 跳过表头
    $1 == -1 { next }                          # 跳过全局统计行
    {
        pos = $1          # 物理位置
        rho = $2          # Mean_rho

        # 记录第一个有效位点
        if (first) {
            min_pos = pos
            prev_pos = pos
            prev_rho = rho
            first = 0
            sum_rho += rho
            rho_count++
            next
        }

        # 检查位置是否严格递增
        if (pos <= prev_pos) {
            print "错误：位置未排序，pos=" pos ", prev_pos=" prev_pos > "/dev/stderr"
            exit 1
        }

        # 计算区间遗传距离
        dist_bp = pos - prev_pos
        avg_rho = (rho + prev_rho) / 2
        cM = (avg_rho / (4 * Ne)) * dist_bp * 100
        cum_cM += cM

        # 输出当前位点累积 cM
        print pos "\t" cum_cM >> "'$outfile'"

        # 更新
        prev_pos = pos
        prev_rho = rho
        sum_rho += rho
        rho_count++
    }
    END {
        max_pos = pos
        total_bp = max_pos - min_pos
        if (rho_count > 0) avg_rho_all = sum_rho / rho_count
        else avg_rho_all = 0

        # 写入汇总文件（注意这里用 >> 追加）
        print chr "\t" total_bp "\t" avg_rho_all "\t" cum_cM >> "'$summary_file'"
        print "  -> 总物理长度: " total_bp " bp, 总遗传长度: " cum_cM " cM"
    }
    ' "$infile"

    echo "染色体 $chr 处理完成，输出文件 $outfile"
done

echo "全部完成！汇总结果保存在 $summary_file"
