import os

files = [f"result_nseq30/chr{i}.stat.res.txt" for i in range(1, 9)]
out_suffix = ".with_span.txt"

summary = []
total_span_all = 0.0
total_gdist_all = 0.0

for fname in files:
    with open(fname, 'r') as f:
        lines = f.readlines()
    header = lines[0].strip().split()
    new_header = header + ["SNP_span", "Genetic_distance"]
    
    data = []
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) < 5:
            continue
        loci = float(parts[0])
        mean_rho = float(parts[1])
        median = float(parts[2])
        l95 = float(parts[3])
        u95 = float(parts[4])
        data.append([loci, mean_rho, median, l95, u95])
    
    real_indices = [i for i, row in enumerate(data) if row[0] != -1.0]
    n_real = len(real_indices)
    
    out_lines = ["\t".join(new_header)]
    total_genetic_dist = 0.0
    prev_real_loci = None
    
    for i, row in enumerate(data):
        loci = row[0]
        mean_rho = row[1]
        span = None
        gdist = None
        
        if loci != -1.0:
            if prev_real_loci is None:
                span = 0.0
                gdist = 0.0
            else:
                span = loci - prev_real_loci
                gdist = mean_rho * span
                total_genetic_dist += gdist
            prev_real_loci = loci
        else:
            span = None
            gdist = None
        
        span_str = f"{span:.2f}" if span is not None else "NA"
        gdist_str = f"{gdist:.2f}" if gdist is not None else "NA"
        row_str = "\t".join([str(x) for x in row])
        out_line = row_str + "\t" + span_str + "\t" + gdist_str
        out_lines.append(out_line)
    
    if n_real >= 2:
        first_loci = data[real_indices[0]][0]
        last_loci = data[real_indices[-1]][0]
        total_physical_span = last_loci - first_loci
    else:
        total_physical_span = 0.0
    
    out_fname = fname + out_suffix
    with open(out_fname, 'w') as f:
        f.write("\n".join(out_lines))
    
    total_span_all += total_physical_span
    total_gdist_all += total_genetic_dist
    
    avg_rho = total_genetic_dist / total_physical_span if total_physical_span > 0 else 0
    span_per_gd = total_physical_span / total_genetic_dist if total_genetic_dist > 0 else float('inf')
    
    summary.append({
        'chrom': fname.split('/')[-1].replace('.stat.res.txt', ''),
        'total_physical_span': total_physical_span,
        'total_genetic_dist': total_genetic_dist,
        'avg_rho': avg_rho,
        'span_per_gd': span_per_gd
    })
    print(f"Processed {fname} -> {out_fname}")

# 准备汇总行列表（包括基因组-wide）
summary_rows = []
for s in summary:
    spd = s['span_per_gd']
    spd_str = "inf" if spd == float('inf') else f"{spd:.2f}"
    avg_rho_str = "inf" if s['avg_rho'] == float('inf') else f"{s['avg_rho']:.2f}"
    summary_rows.append([
        s['chrom'],
        f"{s['total_physical_span']:.2f}",
        f"{s['total_genetic_dist']:.2f}",
        spd_str,
        avg_rho_str
    ])

# 全基因组行
genome_span = total_span_all
genome_gdist = total_gdist_all
genome_avg_rho = genome_gdist / genome_span if genome_span > 0 else 0
genome_span_per_gd = genome_span / genome_gdist if genome_gdist > 0 else float('inf')
summary_rows.append([
    "Genome-wide",
    f"{genome_span:.2f}",
    f"{genome_gdist:.2f}",
    f"{genome_span_per_gd:.2f}",
    f"{genome_avg_rho:.2f}"
])

# 打印到终端
print("\n=== Summary for all chromosomes ===")
print(f"{'Chromosome':<12} {'SNPs span (kb)':<20} {'Genetic distance (cM)':<22} {'Recombination rate (kb/cM)':<25} {'Recombination rate (cM/kb)':<25}")
for row in summary_rows:
    print(f"{row[0]:<12} {row[1]:<20} {row[2]:<22} {row[3]:<25} {row[4]:<25}")

# 保存到文件
with open("RR_summary_results_1391.txt", 'w') as f:
    f.write("Chromosome\tSNPs span (kb)\tGenetic distance (cM)\tRecombination rate (kb/cM)\tRecombination rate (cM/kb)\n")
    for row in summary_rows:
        f.write("\t".join(row) + "\n")

print("\nSummary results have been saved to 'RR_summary_results_1391.txt'")
