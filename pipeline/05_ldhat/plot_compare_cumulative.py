import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ============ 设置路径 ============
dir_1065 = "before_1065_result_nseq30"
dir_1391 = "result_nseq30"
chroms = [f"chr{i}" for i in range(1, 9)]

# ============ 存储数据用于绘图 ============
all_data = {}

for chrom in chroms:
    file1065 = os.path.join(dir_1065, f"{chrom}.stat.res.txt.with_span.txt")
    file1391 = os.path.join(dir_1391, f"{chrom}.stat.res.txt.with_span.txt")
    
    # 检查文件是否存在
    if not os.path.exists(file1065):
        print(f"Warning: {file1065} not found. Skipping.")
        continue
    if not os.path.exists(file1391):
        print(f"Warning: {file1391} not found. Skipping.")
        continue
    
    # 读取数据，跳过头行
    df1065 = pd.read_csv(file1065, sep='\t', skiprows=0)
    df1391 = pd.read_csv(file1391, sep='\t', skiprows=0)
    
    # 过滤掉全局平均行（Loci = -1），并按物理位置排序
    df1065 = df1065[df1065['Loci'] != -1.0].sort_values('Loci').reset_index(drop=True)
    df1391 = df1391[df1391['Loci'] != -1.0].sort_values('Loci').reset_index(drop=True)
    
    # 计算累积遗传距离
    df1065['cum_dist'] = df1065['Genetic_distance'].cumsum()
    df1391['cum_dist'] = df1391['Genetic_distance'].cumsum()
    
    # 存储原始数据用于后续插值
    all_data[chrom] = {
        '1065': df1065[['Loci', 'cum_dist']].copy(),
        '1391': df1391[['Loci', 'cum_dist']].copy()
    }

# ============ 绘制图形 ============

# ---- 图1: 累积遗传距离对比 (2行 x 4列) ----
fig1, axes1 = plt.subplots(2, 4, figsize=(20, 10))
axes1 = axes1.flatten()

for idx, chrom in enumerate(chroms):
    ax = axes1[idx]
    data1065 = all_data[chrom]['1065']
    data1391 = all_data[chrom]['1391']
    
    ax.plot(data1065['Loci'], data1065['cum_dist'], 
            label='1065 SNPs', color='blue', linewidth=2)
    ax.plot(data1391['Loci'], data1391['cum_dist'], 
            label='1391 variants', color='red', linestyle='--', linewidth=2)
    
    ax.set_title(f'{chrom}', fontsize=14)
    ax.set_xlabel('Physical position (kb)', fontsize=10)
    ax.set_ylabel('Cumulative genetic distance (cM)', fontsize=10)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # 显示最终差异值在图上
    final_diff = data1391['cum_dist'].iloc[-1] - data1065['cum_dist'].iloc[-1]
    ax.text(0.05, 0.95, f'Δ终值 = {final_diff:.2f} cM', 
            transform=ax.transAxes, fontsize=10, 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('cumulative_comparison_8chroms.png', dpi=300)
print("Saved: cumulative_comparison_8chroms.png")

# ---- 图2: 差异累积曲线 (1391 - 1065)，显示差距贡献区域 ----
fig2, axes2 = plt.subplots(2, 4, figsize=(20, 10))
axes2 = axes2.flatten()

for idx, chrom in enumerate(chroms):
    ax = axes2[idx]
    df1065 = all_data[chrom]['1065']
    df1391 = all_data[chrom]['1391']
    
    # 取并集Loci用于插值对齐
    union_loci = np.unique(np.concatenate([df1065['Loci'].values, df1391['Loci'].values]))
    # 插值得到同一位置上的累积距离
    cum1065_interp = np.interp(union_loci, df1065['Loci'], df1065['cum_dist'])
    cum1391_interp = np.interp(union_loci, df1391['Loci'], df1391['cum_dist'])
    diff = cum1391_interp - cum1065_interp
    
    # 绘制差异曲线
    ax.plot(union_loci, diff, color='darkgreen', linewidth=2)
    ax.axhline(y=0, color='black', linestyle=':', linewidth=1)  # 零线
    
    # 用不同颜色填充正负区域
    ax.fill_between(union_loci, 0, diff, where=(diff > 0), color='green', alpha=0.3, label='1391 > 1065')
    ax.fill_between(union_loci, 0, diff, where=(diff < 0), color='red', alpha=0.3, label='1391 < 1065')
    
    ax.set_title(f'{chrom} 差异累积 (1391 - 1065)', fontsize=14)
    ax.set_xlabel('Physical position (kb)', fontsize=10)
    ax.set_ylabel('Cumulative difference (cM)', fontsize=10)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # 标记最终差异
    final_diff = diff[-1]
    ax.text(0.05, 0.95, f'总差异 = {final_diff:.2f} cM', 
            transform=ax.transAxes, fontsize=10, 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

plt.tight_layout()
plt.savefig('cumulative_difference_8chroms.png', dpi=300)
print("Saved: cumulative_difference_8chroms.png")

plt.show()
print("\nAll plots generated successfully!")