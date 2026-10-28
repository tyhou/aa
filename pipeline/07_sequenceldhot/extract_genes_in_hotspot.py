#!/usr/bin/env python3
# extract_genes_in_hotspot.py
import glob, os

gff_file = "IId_11730_genome.gff3"
out_detail = "hotspot_covered_genes.txt"
out_ids    = "hotspot_gene_ids.txt"

# ---------- 1. 读取热点区间 ----------
# hotspots[chrom] = [(start, stop), ...]   chrom 用 "1".."8"
hotspots = {}
for f in glob.glob("*_hotspot.txt"):
    chrom = os.path.basename(f).split("_hotspot")[0]   # "chr1"
    chrom_key = chrom.replace("chr", "")               # "1"
    with open(f) as fh:
        fh.readline()  # 跳过表头
        for line in fh:
            line = line.strip()
            if not line:
                continue
            p = line.split()
            try:
                s, e = int(p[0]), int(p[1])   # Hotspot Start/Stop
            except ValueError:
                continue
            if s > e:
                s, e = e, s
            hotspots.setdefault(chrom_key, []).append((s, e))

print("热点区间统计:")
for c in sorted(hotspots, key=lambda x: int(x)):
    print(f"  chr{c}: {len(hotspots[c])} 个区间")

# ---------- 2. 解析 GFF，取 gene 行 ----------
def overlaps(a_s, a_e, b_s, b_e):
    return a_s <= b_e and b_s <= a_e

results = []
gene_count = 0

with open(gff_file) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        chrom = cols[0]
        ftype = cols[2]
        start = int(cols[3])
        end   = int(cols[4])
        strand = cols[6]
        attrs = cols[8]

        if ftype != "gene":
            continue

        # 染色体匹配：GFF 是 "1"，热点是 "chr1" -> 统一成 "1"
        chrom_key = chrom.replace("chr", "")
        if chrom_key not in hotspots:
            continue

        gene_count += 1

        # 解析 ID 和 locus_tag
        gene_id = None
        locus_tag = None
        for kv in attrs.split(";"):
            kv = kv.strip()
            if kv.startswith("ID="):
                gene_id = kv[3:]
            elif kv.startswith("locus_tag="):
                locus_tag = kv[10:]

        # 优先用 locus_tag（CPHLJ_xxx），否则从 ID 里去掉 "gene-" 前缀
        if locus_tag:
            gname = locus_tag
        elif gene_id:
            gname = gene_id.replace("gene-", "")
        else:
            gname = f"{chrom}:{start}-{end}"

        hits = [(s, e) for (s, e) in hotspots[chrom_key]
                if overlaps(start, end, s, e)]
        if hits:
            results.append((gname, chrom, start, end, strand, hits))

# ---------- 3. 输出 ----------
with open(out_detail, "w") as out:
    out.write("GeneID\tChrom\tStart\tEnd\tStrand\tN_hotspots\tHotspot_regions\n")
    for gname, chrom, s, e, strand, hits in results:
        hs = ",".join(f"{a}-{b}" for a, b in hits)
        out.write(f"{gname}\t{chrom}\t{s}\t{e}\t{strand}\t{len(hits)}\t{hs}\n")

with open(out_ids, "w") as out:
    for gname, *_ in results:
        out.write(gname + "\n")

print(f"\n共扫描 gene: {gene_count}")
print(f"被热点覆盖的基因数: {len(results)}")
print(f"详细结果 -> {out_detail}")
print(f"纯基因ID -> {out_ids}")