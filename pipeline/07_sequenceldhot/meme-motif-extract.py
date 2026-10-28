import re

# 读取 meme.txt
with open("meme.txt", "r") as f:
    content = f.read()

# 提取文件头信息
# MEME version
version_match = re.search(r"MEME version (\S+)", content)
version = version_match.group(1) if version_match else "5"

# ALPHABET
alphabet_match = re.search(r"ALPHABET=\s*(\S+)", content)
alphabet = alphabet_match.group(1) if alphabet_match else "ACGT"

# strands（从命令行摘要中获取，默认 +）
strands_match = re.search(r"strands:\s*(.+)", content)
strands = strands_match.group(1).strip() if strands_match else "+"

# Background letter frequencies
bg_match = re.search(
    r"Background letter frequencies.*?\nA\s+([\d.]+)\s+C\s+([\d.]+)\s+G\s+([\d.]+)\s+T\s+([\d.]+)",
    content, re.DOTALL
)
if bg_match:
    bg_line = f"A {bg_match.group(1)} C {bg_match.group(2)} G {bg_match.group(3)} T {bg_match.group(4)}"
else:
    bg_line = "A 0.25 C 0.25 G 0.25 T 0.25"

# 提取每个 motif
# 匹配 MOTIF 行 + letter-probability matrix 行 + 矩阵数据
motif_pattern = re.compile(
    r"MOTIF\s+(\S+)\s+(MEME-\d+).*?"
    r"letter-probability matrix:\s*alength=\s*(\d+)\s+w=\s*(\d+)\s+nsites=\s*(\d+)\s+E=\s*(\S+)\s*\n"
    r"((?:\s*[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*\n)+)",
    re.DOTALL
)

motifs = motif_pattern.findall(content)

for idx, (name, meme_id, alength, w, nsites, evalue, matrix_block) in enumerate(motifs, 1):
    # 清理矩阵行
    matrix_lines = []
    for line in matrix_block.strip().split("\n"):
        vals = line.split()
        if len(vals) == 4:
            matrix_lines.append("  " + "  ".join(f"{float(v):.6f}" for v in vals) + " ")
    
    out = f"""MEME version {version}

ALPHABET= {alphabet}

strands: {strands}

Background letter frequencies (from dataset with add-one prior applied):
{bg_line}

MOTIF {name} {meme_id}
letter-probability matrix: alength= {alength} w= {w} nsites= {nsites} E= {evalue}
""" + "\n".join(matrix_lines) + "\n"

    filename = f"meme_motif{idx}.meme"
    with open(filename, "w") as out_f:
        out_f.write(out)
    print(f"已生成 {filename}  ({name} {meme_id})")

print(f"\n共提取 {len(motifs)} 个 motif")