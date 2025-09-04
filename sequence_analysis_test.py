from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import Counter

# 读取FASTA格式的基因序列
def read_fasta(file_path):
    records = list(SeqIO.parse(file_path, "fasta"))
    return records[0].seq if records else ""

# 计算GC含量
def gc_content(seq):
    g = seq.count('G') + seq.count('g')
    c = seq.count('C') + seq.count('c')
    return (g + c) / len(seq) * 100 if len(seq) > 0 else 0

# 检测所有可能的ORF（以ATG为起始密码子，TAG/TAA/TGA为终止密码子）
def find_orfs(seq):
    orfs = []
    for frame in range(3):
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3].upper()
            if codon == "ATG":  # 找到起始密码子
                for j in range(i+3, len(seq) - 2, 3):
                    stop_codon = seq[j:j+3].upper()
                    if stop_codon in ["TAG", "TAA", "TGA"]:
                        orf = seq[i:j+3]
                        orfs.append((orf, len(orf)))
                        break
    return orfs

if __name__ == "__main__":

    sequence = read_fasta("example.fasta") or \
        "ATGGCGCTAGTAACTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG