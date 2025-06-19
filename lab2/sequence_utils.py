"""
序列处理工具模块
包含DNA序列的基本操作函数
"""

def reverse_complement(seq):
    """计算 DNA 序列的反向互补序列。"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'} # N表示未知碱基
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def rc(seq):
    """计算反向互补序列的简写函数。"""
    return reverse_complement(seq)
