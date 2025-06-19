"""
K-mer处理模块
包含k-mer哈希计算和相似性评估功能
"""

def get_hashes(s, k):
    """计算序列中所有 k-mer 的滚动哈希值集合。"""
    hashes = set()
    n = len(s)
    if n < k:
        return hashes
    
    base = 256  # 哈希计算的基数
    mod = 1000000007  # 大素数模数，减少哈希冲突

    # 计算第一个 k-mer 的哈希值
    current_hash = 0
    for i in range(k):
        current_hash = (current_hash * base + ord(s[i])) % mod
    hashes.add(current_hash)

    # 预计算 base^(k-1) % mod 用于滚动哈希
    base_k_minus_1 = pow(base, k - 1, mod)

    # 使用滚动哈希计算后续 k-mer 的哈希值
    for i in range(1, n - k + 1):
        current_hash = (current_hash - ord(s[i - 1]) * base_k_minus_1) % mod
        current_hash = (current_hash * base + ord(s[i + k - 1])) % mod
        hashes.add(current_hash)
    return hashes

def count_common_kmers(seq1, seq2, k=9):
    """计算两个序列之间共同 k-mer 的数量。"""
    if not seq1 or not seq2 or len(seq1) < k or len(seq2) < k:
        return 0
    hashes1 = get_hashes(seq1, k)
    hashes2 = get_hashes(seq2, k)
    common_hashes = hashes1.intersection(hashes2)
    return len(common_hashes)
