"""
序列匹配处理模块
包含序列块匹配和并行处理功能
"""

import edlib
from sequence_utils import reverse_complement

def _process_query_chunk_matching_worker(args):
    """多进程并行处理的核心工作函数"""
    q_idx, query_s, ref_s, chunk_sz, k, kmer_thresh, query_s_len, ref_s_len = args

    q_s_offset = q_idx * chunk_sz
    q_e_offset = min((q_idx + 1) * chunk_sz, query_s_len)
    q_chunk = query_s[q_s_offset:q_e_offset]
    chunk_sz = q_e_offset - q_s_offset  # 实际查询块大小

    fwd_matches = []
    rc_matches = []

    if len(q_chunk) < k:
        return fwd_matches, rc_matches

    # 正向匹配
    for r_start in range(ref_s_len - chunk_sz + 1):
        ref_chunk = ref_s[r_start : r_start + chunk_sz]
        # 使用 edlib 计算编辑距离作为匹配分数
        edit_distance = edlib.align(q_chunk, ref_chunk)['editDistance']
        score = chunk_sz - edit_distance # 分数越高，匹配越好
        if score >= kmer_thresh:
            fwd_matches.append((q_idx, r_start, score))

    # 反向互补匹配
    rc_q_chunk = reverse_complement(q_chunk)

    for r_start in range(ref_s_len - chunk_sz + 1):
        ref_chunk = ref_s[r_start : r_start + chunk_sz]
        edit_distance = edlib.align(rc_q_chunk, ref_chunk)['editDistance']
        score = chunk_sz - edit_distance # 分数越高，匹配越好
        if score >= kmer_thresh:
            rc_matches.append((q_idx, r_start, score))
            
    if not(fwd_matches or rc_matches):
        # 如果没有找到任何匹配，打印提示信息
        print(f"查询块 {q_idx} (长度 {len(q_chunk)}) 在参考序列中未找到匹配。")
        print(f"序列坐标：{q_s_offset}，参考序列长度: {ref_s_len}bp。")
    
    return fwd_matches, rc_matches

def process_matches_for_visualization(raw_matches, match_type):
    """
    处理原始匹配列表，识别共线段和孤立点，用于可视化。

    参数:
        raw_matches (list): 原始匹配点列表，每个元素为 (q_idx, r_idx, score)。
        match_type (str): 匹配类型，"forward" 或 "rc"。

    返回:
        tuple: (segments, isolated_points_data)
               segments: 共线段列表，每个元素为 (q_start, r_start, q_end, r_end)。
               isolated_points_data: 孤立点列表，每个元素为 (q, r, score)。
    """
    if not raw_matches:
        return [], []
    
    # 将匹配点及其分数存储在字典中，方便查找
    points_with_scores = {(m[0], m[1]): m[2] for m in raw_matches}
    available_points = set(points_with_scores.keys()) # 所有可用匹配点
    
    segments = [] # 存储识别出的共线段
    isolated_points_data = [] # 存储识别出的孤立点
    
    # 按 q_idx, r_idx 排序，确保处理顺序一致
    sorted_unique_points = sorted(list(available_points), key=lambda p: (p[0], p[1]))
    
    visited_in_segment = set() # 记录已包含在某个段中的点

    for q_start, r_start in sorted_unique_points:
        if (q_start, r_start) in visited_in_segment: # 如果点已处理过，则跳过
            continue
        
        current_segment_path = [(q_start, r_start)] # 当前段的路径
        q_curr, r_curr = q_start, r_start

        # 尝试扩展当前段
        while True:
            if match_type == "forward": # 正向匹配，下一个点应该是 (q+1, r+1)
                q_next, r_next = q_curr + 1, r_curr + 1
            else: # 反向互补匹配，下一个点应该是 (q+1, r-1)
                q_next, r_next = q_curr + 1, r_curr - 1
            
            # 如果下一个点存在且未被访问过
            if (q_next, r_next) in available_points and (q_next, r_next) not in visited_in_segment:
                current_segment_path.append((q_next, r_next))
                q_curr, r_curr = q_next, r_next
            else:
                break # 无法扩展或下一点已被使用
        
        if len(current_segment_path) > 1: # 如果路径长度大于1，则认为是一个段
            seg_q_start_coord, seg_r_start_coord = current_segment_path[0]
            seg_q_end_coord, seg_r_end_coord = current_segment_path[-1]
            segments.append((seg_q_start_coord, seg_r_start_coord, seg_q_end_coord, seg_r_end_coord))
            for p_item in current_segment_path: # 将段中的所有点标记为已访问
                visited_in_segment.add(p_item)
                
    # 收集孤立点：任何未包含在段中的原始匹配点
    for q, r in sorted_unique_points:
        if (q, r) not in visited_in_segment:
            isolated_points_data.append((q, r, points_with_scores[(q,r)]))
            
    return segments, isolated_points_data
