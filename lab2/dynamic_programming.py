"""
动态规划路径查找模块
实现最长对齐路径的计算和路径回溯
"""

from tqdm import tqdm

def calculate_longest_path(nodes, switching_penalty=1000, collinear_jump_penalty=10, chunk_sz=50):
    """
    使用动态规划计算最长对齐路径
    
    参数:
        nodes: 所有匹配节点的列表
        switching_penalty: 匹配类型切换惩罚
        collinear_jump_penalty: 共线跳跃惩罚
        chunk_sz: 序列块大小
    
    返回:
        tuple: (longest_path_nodes, max_overall_score)
    """
    longest_path_nodes = []
    max_overall_score = -float('inf')
    
    if not nodes:
        print("没有找到匹配节点，无法计算最长路径。")
        return longest_path_nodes, max_overall_score
    
    # 按 q_idx, r_idx, type 排序 (拓扑排序)
    nodes.sort(key=lambda n: (n['q'], n['r'], n['type']))
    
    N_nodes = len(nodes)
    # dp[i] 存储以 nodes[i] 结尾的路径的最大分数
    dp = [node['original_score'] for node in nodes] 
    predecessor = [-1] * N_nodes  # 用于路径回溯的前驱节点索引

    # 动态规划计算最长路径
    for i in tqdm(range(N_nodes), desc="计算最长路径 (DP)"):
        n_i = nodes[i] # 当前节点 i
        q_i_idx, r_i_pos, type_i = n_i['q'], n_i['r'], n_i['type']
        
        # 遍历所有之前的节点 j 作为节点 i 的潜在前驱
        for j in range(i):
            n_j = nodes[j] # 前驱节点 j
            q_j_idx, r_j_pos, type_j = n_j['q'], n_j['r'], n_j['type']

            # 有效转换要求查询块索引增加
            if q_i_idx > q_j_idx:
                current_penalty = 0 # 当前转换的惩罚
                # 将块索引转换为碱基对坐标
                q_i_bp = q_i_idx * chunk_sz 
                q_j_bp = q_j_idx * chunk_sz 
                delta_q_bp = q_i_bp - q_j_bp # 查询序列上碱基对的差值

                if type_i != type_j: # 如果匹配类型切换（例如，正向到反向互补）
                    current_penalty = switching_penalty
                else: # 匹配类型相同
                    is_collinear = False # 是否共线
                    if type_i == 'forward': # 正向匹配
                        expected_r_i_pos = r_j_pos + delta_q_bp # 期望的参考序列起始位置
                        if r_i_pos == expected_r_i_pos:
                            is_collinear = True
                    else: # 反向互补匹配
                        expected_r_i_pos = r_j_pos - delta_q_bp # 期望的参考序列起始位置
                        if r_i_pos == expected_r_i_pos:
                            is_collinear = True
                    
                    is_q_adj = (q_i_idx == q_j_idx + 1) # 查询块是否相邻

                    if is_q_adj: # 如果查询块相邻
                        if not is_collinear: # 但不共线
                            current_penalty = switching_penalty # 视为类型切换
                    else: # 如果查询块不相邻 (q_i_idx > q_j_idx + 1)
                        if is_collinear: # 但共线
                            current_penalty = collinear_jump_penalty # 共线跳跃惩罚
                        else: # 既不相邻也不共线
                            current_penalty = switching_penalty # 视为类型切换
            
                # 候选分数 = 路径到j的分数 + 当前节点i的原始分数 - 惩罚
                potential_score = dp[j] + n_i['original_score'] - current_penalty
                
                if potential_score > dp[i]: # 如果找到更高分数的路径
                    dp[i] = potential_score
                    predecessor[i] = j # 更新前驱
    
    # 找到 DP 表中的最大分数及其对应的结束节点
    max_dp_s = -float('inf')
    max_s_idx = -1
    if N_nodes > 0:
        for dp_i in range(N_nodes):
            if dp[dp_i] > max_dp_s:
                max_dp_s = dp[dp_i]
                max_s_idx = dp_i
        max_overall_score = max_dp_s

        # 回溯路径
        if max_s_idx != -1:
            curr_idx = max_s_idx
            while curr_idx != -1:
                p_node = nodes[curr_idx].copy()
                p_node['path_score_at_node'] = dp[curr_idx] # 存储节点处的累积路径分数
                longest_path_nodes.append(p_node)
                curr_idx = predecessor[curr_idx]
            longest_path_nodes.reverse() # 反转得到正确的路径顺序
    
    return longest_path_nodes, max_overall_score

def merge_path_segments(longest_path_nodes, query_s_len, ref_s_len, chunk_sz, allowed_chunk_gap=0):
    """
    合并最长路径中的连续片段
    
    参数:
        longest_path_nodes: 最长路径上的节点列表
        query_s_len: 查询序列长度
        ref_s_len: 参考序列长度  
        chunk_sz: 序列块大小
        allowed_chunk_gap: 允许合并的最大查询块间隔
    
    返回:
        list: 合并后的片段坐标列表
    """
    merged_segments_output = []
    
    if not longest_path_nodes: 
        return merged_segments_output
    
    curr_seg = None # 当前正在合并的片段
    for i in range(len(longest_path_nodes)):
        node = longest_path_nodes[i]
        q_idx_seg = node['q']
        r_pos_seg = node['r'] # 参考块的0-based起始位置
        node_type = node['type']

        # 计算当前块的0-based开区间坐标
        q_s_bp = q_idx_seg * chunk_sz
        q_e_bp = min((q_idx_seg + 1) * chunk_sz, query_s_len)
        r_s_bp = r_pos_seg 
        r_e_bp = min(r_pos_seg + chunk_sz, ref_s_len) 

        if curr_seg is None: # 如果是第一个节点，开始新的合并片段
            curr_seg = {
                'q_start_merged_0based': q_s_bp,
                'q_end_merged_0based_open': q_e_bp,
                'ref_s_merged_0based': r_s_bp, 
                'ref_e_merged_0based_open': r_e_bp,
                'type': node_type,
                'last_q_idx': q_idx_seg,
                'last_r_idx': r_pos_seg 
            }
        else: # 尝试与前一个节点合并
            is_collinear_cont = False # 是否为共线延续
            delta_q_c = q_idx_seg - curr_seg['last_q_idx'] # 查询块索引的差值
            
            # 检查是否可以合并：类型相同，且查询块间隔在允许范围内
            if (node_type == curr_seg['type'] and
                0 < delta_q_c <= (allowed_chunk_gap + 1)):
                
                exp_r_pos = -1 # 期望的参考块起始位置
                if node_type == 'forward':
                    exp_r_pos = curr_seg['last_r_idx'] + delta_q_c * chunk_sz
                elif node_type == 'rc':
                    exp_r_pos = curr_seg['last_r_idx'] - delta_q_c * chunk_sz
                
                if r_pos_seg == exp_r_pos: # 如果参考块位置符合预期，则为共线延续
                    is_collinear_cont = True
            
            if is_collinear_cont: # 如果可以合并
                curr_seg['q_end_merged_0based_open'] = q_e_bp # 更新查询结束位置
                # 更新参考序列的起始和结束位置，取合并片段的最小起始和最大结束
                curr_seg['ref_s_merged_0based'] = min(curr_seg['ref_s_merged_0based'], r_s_bp)
                curr_seg['ref_e_merged_0based_open'] = max(curr_seg['ref_e_merged_0based_open'], r_e_bp)
                curr_seg['last_q_idx'] = q_idx_seg
                curr_seg['last_r_idx'] = r_pos_seg
            else: # 如果不能合并，则完成前一个片段，开始新的片段
                merged_segments_output.append(
                    (curr_seg['q_start_merged_0based'], 
                     curr_seg['q_end_merged_0based_open'],
                     curr_seg['ref_s_merged_0based'], 
                     curr_seg['ref_e_merged_0based_open']))
                curr_seg = { # 开始新的合并片段
                    'q_start_merged_0based': q_s_bp,
                    'q_end_merged_0based_open': q_e_bp,
                    'ref_s_merged_0based': r_s_bp,
                    'ref_e_merged_0based_open': r_e_bp,
                    'type': node_type,
                    'last_q_idx': q_idx_seg,
                    'last_r_idx': r_pos_seg
                }
    
    if curr_seg: # 添加最后一个处理的片段
        merged_segments_output.append(
            (curr_seg['q_start_merged_0based'], 
             curr_seg['q_end_merged_0based_open'],
             curr_seg['ref_s_merged_0based'], 
             curr_seg['ref_e_merged_0based_open']))
    
    # 扩展最后一个片段到序列末尾
    if merged_segments_output:
        end_idx = len(merged_segments_output) - 1
        begin_q = merged_segments_output[end_idx][0]
        end_q   = merged_segments_output[end_idx][1]
        begin_r = merged_segments_output[end_idx][2]
        end_r   = merged_segments_output[end_idx][3]

        dl = min(query_s_len - end_q, ref_s_len - end_r)
        merged_segments_output = [list(segment) for segment in merged_segments_output]
        merged_segments_output[end_idx][1] = end_q + dl
        merged_segments_output[end_idx][3] = end_r + dl
    
    return merged_segments_output

def format_segments_output(merged_segments_output):
    """
    格式化片段输出为字符串形式
    
    参数:
        merged_segments_output: 合并后的片段列表
    
    返回:
        str: 格式化的坐标字符串
    """
    if not merged_segments_output:
        return "[]"
    
    fmt_tuples = [f" ( {s[0]}, {s[1]}, {s[2]}, {s[3]} )" for s in merged_segments_output]
    output_str = "["
    for i, t_str in enumerate(fmt_tuples):
        if i == 0:
            output_str += t_str
        else:
            output_str += " ,  " 
            output_str += t_str
    output_str += " ]"
    return output_str
