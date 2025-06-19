import os
import multiprocessing
from tqdm import tqdm

# 导入模块化的组件
from sequence_utils import reverse_complement, rc
from kmer_utils import get_hashes, count_common_kmers
from matching import _process_query_chunk_matching_worker, process_matches_for_visualization
from scoring import get_points, calculate_distance, calculate_value
from visualization import generate_svg_grid
from dynamic_programming import calculate_longest_path, merge_path_segments, format_segments_output

# --- 主要处理逻辑 ---
def find_and_visualize_matches(query_f, ref_f, svg_out_f,
                               k=7, chunk_sz=50, kmer_thresh=5,
                               switching_penalty=1000, collinear_jump_penalty=10,
                               allowed_chunk_gap=0):
    try:
        with open(query_f, 'r') as f_in:
            query_s = f_in.read().replace('\\\\n', '').strip().upper() 
        with open(ref_f, 'r') as f_in:
            ref_s = f_in.read().replace('\\\\n', '').strip().upper()
    except FileNotFoundError as e:
        print(f"错误: {e}。请确保输入文件存在。")
        return

    if not query_s or not ref_s:
        print("错误: 查询序列或参考序列为空。")
        return

    query_s_len = len(query_s)
    ref_s_len = len(ref_s)

    fwd_matches_raw = [] # 存储原始正向匹配
    rc_matches_raw = []  # 存储原始反向互补匹配

    num_q_chunks = (query_s_len + chunk_sz - 1) // chunk_sz # 查询序列的块数
    ref_chunks_svg = (ref_s_len + chunk_sz - 1) // chunk_sz # 参考序列的块数（用于SVG）

    print(f"正在处理查询序列 ({query_s_len}bp) 和参考序列 ({ref_s_len}bp)...")
    print(f"K-mer 大小: {k}, 块大小: {chunk_sz}, K-mer 阈值: {kmer_thresh}, 切换惩罚: {switching_penalty}, 共线跳跃惩罚: {collinear_jump_penalty}, 允许的块间隔: {allowed_chunk_gap}") 
    print("注意: 块匹配分数基于edlib编辑距离。")

    # 准备并行处理的任务参数
    tasks = []
    for i in range(num_q_chunks):
        tasks.append((i, query_s, ref_s, chunk_sz, k, kmer_thresh, query_s_len, ref_s_len))

    # 尝试使用多进程并行处理查询块，如果失败则使用单线程处理
    try:
        with multiprocessing.Pool(processes=min(4, multiprocessing.cpu_count())) as pool:
            pool_res = pool.imap_unordered(_process_query_chunk_matching_worker, tasks)
            # 使用 tqdm 显示进度条
            for f_match_list, rc_match_list in tqdm(pool_res, total=len(tasks), desc="并行处理查询块"):
                if f_match_list:
                    fwd_matches_raw.extend(f_match_list)
                if rc_match_list:
                    rc_matches_raw.extend(rc_match_list)
    except (RuntimeError, OSError) as e:
        print(f"多进程处理失败 ({e})，切换到单线程处理...")
        # 单线程处理
        for task in tqdm(tasks, desc="单线程处理查询块"):
            f_match_list, rc_match_list = _process_query_chunk_matching_worker(task)
            if f_match_list:
                fwd_matches_raw.extend(f_match_list)
            if rc_match_list:
                rc_matches_raw.extend(rc_match_list)
    
    print(f"找到 {len(fwd_matches_raw)} 个正向匹配和 {len(rc_matches_raw)} 个反向互补匹配 (基于阈值)。")

    # --- 最长路径计算 ---
    nodes = [] # 存储所有匹配节点信息
    # r_node_pos 是参考块的起始位置 (0-based)
    for q_node_idx, r_node_pos, score in fwd_matches_raw:
        nodes.append({'q': q_node_idx, 'r': r_node_pos, 'score': float(score), 'type': 'forward', 'original_score': float(score)})
    for q_node_idx, r_node_pos, score in rc_matches_raw:
        nodes.append({'q': q_node_idx, 'r': r_node_pos, 'score': float(score), 'type': 'rc', 'original_score': float(score)})

    longest_path_nodes, max_overall_score = calculate_longest_path(
        nodes, switching_penalty, collinear_jump_penalty, chunk_sz
    )

    print(f"最长对齐路径分数: {max_overall_score if nodes else 'N/A'}")
    
    calculated_value_score = 0
    if longest_path_nodes:
        merged_segments_output = merge_path_segments(
            longest_path_nodes, query_s_len, ref_s_len, chunk_sz, allowed_chunk_gap
        )
        
        if not merged_segments_output:
            print("[]") # 如果没有合并的片段，打印空列表
        else:
            output_str = format_segments_output(merged_segments_output)
            print(output_str)
            # 使用 calculate_value 计算最终分数
            calculated_value_score = calculate_value(output_str, ref_s, query_s)
            print(f"最长路径的计算分值: {calculated_value_score}")
    else: # 如果没有找到最长路径
        print("未找到有效的最长对齐路径。")
        print("[]")
        calculated_value_score = calculate_value("[]", ref_s, query_s) # 对空路径计算分数应为0
        print(f"计算分值 (无路径): {calculated_value_score}")

    # --- 为SVG可视化准备数据 ---
    # 将原始匹配中的 ref_chunk_start (碱基对位置) 转换为 ref_chunk_idx (块索引)
    fwd_matches_grid = []
    if fwd_matches_raw:
        fwd_matches_grid = [(q_idx, r_start // chunk_sz, score) 
                                        for q_idx, r_start, score in fwd_matches_raw]

    rc_matches_grid = []
    if rc_matches_raw:
        rc_matches_grid = [(q_idx, r_start // chunk_sz, score) 
                                   for q_idx, r_start, score in rc_matches_raw]

    # 处理转换后的匹配以识别线段和孤立点 (用于SVG网格)
    fwd_segs_grid, fwd_pts_grid = process_matches_for_visualization(fwd_matches_grid, "forward")
    rc_segs_grid, rc_pts_grid = process_matches_for_visualization(rc_matches_grid, "rc")

    print(f"处理后用于分段显示: {len(fwd_segs_grid)} 个正向线段, {len(fwd_pts_grid)} 个正向孤立点")
    print(f"处理后用于分段显示: {len(rc_segs_grid)} 个反向互补线段, {len(rc_pts_grid)} 个反向互补孤立点")

    # 为最长路径高亮转换坐标 (将参考序列的碱基对位置转换为块索引)
    long_path_grid = []
    if longest_path_nodes: # longest_path_nodes 中的 'r' 是 ref_chunk_start (碱基对位置)
        long_path_grid = [
            {
                'q': node['q'], 
                'r': node['r'] // chunk_sz, # 转换为 ref_chunk_idx
                'type': node['type'],
                'original_score': node['original_score'],
                'path_score_at_node': node.get('path_score_at_node', 'N/A')
            } for node in longest_path_nodes
        ]

    generate_svg_grid(query_s_len, ref_s_len, 
                      fwd_segs_grid, fwd_pts_grid,
                      rc_segs_grid, rc_pts_grid,
                      chunk_sz,
                      num_q_chunks, ref_chunks_svg,
                      svg_out_f,
                      long_path=long_path_grid) # 传递转换后的最长路径数据

def main():
    """默认主函数，可以直接运行进行测试"""
    # 定义实际使用的文件路径
    use_query_f = "query.txt"
    use_ref_f = "ref.txt"
    svg_path = "matches_grid_visualization.svg"

    find_and_visualize_matches(
        query_f=use_query_f,
        ref_f=use_ref_f,
        svg_out_f=svg_path,
        k=6,                # K-mer 大小
        chunk_sz=25,        # 比较的序列块大小
        kmer_thresh=22,     # K-mer 匹配阈值 (基于编辑距离的分数)
        switching_penalty=10, # 切换匹配类型的惩罚
        collinear_jump_penalty=9, # 共线跳跃的惩罚
        allowed_chunk_gap=18    # 允许合并的共线段之间的最大查询块间隔
    )

def run_default_analysis():
    """
    运行默认的DNA序列匹配分析。
    使用预设参数分析query.txt和ref.txt文件。
    """
    # 定义实际使用的文件路径
    use_query_f = "query.txt"
    use_ref_f = "ref.txt"
    svg_path = "matches_grid_visualization.svg"

    # 检查文件是否存在
    import os
    if os.path.exists(use_query_f) and os.path.exists(use_ref_f):
        find_and_visualize_matches(
            query_f=use_query_f,
            ref_f=use_ref_f,
            svg_out_f=svg_path,
            k=6,                # K-mer 大小
            chunk_sz=25,        # 比较的序列块大小
            kmer_thresh=22,     # K-mer 匹配阈值 (基于编辑距离的分数)
            switching_penalty=10, # 切换匹配类型的惩罚
            collinear_jump_penalty=9, # 共线跳跃的惩罚
            allowed_chunk_gap=18    # 允许合并的共线段之间的最大查询块间隔
        )
    else:
        print("错误: 找不到 query.txt 或 ref.txt 文件")
        print("请使用 run.py 脚本来分析 data/ 目录中的数据")

def run_analysis_with_files(query_file, ref_file, output_svg, **kwargs):
    """
    使用指定文件运行DNA序列匹配分析。
    
    参数:
        query_file (str): 查询序列文件路径
        ref_file (str): 参考序列文件路径
        output_svg (str): 输出SVG文件路径
        **kwargs: 其他可选参数，包括：
            k (int): K-mer大小，默认6
            chunk_sz (int): 序列块大小，默认25
            kmer_thresh (int): K-mer匹配阈值，默认22
            switching_penalty (int): 切换匹配类型的惩罚，默认10
            collinear_jump_penalty (int): 共线跳跃的惩罚，默认9
            allowed_chunk_gap (int): 允许合并的共线段之间的最大查询块间隔，默认18
    
    返回:
        bool: 分析是否成功完成
    """
    # 设置默认参数
    default_params = {
        'k': 6,
        'chunk_sz': 25,
        'kmer_thresh': 22,
        'switching_penalty': 10,
        'collinear_jump_penalty': 9,
        'allowed_chunk_gap': 18
    }
    
    # 更新参数
    params = default_params.copy()
    params.update(kwargs)
    
    try:
        print(f"\n开始分析:")
        print(f"  查询文件: {query_file}")
        print(f"  参考文件: {ref_file}")
        print(f"  输出SVG: {output_svg}")
        print(f"  参数: {params}")
        
        find_and_visualize_matches(
            query_f=query_file,
            ref_f=ref_file,
            svg_out_f=output_svg,
            **params
        )
        print(f"分析完成，结果保存到: {output_svg}")
        return True
    except Exception as e:
        print(f"分析过程中出现错误: {e}")
        return False

def get_default_parameters():
    """
    获取默认参数设置。
    
    返回:
        dict: 包含所有默认参数的字典
    """
    return {
        'k': 6,                # K-mer 大小
        'chunk_sz': 25,        # 比较的序列块大小
        'kmer_thresh': 22,     # K-mer 匹配阈值 (基于编辑距离的分数)
        'switching_penalty': 10, # 切换匹配类型的惩罚
        'collinear_jump_penalty': 9, # 共线跳跃的惩罚
        'allowed_chunk_gap': 18    # 允许合并的共线段之间的最大查询块间隔
    }

def analyze_sequence_pair(query_file, ref_file, output_dir="result", pair_name="pair"):
    """
    分析一对序列文件并将结果保存到指定目录。
    
    参数:
        query_file (str): 查询序列文件路径
        ref_file (str): 参考序列文件路径
        output_dir (str): 输出目录，默认为"result"
        pair_name (str): 文件对的名称，用于命名输出文件
    
    返回:
        bool: 分析是否成功完成
    """
    import os
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 构造输出文件路径
    output_svg = os.path.join(output_dir, f"{pair_name}_matches_visualization.svg")
    
    # 运行分析
    return run_analysis_with_files(query_file, ref_file, output_svg)

if __name__ == '__main__':
    # 当直接运行此文件时，执行默认分析
    run_default_analysis()



