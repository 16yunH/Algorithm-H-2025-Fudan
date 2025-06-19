"""
分数计算模块
包含编辑距离计算和最终分数评估功能
"""

import re
import numpy as np
import edlib
from sequence_utils import reverse_complement

def get_points(tuples_str):
    """
    从表示坐标元组的字符串中提取数字点。
    例如，从 " ( 0, 300, 0, 300 ) , ( 300, 400, 400, 500 ) " 提取数字。
    """
    numbers = re.findall(r"\d+", tuples_str)
    return [int(num) for num in numbers]

def calculate_distance(ref_s, query_s, r_s_1b, r_e_1b, q_s_1b, q_e_1b):
    """
    计算两个序列片段之间的编辑距离（考虑正向和反向互补）。
    坐标为1-based闭区间。
    """
    # 将1-based闭区间坐标转换为0-based开区间用于Python切片
    A = ref_s[r_s_1b - 1 : r_e_1b] 
    a = query_s[q_s_1b - 1 : q_e_1b]
    _a = reverse_complement(a) # 计算查询片段的反向互补序列
    # 返回正向和反向互补匹配中较小的编辑距离
    return min(edlib.align(A, a)['editDistance'], edlib.align(A, _a)['editDistance'])

def calculate_value(tuples_str, ref_s, query_s):  
    """
    根据给定的坐标元组字符串、参考序列和查询序列计算总分。
    分数计算考虑了编辑距离和对齐长度，并对重叠和片段长度不足进行惩罚。
    """
    try:
        points_data = np.array(get_points(tuples_str)) 
        
        if len(points_data) == 0: # 如果没有点，分数为0
            return 0
        if len(points_data) % 4 != 0: # 点的数量必须是4的倍数
            print(f"警告: 坐标点数量 ({len(points_data)}) 无效，不是4的倍数。返回0。")
            return 0

        total_edit_dist = 0 # 总编辑距离
        total_aligned_len = 0 # 总对齐长度
        prev_q_e = 0 # 前一个查询片段的0-based开区间结束位置
        
        points_arr = points_data.reshape((-1, 4)) # 将点数据重塑为 (N, 4) 的数组
        points_arr = points_arr[points_arr[:, 0].argsort()] # 按查询起始位置排序

        for idx, onetuple in enumerate(points_arr):
            # onetuple 是 (q_s, q_e, r_s, r_e)，均为0-based开区间
            q_s, q_e, r_s, r_e = map(int, onetuple)
            
            # 检查重叠：如果前一个查询片段的结束位置大于当前查询片段的起始位置，则存在重叠
            if idx > 0 and prev_q_e > q_s: 
                print(f"在索引 {idx - 1} 和 {idx} 之间的片段检测到重叠:")
                print(f"前一个查询结束位置 (0-based open): {prev_q_e}, 当前查询起始位置 (0-based open): {q_s}")
                return 0 # 重叠则分数为0
            
            prev_q_e = q_e # 更新前一个查询片段的结束位置

            current_edit_dist = 0
            current_aligned_len = 0

            # 确保片段非空
            if q_s < q_e and r_s < r_e:
                # 将0-based开区间转换为1-based闭区间以调用 calculate_distance
                q_s_1b = q_s + 1
                q_e_1b = q_e 
                r_s_1b = r_s + 1
                r_e_1b = r_e

                current_edit_dist = calculate_distance(ref_s, query_s, r_s_1b, r_e_1b, q_s_1b, q_e_1b)
                current_aligned_len = (q_e - q_s) # 0-based开区间的长度为 e - s
            
            total_edit_dist += current_edit_dist
            total_aligned_len += current_aligned_len
            if current_aligned_len > 0: # 只有当片段贡献了对齐长度时才应用惩罚
                total_aligned_len -= 30  # 每个片段的惩罚

        final_score = max(total_aligned_len - total_edit_dist, 0) # 最终分数不能为负
        return final_score
        
    except Exception as e:
        print(f"calculate_value 函数出错: {e}")
        return 0
