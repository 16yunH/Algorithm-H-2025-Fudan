"""
公共数据结构和工具函数
"""

from dataclasses import dataclass
from typing import List, Tuple
import numpy as np


@dataclass
class Segment:
    """
    表示查询序列和参考序列之间的匹配片段信息
    
    Attributes:
        x0: 片段在查询序列中的起始坐标
        y0: 片段在参考序列中的起始坐标  
        length: 片段的长度
        direction: 方向（1表示主对角线直接匹配，0表示反对角线互补匹配）
        distance: 片段中的错误数（不匹配的位点数）
    """
    x0: int
    y0: int
    length: int
    direction: int  # 1 for main diagonal, 0 for anti-diagonal
    distance: int = 0
    
    @property
    def x1(self) -> int:
        """片段结束的x坐标"""
        return self.x0 + self.length - 1
    
    @property
    def y1(self) -> int:
        """片段结束的y坐标"""
        if self.direction == 1:
            return self.y0 + self.length - 1
        else:
            return self.y0 - (self.length - 1)
    
    @property
    def error_rate(self) -> float:
        """错误率"""
        return self.distance / self.length if self.length > 0 else 0.0
    
    def __repr__(self):
        return f"Segment(x0={self.x0}, y0={self.y0}, length={self.length}, " \
               f"direction={self.direction}, distance={self.distance})"


class SequenceUtils:
    """序列处理工具类"""
    
    @staticmethod
    def is_complement(base_a: str, base_b: str) -> bool:
        """
        检查两个碱基是否互补（A-T 或 C-G）
        
        Args:
            base_a: 第一个碱基
            base_b: 第二个碱基
            
        Returns:
            bool: 如果互补返回True，否则返回False
        """
        complement_pairs = {
            ('A', 'T'), ('T', 'A'),
            ('C', 'G'), ('G', 'C')
        }
        return (base_a, base_b) in complement_pairs
    
    @staticmethod
    def get_match_value(query_base: str, ref_base: str) -> int:
        """
        获取两个碱基的匹配值
        
        Args:
            query_base: 查询序列的碱基
            ref_base: 参考序列的碱基
            
        Returns:
            int: 1表示直接匹配，2表示互补匹配，0表示不匹配
        """
        if query_base == ref_base:
            return 1  # 直接匹配
        elif SequenceUtils.is_complement(query_base, ref_base):
            return 2  # 互补匹配
        else:
            return 0  # 不匹配
    
    @staticmethod
    def extract_runs(line: List[int]) -> Tuple[List[int], List[int]]:
        """
        从对角线提取连续匹配段
        
        Args:
            line: 对角线上的匹配值列表
            
        Returns:
            Tuple[List[int], List[int]]: 匹配段的起点和终点列表
        """
        starts, ends = [], []
        mask = [0] + [1 if val != 0 else 0 for val in line] + [0]
        
        for i in range(1, len(mask) - 1):
            if mask[i] == 1 and mask[i-1] == 0:
                starts.append(i-1)
            if mask[i] == 1 and mask[i+1] == 0:
                ends.append(i)
        
        return starts, ends
    
    @staticmethod
    def format_output(segments: List[Segment], min_length: int = 30) -> str:
        """
        格式化输出匹配结果
        
        Args:
            segments: 匹配片段列表
            min_length: 最小输出长度阈值
            
        Returns:
            str: 格式化的输出字符串
        """
        result = "Remain segments:\n"
        for seg in segments:
            if seg.length >= min_length:
                if seg.direction == 1:  # 主对角线
                    result += f"({seg.x0},{seg.x1 + 1},{seg.y0},{seg.y1 + 1}),"
                else:  # 反对角线
                    result += f"({seg.x0},{seg.x1 + 1},{seg.y1},{seg.y0 + 1}),"
        return result + "\n"
    
    @staticmethod
    def calculate_score(segments: List[Segment], query: str, reference: str) -> int:
        """
        计算最终得分
        
        Args:
            segments: 匹配片段列表
            query: 查询序列
            reference: 参考序列
            
        Returns:
            int: 总得分
        """
        total_score = 0
        for seg in segments:
            if seg.length <= 0:
                continue
                
            # 计算有效匹配数
            effective_matches = 0
            for i in range(seg.length):
                x = seg.x0 + i
                if seg.direction == 1:
                    y = seg.y0 + i
                else:
                    y = seg.y0 - i
                
                if 0 <= x < len(query) and 0 <= y < len(reference):
                    match_val = SequenceUtils.get_match_value(query[x], reference[y])
                    if (seg.direction == 1 and match_val == 1) or \
                       (seg.direction == 0 and match_val == 2):
                        effective_matches += 1
            
            # 只有匹配率超过90%的片段才计入得分
            if effective_matches / seg.length > 0.9:
                total_score += effective_matches
        
        return total_score


class IntervalScheduler:
    """区间调度算法实现"""
    
    @staticmethod
    def weighted_interval_scheduling(segments: List[Segment]) -> List[Segment]:
        """
        加权区间调度算法，选择不重叠的片段以最大化总长度
        
        Args:
            segments: 输入片段列表
            
        Returns:
            List[Segment]: 最优片段集合
        """
        if not segments:
            return []
        
        # 按片段结束位置排序
        sorted_indices = sorted(range(len(segments)), 
                              key=lambda i: segments[i].x1)
        
        # 计算每个片段的前驱（最后一个不重叠的片段）
        predecessors = []
        end_positions = [segments[sorted_indices[i]].x1 for i in range(len(segments))]
        
        for i in range(len(segments)):
            start_pos = segments[sorted_indices[i]].x0
            # 二分查找最后一个结束位置小于当前起始位置的片段
            left, right = 0, i
            pred = -1
            while left <= right:
                mid = (left + right) // 2
                if end_positions[mid] < start_pos:
                    pred = mid
                    left = mid + 1
                else:
                    right = mid - 1
            predecessors.append(pred)
        
        # 动态规划选择最优片段
        dp = [0] * (len(segments) + 1)
        choose = [False] * len(segments)
        
        for j in range(1, len(segments) + 1):
            weight = segments[sorted_indices[j-1]].length
            pred_value = dp[predecessors[j-1] + 1] if predecessors[j-1] >= 0 else 0
            
            if weight + pred_value > dp[j-1]:
                dp[j] = weight + pred_value
                choose[j-1] = True
            else:
                dp[j] = dp[j-1]
        
        # 回溯构建结果
        result = []
        j = len(segments)
        while j > 0:
            if choose[j-1]:
                result.append(segments[sorted_indices[j-1]])
                j = predecessors[j-1] + 1
            else:
                j -= 1
        
        # 按起始位置排序
        result.sort(key=lambda seg: seg.x0)
        return result
