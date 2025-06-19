"""
预计算网格的DNA序列比对算法实现
通过预先构建完整的点阵图和累积表实现高速查询
"""

from typing import List, Tuple, Dict
import numpy as np
from collections import defaultdict

from common import Segment, SequenceUtils, IntervalScheduler


class PrecomputedAligner:
    """
    预计算网格的DNA序列比对器
    
    特点：
    - 预先构建完整的点阵图
    - 使用累积表快速计算区间匹配数
    - 适用于需要频繁查询的中等规模序列
    - 内存使用量较高但查询速度快
    """
    
    def __init__(self):
        self.query = ""
        self.reference = ""
        self.M = 0  # 查询序列长度
        self.N = 0  # 参考序列长度
        self.grid = None  # 点阵图
        self.main_diagonal_table = None  # 主对角线累积表
        self.anti_diagonal_table = None  # 反对角线累积表
    
    def build_dotplot(self, query: str, reference: str) -> np.ndarray:
        """
        构建点阵图，比较查询序列和参考序列
        
        Args:
            query: 查询序列
            reference: 参考序列
            
        Returns:
            np.ndarray: M x N 的点阵图，1表示直接匹配，2表示互补匹配，0表示不匹配
        """
        M, N = len(query), len(reference)
        grid = np.zeros((M, N), dtype=np.uint8)
        
        for x in range(M):
            for y in range(N):
                grid[x, y] = SequenceUtils.get_match_value(query[x], reference[y])
        
        return grid
    
    def init_diagonal_tables(self, grid: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        初始化主对角线和反对角线匹配表
        
        Args:
            grid: 点阵图
            
        Returns:
            Tuple[np.ndarray, np.ndarray]: 主对角线表和反对角线表
        """
        M, N = grid.shape
        main_diagonal = np.zeros((M, N), dtype=np.int16)
        anti_diagonal = np.zeros((M, N), dtype=np.int16)
        
        # 主对角线：从左上到右下，记录连续直接匹配的累计长度
        for i in range(M):
            for j in range(N):
                count = 1 if grid[i, j] == 1 else 0
                if i > 0 and j > 0:
                    count += main_diagonal[i-1, j-1]
                main_diagonal[i, j] = count
        
        # 反对角线：从右上到左下，记录连续互补匹配的累计长度
        for i in range(M):
            for j in range(N-1, -1, -1):
                count = 1 if grid[i, j] == 2 else 0
                if i > 0 and j + 1 < N:
                    count += anti_diagonal[i-1, j+1]
                anti_diagonal[i, j] = count
        
        return main_diagonal, anti_diagonal
    
    def find_diagonal_segments(self, grid: np.ndarray, min_length: int = 5) -> List[Segment]:
        """
        寻找对角线上的匹配片段
        
        Args:
            grid: 点阵图
            min_length: 最小片段长度
            
        Returns:
            List[Segment]: 包含所有匹配片段的列表
        """
        M, N = grid.shape
        segments = []
        
        # 主对角线
        for d in range(-M + 1, N):
            start_x = max(0, -d)
            start_y = start_x + d
            
            # 提取对角线
            diagonal = []
            x, y = start_x, start_y
            while x < M and y < N:
                diagonal.append(grid[x, y])
                x += 1
                y += 1
            
            # 提取连续匹配段（直接匹配）
            starts, ends = SequenceUtils.extract_runs([1 if val == 1 else 0 for val in diagonal])
            for start, end in zip(starts, ends):
                length = end - start
                if length >= min_length:
                    x0 = start_x + start
                    y0 = x0 + d
                    segments.append(Segment(x0, y0, length, 1, 0))
        
        # 反对角线（通过水平翻转网格处理）
        flipped_grid = np.fliplr(grid)
        for d in range(-M + 1, N):
            start_x = max(0, -d)
            start_y = N - 1 - (start_x + d)
            
            # 提取反对角线
            diagonal = []
            x, y = start_x, start_y
            while x < M and y >= 0:
                diagonal.append(flipped_grid[x, N-1-y])
                x += 1
                y -= 1
            
            # 提取连续匹配段（互补匹配）
            starts, ends = SequenceUtils.extract_runs([1 if val == 2 else 0 for val in diagonal])
            for start, end in zip(starts, ends):
                length = end - start
                if length >= min_length:
                    x0 = start_x + start
                    y0 = N - 1 - (start_x + start + d)
                    segments.append(Segment(x0, y0, length, 0, 0))
        
        return segments
    
    def compute_diagonal_matches_fast(self, x0: int, y0: int, length: int, direction: int) -> int:
        """
        使用累积表快速计算对角线匹配数
        
        Args:
            x0: 起始x坐标
            y0: 起始y坐标
            length: 长度
            direction: 方向
            
        Returns:
            int: 匹配数
        """
        if length <= 0:
            return 0
        
        x1 = x0 + length - 1
        if direction == 1:
            y1 = y0 + length - 1
            if x1 >= self.M or y1 >= self.N:
                return 0
            
            prev = 0
            if x0 > 0 and y0 > 0:
                prev = self.main_diagonal_table[x0-1, y0-1]
            return self.main_diagonal_table[x1, y1] - prev
        else:
            y1 = y0 - (length - 1)
            if x1 >= self.M or y1 < 0:
                return 0
            
            prev = 0
            if x0 > 0 and y0 + 1 < self.N:
                prev = self.anti_diagonal_table[x0-1, y0+1]
            return self.anti_diagonal_table[x1, y1] - prev
    
    def merge_in_blanks(self, segments: List[Segment], error_rate: float) -> List[Segment]:
        """
        在片段之间的空白区域合并片段
        
        Args:
            segments: 输入片段
            error_rate: 最大错误率
            
        Returns:
            List[Segment]: 合并后的片段集合
        """
        if not segments:
            return []
        
        segments.sort(key=lambda s: s.x0)
        merged = []
        prev = segments[0]
        
        for i in range(1, len(segments)):
            curr = segments[i]
            
            if (curr.direction == prev.direction and
                self._same_diagonal(prev, curr)):
                
                new_length = curr.x1 - prev.x0 + 1
                # 使用快速计算方法
                effective = self.compute_diagonal_matches_fast(
                    prev.x0, prev.y0, new_length, prev.direction)
                distance = new_length - effective
                
                if distance / new_length < error_rate:
                    prev.length = new_length
                    prev.distance = distance
                    continue
            
            merged.append(prev)
            prev = curr
        
        merged.append(prev)
        return merged
    
    def _same_diagonal(self, seg1: Segment, seg2: Segment) -> bool:
        """检查两个片段是否在同一对角线上"""
        if seg1.direction == 1:
            return (seg2.y0 - seg2.x0) == (seg1.y0 - seg1.x0)
        else:
            return (seg1.x0 + seg1.y0) == (seg2.x0 + seg2.y0)
    
    def find_large_segments_in_range(self, x_start: int, x_end: int,
                                   error_rate: float, min_length: int) -> List[Segment]:
        """
        在指定范围内寻找满足错误率要求的大片段
        
        Args:
            x_start: 范围起始
            x_end: 范围结束
            error_rate: 最大错误率
            min_length: 最小长度
            
        Returns:
            List[Segment]: 找到的片段集合
        """
        segments = []
        diagonal_groups = defaultdict(list)
        
        # 主对角线
        for d in range(-self.M + 1, self.N):
            lo = max(0, -d, x_start)
            hi = min(self.M - 1, self.N - 1 - d, x_end)
            
            if hi - lo + 1 >= 1:
                points = [(x, x + d) for x in range(lo, hi + 1)]
                diagonal_groups[(1, d)] = points
        
        # 反对角线
        for d in range(self.M + self.N - 1):
            lo = max(0, d - (self.N - 1), x_start)
            hi = min(self.M - 1, d, x_end)
            
            if hi - lo + 1 >= 1:
                points = [(x, d - x) for x in range(lo, hi + 1)]
                diagonal_groups[(0, d)] = points
        
        # 处理每条对角线
        for (direction, _), points in diagonal_groups.items():
            segments.extend(self._process_diagonal_dp(
                points, direction, error_rate, min_length))
        
        return segments
    
    def _process_diagonal_dp(self, points: List[Tuple[int, int]], direction: int,
                           error_rate: float, min_length: int) -> List[Segment]:
        """使用动态规划处理单条对角线"""
        L = len(points)
        if L == 0:
            return []
        
        # 计算累积错误数
        errors = [0] * (L + 1)
        for i in range(L):
            x, y = points[i]
            if 0 <= x < self.M and 0 <= y < self.N:
                is_error = ((direction == 1 and self.grid[x, y] != 1) or
                           (direction == 0 and self.grid[x, y] != 2))
                errors[i + 1] = errors[i] + (1 if is_error else 0)
            else:
                errors[i + 1] = errors[i] + 1
        
        # 动态规划寻找最优片段
        dp = [0] * (L + 1)
        prev_idx = list(range(L + 1))
        left_ptr = 0
        
        for right in range(1, L + 1):
            while (left_ptr < right and
                   (errors[right] - errors[left_ptr]) > error_rate * (right - left_ptr)):
                left_ptr += 1
            
            dp[right] = dp[right - 1]
            prev_idx[right] = right - 1
            
            if right - left_ptr >= min_length:
                length = right - left_ptr
                if dp[left_ptr] + length > dp[right]:
                    dp[right] = dp[left_ptr] + length
                    prev_idx[right] = left_ptr
        
        # 回溯提取片段
        segments = []
        idx = L
        while idx > 0:
            pi = prev_idx[idx]
            if pi < idx - 1:
                start, end = pi, idx - 1
                x0, y0 = points[start]
                length = end - start + 1
                error_count = errors[end + 1] - errors[start]
                segments.append(Segment(x0, y0, length, direction, error_count))
                idx = pi
            else:
                idx -= 1
        
        return segments
    
    def fill_in_blanks_global(self, segments: List[Segment], error_rate: float,
                            min_gap: int) -> List[Segment]:
        """
        全局填充空白区域
        
        Args:
            segments: 输入片段
            error_rate: 最大错误率
            min_gap: 最小间隙长度
            
        Returns:
            List[Segment]: 填充后的片段集合
        """
        if not segments:
            return []
        
        segments.sort(key=lambda s: s.x0)
        result = []
        
        for i in range(len(segments) - 1):
            result.append(segments[i])
            
            gap_start = segments[i].x1 + 1
            gap_end = segments[i + 1].x0 - 1
            gap_length = gap_end - gap_start + 1
            
            if gap_length >= min_gap:
                extras = self.find_large_segments_in_range(
                    gap_start, gap_end, error_rate, min_gap)
                result.extend(extras)
        
        if segments:
            result.append(segments[-1])
        
        # 合并和优化
        result = self.merge_in_blanks(result, error_rate)
        result = IntervalScheduler.weighted_interval_scheduling(result)
        
        return result
    
    def extend_end_backward(self, segments: List[Segment], error_rate: float) -> List[Segment]:
        """
        向后扩展片段的结束位置
        
        Args:
            segments: 输入片段
            error_rate: 最大错误率
            
        Returns:
            List[Segment]: 扩展后的片段集合
        """
        if not segments:
            return []
        
        # 按x0降序排序
        segments.sort(key=lambda s: s.x0, reverse=True)
        result = []
        prev_end = segments[0].x1 + 1
        
        for seg in segments:
            target_end = prev_end - 1
            space = target_end - seg.x1
            
            new_seg = seg
            if seg.x1 < target_end:
                # 尝试扩展
                best_length = seg.length
                for ext in range(1, space + 1):
                    new_length = seg.length + ext
                    matches = self.compute_diagonal_matches_fast(
                        seg.x0, seg.y0, new_length, seg.direction)
                    distance = new_length - matches
                    
                    if distance / new_length < error_rate:
                        best_length = new_length
                    else:
                        break
                
                if best_length > seg.length:
                    matches = self.compute_diagonal_matches_fast(
                        seg.x0, seg.y0, best_length, seg.direction)
                    new_seg = Segment(seg.x0, seg.y0, best_length, seg.direction,
                                    best_length - matches)
            
            result.append(new_seg)
            prev_end = seg.x0
        
        # 按x0升序重新排序
        result.sort(key=lambda s: s.x0)
        return result
    
    def extend_start_backward(self, segments: List[Segment], error_rate: float) -> List[Segment]:
        """
        向前扩展片段的起始位置
        
        Args:
            segments: 输入片段
            error_rate: 最大错误率
            
        Returns:
            List[Segment]: 扩展后的片段集合
        """
        if not segments:
            return []
        
        segments.sort(key=lambda s: s.x0)
        result = []
        prev_end = -1
        
        for seg in segments:
            target_start = prev_end + 1
            space = seg.x0 - target_start
            
            new_seg = seg
            if seg.x0 > target_start:
                # 尝试扩展
                best_x0 = seg.x0
                best_length = seg.length
                
                for ext in range(1, space + 1):
                    new_x0 = seg.x0 - ext
                    new_length = seg.length + ext
                    step = 1 if seg.direction == 1 else -1
                    new_y0 = seg.y0 - step * ext
                    
                    if new_x0 < 0 or new_y0 < 0 or new_y0 >= self.N:
                        break
                    
                    matches = self.compute_diagonal_matches_fast(
                        new_x0, new_y0, new_length, seg.direction)
                    distance = new_length - matches
                    
                    if distance / new_length < error_rate:
                        best_x0 = new_x0
                        best_length = new_length
                    else:
                        break
                
                if best_x0 < seg.x0:
                    step = 1 if seg.direction == 1 else -1
                    best_y0 = seg.y0 - step * (seg.x0 - best_x0)
                    matches = self.compute_diagonal_matches_fast(
                        best_x0, best_y0, best_length, seg.direction)
                    new_seg = Segment(best_x0, best_y0, best_length, seg.direction,
                                    best_length - matches)
            
            result.append(new_seg)
            prev_end = result[-1].x1
        
        return result
    
    def choose_segments(self, segments: List[Segment], min_length: int) -> List[Segment]:
        """选择长度大于等于min_length的片段"""
        return [seg for seg in segments if seg.length >= min_length]
    
    def minimal_interval_cover(self, segments: List[Segment], error_rate: float,
                             length_thresh: int) -> List[Segment]:
        """
        最小区间覆盖算法，选择覆盖查询序列的最少片段
        
        Args:
            segments: 输入片段
            error_rate: 最大错误率
            length_thresh: 最小长度阈值
            
        Returns:
            List[Segment]: 覆盖片段集合
        """
        if not segments:
            return []
        
        # 确定覆盖范围
        min_x = min(seg.x0 for seg in segments)
        max_x = max(seg.x1 for seg in segments)
        
        # 构建最佳结束位置表
        best_end = [-1] * (max_x + 1)
        best_idx = [-1] * (max_x + 1)
        
        for i, seg in enumerate(segments):
            if seg.x1 > best_end[seg.x0]:
                best_end[seg.x0] = seg.x1
                best_idx[seg.x0] = i
        
        result = []
        covered_end = min_x - 1
        
        while covered_end < max_x:
            start_pos = covered_end + 1
            candidate_end = -1
            candidate_idx = -1
            
            for x0 in range(min_x, start_pos + 1):
                if best_end[x0] > candidate_end:
                    candidate_end = best_end[x0]
                    candidate_idx = best_idx[x0]
            
            if candidate_idx < 0 or candidate_end < start_pos:
                break
            
            best_seg = segments[candidate_idx]
            new_length = min(best_seg.x1, candidate_end) - start_pos + 1
            
            if best_seg.direction == 1:
                new_y0 = best_seg.y0 + (start_pos - best_seg.x0)
            else:
                new_y0 = best_seg.y0 - (start_pos - best_seg.x0)
            
            matches = self.compute_diagonal_matches_fast(
                start_pos, new_y0, new_length, best_seg.direction)
            distance = new_length - matches
            
            if (new_length >= length_thresh and
                distance / new_length < error_rate):
                result.append(Segment(start_pos, new_y0, new_length,
                                    best_seg.direction, distance))
            
            covered_end = candidate_end
        
        return result
    
    def merge_with_tolerance(self, segments: List[Segment], max_gap: int,
                           max_error_rate: float) -> List[Segment]:
        """
        合并具有容差的片段
        
        Args:
            segments: 输入片段
            max_gap: 最大间隙
            max_error_rate: 最大错误率
            
        Returns:
            List[Segment]: 合并后的片段集合
        """
        if not segments:
            return []
        
        # 按对角线分组
        groups = defaultdict(list)
        for i, seg in enumerate(segments):
            if seg.direction == 1:
                diag_id = seg.y0 - seg.x0
            else:
                diag_id = seg.y0 + seg.x0
            groups[(seg.direction, diag_id)].append(i)
        
        merged = []
        for indices in groups.values():
            # 按x0排序
            indices.sort(key=lambda i: segments[i].x0)
            
            current = segments[indices[0]]
            for i in range(1, len(indices)):
                next_seg = segments[indices[i]]
                gap = next_seg.x0 - current.x1 - 1
                merged_length = next_seg.x1 - current.x0 + 1
                
                # 检查对角线对齐
                if current.direction == 1:
                    cond_align = (next_seg.y0 - current.y1 == gap + 1)
                else:
                    cond_align = (current.y1 - next_seg.y0 == gap + 1)
                
                # 检查是否可以合并
                if (gap <= max_gap and cond_align and
                    (gap + current.distance + next_seg.distance) / merged_length <= max_error_rate):
                    current.length = merged_length
                    current.distance = gap + current.distance + next_seg.distance
                else:
                    merged.append(current)
                    current = next_seg
            
            merged.append(current)
        
        return merged
    
    def align(self, query: str, reference: str) -> List[Segment]:
        """
        主要的序列比对方法
        
        Args:
            query: 查询序列
            reference: 参考序列
            
        Returns:
            List[Segment]: 最佳匹配片段集合
        """
        self.query = query
        self.reference = reference
        self.M = len(query)
        self.N = len(reference)
        
        # 构建点阵图和累积表
        self.grid = self.build_dotplot(query, reference)
        self.main_diagonal_table, self.anti_diagonal_table = \
            self.init_diagonal_tables(self.grid)
        
        # 执行一系列处理步骤
        segments = self.find_diagonal_segments(self.grid)
        segments = self.merge_with_tolerance(segments, 1, 0.039)
        segments = self.minimal_interval_cover(segments, 0.08, 20)
        segments = self.fill_in_blanks_global(segments, 0.065, 25)
        segments = self.choose_segments(segments, 25)
        segments = self.extend_start_backward(segments, 0.1)
        
        # 向后扩展需要重新排序
        segments.sort(key=lambda s: s.x0, reverse=True)
        segments = self.extend_end_backward(segments, 0.1)
        
        # 最终排序
        segments.sort(key=lambda s: s.x0)
        return segments
