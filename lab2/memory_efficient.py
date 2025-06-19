"""
内存高效的DNA序列比对算法实现
采用动态计算网格值的策略以节省内存使用
"""

from typing import List, Tuple, Dict, Set
import bisect
from collections import defaultdict

from common import Segment, SequenceUtils, IntervalScheduler


class MemoryEfficientAligner:
    """
    内存高效的DNA序列比对器
    
    特点：
    - 动态计算网格值，不预存储整个点阵图
    - 适用于大序列比对
    - 内存使用量低
    """
    
    def __init__(self):
        self.query = ""
        self.reference = ""
        self.M = 0  # 查询序列长度
        self.N = 0  # 参考序列长度
    
    def get_grid_value(self, x: int, y: int) -> int:
        """
        动态计算网格中指定位置(x, y)的值
        
        Args:
            x: 查询序列坐标
            y: 参考序列坐标
            
        Returns:
            int: 1表示直接匹配，2表示互补匹配，0表示不匹配
        """
        if x >= self.M or y >= self.N or x < 0 or y < 0:
            return 0
        
        return SequenceUtils.get_match_value(self.query[x], self.reference[y])
    
    def compute_diagonal_matches(self, x_start: int, y_start: int, length: int, 
                               direction: int) -> int:
        """
        计算指定对角线范围内的匹配数
        
        Args:
            x_start: 起始x坐标
            y_start: 起始y坐标
            length: 对角线长度
            direction: 方向（1为主对角线，0为反对角线）
            
        Returns:
            int: 匹配点的数量
        """
        matches = 0
        for i in range(length):
            x = x_start + i
            y = y_start + i if direction == 1 else y_start - i
            
            if x >= self.M or y < 0 or y >= self.N:
                break
                
            val = self.get_grid_value(x, y)
            if (direction == 1 and val == 1) or (direction == 0 and val == 2):
                matches += 1
        
        return matches
    
    def find_diagonal_segments(self, min_length: int = 5) -> List[Segment]:
        """
        寻找对角线上的匹配片段，不存储整个网格以节省内存
        
        Args:
            min_length: 最小片段长度
            
        Returns:
            List[Segment]: 匹配片段列表
        """
        segments = []
        
        # 处理主对角线（direction = 1，直接匹配）
        for d in range(-self.M + 1, self.N):
            start_x = max(0, -d)
            start_y = start_x + d
            
            # 收集对角线上的值
            diagonal_values = []
            x, y = start_x, start_y
            while x < self.M and y < self.N:
                diagonal_values.append(self.get_grid_value(x, y))
                x += 1
                y += 1
            
            # 提取连续匹配段
            starts, ends = SequenceUtils.extract_runs(diagonal_values)
            for start, end in zip(starts, ends):
                length = end - start
                if length >= min_length:
                    x0 = start_x + start
                    y0 = x0 + d
                    segments.append(Segment(x0, y0, length, 1, 0))
        
        # 处理反对角线（direction = 0，互补匹配）
        for d in range(-self.M + 1, self.N):
            start_x = max(0, -d)
            start_y = self.N - 1 - (start_x + d)
            
            # 收集反对角线上的值
            diagonal_values = []
            x, y = start_x, start_y
            while x < self.M and y >= 0:
                # 注意：反对角线需要检查互补匹配
                val = self.get_grid_value(x, self.N - 1 - y)
                diagonal_values.append(2 if val == 2 else 0)
                x += 1
                y -= 1
            
            # 提取连续匹配段
            starts, ends = SequenceUtils.extract_runs(diagonal_values)
            for start, end in zip(starts, ends):
                length = end - start
                if length >= min_length:
                    x0 = start_x + start
                    y0 = self.N - 1 - (start_x + start + d)
                    segments.append(Segment(x0, y0, length, 0, 0))
        
        return segments
    
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
        
        # 按x0升序排序
        segments.sort(key=lambda s: s.x0)
        
        merged = []
        prev = segments[0]
        
        for i in range(1, len(segments)):
            curr = segments[i]
            
            # 检查是否可以合并
            if (curr.direction == prev.direction and
                self._same_diagonal(prev, curr)):
                
                new_length = curr.x1 - prev.x0 + 1
                effective_matches = self.compute_diagonal_matches(
                    prev.x0, prev.y0, new_length, prev.direction)
                distance = new_length - effective_matches
                
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
            
            if hi - lo + 1 >= min_length:
                points = [(x, x + d) for x in range(lo, hi + 1)]
                diagonal_groups[(1, d)] = points
        
        # 反对角线
        for d in range(self.M + self.N - 1):
            lo = max(0, d - (self.N - 1), x_start)
            hi = min(self.M - 1, d, x_end)
            
            if hi - lo + 1 >= min_length:
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
        
        # 计算累积错误数
        errors = [0] * (L + 1)
        for i in range(L):
            x, y = points[i]
            val = self.get_grid_value(x, y)
            is_error = (direction == 1 and val != 1) or (direction == 0 and val != 2)
            errors[i + 1] = errors[i] + (1 if is_error else 0)
        
        # 动态规划寻找最优片段
        dp = [0] * (L + 1)
        prev_idx = list(range(L + 1))
        left_ptr = 0
        
        for right in range(1, L + 1):
            # 移动左指针确保错误率
            while left_ptr < right and (errors[right] - errors[left_ptr]) > error_rate * (right - left_ptr):
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
    
    def extend_segments(self, segments: List[Segment], error_rate: float, 
                       extend_start: bool = True) -> List[Segment]:
        """
        扩展片段的起始或结束位置
        
        Args:
            segments: 输入片段
            error_rate: 最大错误率
            extend_start: True扩展起始位置，False扩展结束位置
            
        Returns:
            List[Segment]: 扩展后的片段集合
        """
        if not segments:
            return []
        
        if extend_start:
            segments.sort(key=lambda s: s.x0)
            prev_end = -1
        else:
            segments.sort(key=lambda s: s.x0, reverse=True)
            prev_end = segments[0].x1 + 1
        
        result = []
        
        for seg in segments:
            if extend_start:
                target = prev_end + 1
                space = seg.x0 - target
                if space > 0:
                    # 尝试向前扩展
                    extended = self._try_extend_segment(seg, error_rate, space, True)
                    result.append(extended)
                else:
                    result.append(seg)
                prev_end = result[-1].x1
            else:
                target = prev_end - 1
                space = target - seg.x1
                if space > 0:
                    # 尝试向后扩展
                    extended = self._try_extend_segment(seg, error_rate, space, False)
                    result.append(extended)
                else:
                    result.append(seg)
                prev_end = seg.x0
        
        if not extend_start:
            result.sort(key=lambda s: s.x0)
        
        return result
    
    def _try_extend_segment(self, segment: Segment, error_rate: float, 
                          max_extension: int, extend_start: bool) -> Segment:
        """尝试扩展单个片段"""
        best_segment = segment
        
        for ext in range(1, max_extension + 1):
            if extend_start:
                new_x0 = segment.x0 - ext
                new_length = segment.length + ext
                step = 1 if segment.direction == 1 else -1
                new_y0 = segment.y0 - step * ext
                
                if new_x0 < 0 or new_y0 < 0 or new_y0 >= self.N:
                    break
                    
                matches = self.compute_diagonal_matches(new_x0, new_y0, new_length, segment.direction)
            else:
                new_length = segment.length + ext
                matches = self.compute_diagonal_matches(segment.x0, segment.y0, new_length, segment.direction)
            
            distance = new_length - matches
            if distance / new_length < error_rate:
                if extend_start:
                    best_segment = Segment(new_x0, new_y0, new_length, segment.direction, distance)
                else:
                    best_segment = Segment(segment.x0, segment.y0, new_length, segment.direction, distance)
            else:
                break
        
        return best_segment
    
    def choose_segments(self, segments: List[Segment], min_length: int) -> List[Segment]:
        """选择长度大于等于min_length的片段"""
        return [seg for seg in segments if seg.length >= min_length]
    
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
        
        # 执行一系列处理步骤
        segments = self.find_diagonal_segments(5)
        segments = self._merge_with_tolerance(segments, 3, 0.08)
        segments = self._minimal_interval_cover(segments, 0.085, 30)
        segments = self.fill_in_blanks_global(segments, 0.075, 15)
        segments = self.choose_segments(segments, 30)
        segments = self.extend_segments(segments, 0.085, extend_start=True)
        segments = self.extend_segments(segments, 0.085, extend_start=False)
        
        # 最终排序
        segments.sort(key=lambda s: s.x0)
        return segments
    
    def _merge_with_tolerance(self, segments: List[Segment], max_gap: int, 
                            max_error_rate: float) -> List[Segment]:
        """合并具有容差的片段"""
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
                
                # 检查是否可以合并
                if (gap <= max_gap and 
                    (gap + current.distance + next_seg.distance) / merged_length <= max_error_rate):
                    current.length = merged_length
                    current.distance = gap + current.distance + next_seg.distance
                else:
                    merged.append(current)
                    current = next_seg
            
            merged.append(current)
        
        return merged
    
    def _minimal_interval_cover(self, segments: List[Segment], error_rate: float, 
                              length_thresh: int) -> List[Segment]:
        """最小区间覆盖算法"""
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
            
            matches = self.compute_diagonal_matches(start_pos, new_y0, new_length, best_seg.direction)
            distance = new_length - matches
            
            if new_length >= length_thresh and distance / new_length < error_rate:
                result.append(Segment(start_pos, new_y0, new_length, best_seg.direction, distance))
            
            covered_end = candidate_end
        
        return result
