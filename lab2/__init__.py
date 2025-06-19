"""
DNA序列比对算法Python实现
作者: Algorithm H Lab
日期: 2025年春季

本包提供了高效的DNA序列比对算法实现：
1. MemoryEfficientAligner: 内存高效的实现，适用于大序列
2. PrecomputedAligner: 预计算网格的实现，适用于频繁查询
"""

from memory_efficient import MemoryEfficientAligner
from precomputed import PrecomputedAligner
from common import Segment, SequenceUtils

__version__ = "1.0.0"
__author__ = "Algorithm H Lab"

__all__ = [
    'MemoryEfficientAligner',
    'PrecomputedAligner', 
    'Segment',
    'SequenceUtils'
]
