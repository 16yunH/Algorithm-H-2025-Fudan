# Lab 2 实验报告：DNA序列匹配与动态规划

## 1. 实验概述

本实验实现了一个基于动态规划的DNA序列匹配算法，用于在参考序列中寻找查询序列的最佳对齐路径。该算法支持正向匹配和反向互补匹配，并通过动态规划优化找到全局最优的匹配路径。

## 2. 核心算法设计

### 2.1 序列块匹配算法

#### 算法伪代码

```pseudocode
function ProcessQueryChunk(q_idx, query_s, ref_s, chunk_sz, k, kmer_thresh):
    q_start = q_idx * chunk_sz
    q_end = min((q_idx + 1) * chunk_sz, len(query_s))
    q_chunk = query_s[q_start:q_end]
    
    fwd_matches = []
    rc_matches = []
    
    // 正向匹配
    for r_start = 0 to len(ref_s) - chunk_sz:
        ref_chunk = ref_s[r_start:r_start + chunk_sz]
        edit_distance = edlib.align(q_chunk, ref_chunk)
        score = chunk_sz - edit_distance
        if score >= kmer_thresh:
            fwd_matches.append((q_idx, r_start, score))
    
    // 反向互补匹配
    rc_q_chunk = reverse_complement(q_chunk)
    for r_start = 0 to len(ref_s) - chunk_sz:
        ref_chunk = ref_s[r_start:r_start + chunk_sz]
        edit_distance = edlib.align(rc_q_chunk, ref_chunk)
        score = chunk_sz - edit_distance
        if score >= kmer_thresh:
            rc_matches.append((q_idx, r_start, score))
    
    return fwd_matches, rc_matches
```

### 2.2 动态规划最长路径算法

#### 动态规划伪代码

```pseudocode
function CalculateLongestPath(nodes, switching_penalty, collinear_jump_penalty):
    // 按查询位置、参考位置、匹配类型排序
    sort(nodes by (q, r, type))
    
    N = len(nodes)
    dp[N] = [node.original_score for node in nodes]
    predecessor[N] = [-1] * N
    
    // 动态规划
    for i = 0 to N-1:
        for j = 0 to i-1:
            if nodes[i].q > nodes[j].q:  // 有效转换
                penalty = CalculatePenalty(nodes[j], nodes[i])
                potential_score = dp[j] + nodes[i].original_score - penalty
                
                if potential_score > dp[i]:
                    dp[i] = potential_score
                    predecessor[i] = j
    
    // 回溯最优路径
    max_idx = argmax(dp)
    path = BacktrackPath(predecessor, max_idx)
    
    return path, dp[max_idx]

function CalculatePenalty(node_j, node_i):
    if node_i.type != node_j.type:
        return switching_penalty
    
    // 检查是否共线
    delta_q = (node_i.q - node_j.q) * chunk_sz
    
    if node_i.type == "forward":
        expected_r = node_j.r + delta_q
    else:  // reverse complement
        expected_r = node_j.r - delta_q
    
    is_collinear = (node_i.r == expected_r)
    is_adjacent = (node_i.q == node_j.q + 1)
    
    if is_adjacent and not is_collinear:
        return switching_penalty
    elif not is_adjacent and is_collinear:
        return collinear_jump_penalty
    elif not is_adjacent and not is_collinear:
        return switching_penalty
    else:
        return 0
```

### 2.3 路径合并算法

#### 路径合并伪代码

```pseudocode
function MergePathSegments(path_nodes, allowed_gap):
    merged_segments = []
    current_segment = null
    
    for node in path_nodes:
        if current_segment is null:
            current_segment = CreateSegment(node)
        else:
            if CanMerge(current_segment, node, allowed_gap):
                MergeIntoSegment(current_segment, node)
            else:
                merged_segments.append(current_segment)
                current_segment = CreateSegment(node)
    
    if current_segment is not null:
        merged_segments.append(current_segment)
    
    // 扩展最后一个片段到序列末尾
    if merged_segments:
        ExtendLastSegment(merged_segments[-1])
    
    return merged_segments
```

## 3. 时空复杂度分析

### 3.1 时间复杂度

1. **序列块匹配阶段**：
   - 查询块数量：O(m/c)，其中m为查询序列长度，c为块大小
   - 每个查询块对参考序列的比较：O(n-c+1)，其中n为参考序列长度
   - 编辑距离计算：O(c²)
   - **总时间复杂度**：O((m/c) × n × c²) = O(m × n × c)

2. **动态规划阶段**：
   - 匹配节点数量：O(k)，其中k为找到的匹配数
   - 双重循环计算最优路径：O(k²)
   - **总时间复杂度**：O(k²)

3. **路径合并阶段**：
   - 线性扫描路径节点：O(p)，其中p为最优路径长度
   - **总时间复杂度**：O(p)

**算法整体时间复杂度**：O(m × n × c + k²)

### 3.2 空间复杂度

1. **匹配存储**：O(k)存储所有匹配节点
2. **动态规划表**：O(k)存储dp数组和前驱数组
3. **路径存储**：O(p)存储最优路径
4. **序列存储**：O(m + n)存储输入序列

**算法整体空间复杂度**：O(m + n + k)

## 4. 算法参数说明

| 参数                   | 默认值 | 说明                                 |
| ---------------------- | ------ | ------------------------------------ |
| k                      | 6      | K-mer大小                            |
| chunk_sz               | 25     | 序列块大小                           |
| kmer_thresh            | 22     | K-mer匹配阈值（基于编辑距离的分数）  |
| switching_penalty      | 10     | 匹配类型切换惩罚                     |
| collinear_jump_penalty | 9      | 共线跳跃惩罚                         |
| allowed_chunk_gap      | 18     | 允许合并的共线段之间的最大查询块间隔 |

## 5. 实验结果

### 5.1 数据集1 (pair1)

- **查询序列长度**：29,845 bp
- **参考序列长度**：29,830 bp
- **运行时间**：30.90 秒
- **找到匹配**：
  - 正向匹配：2,233 个
  - 反向互补匹配：3,196 个
- **最长对齐路径分数**：29,787.0
- **计算分值**：29,701

**最优对齐结果**：

```text
[ ( 0, 6800, 0, 6800 ) ,   
  ( 6800, 12175, 17654, 23029 ) ,   
  ( 12175, 23725, 6103, 17653 ) ,   
  ( 23725, 29829, 23726, 29830 ) ]
```

**结果分析**：

- 找到4个主要对齐片段
- 第一个片段是完全正向匹配（0-6800）
- 中间两个片段表现出反向互补匹配的特征
- 最后一个片段重新回到正向匹配

### 5.2 数据集2 (pair2)

- **查询序列长度**：2,500 bp
- **参考序列长度**：1,700 bp
- **运行时间**：0.75 秒
- **找到匹配**：
  - 正向匹配：271 个
  - 反向互补匹配：46 个
- **最长对齐路径分数**：1,947.0
- **计算分值**：1,597

**最优对齐结果**：

```text
[ ( 0, 300, 0, 300 ) ,   ( 300, 400, 400, 500 ) ,   ( 400, 500, 500, 600 ) ,   
  ( 500, 700, 600, 800 ) ,   ( 700, 800, 700, 800 ) ,   ( 800, 900, 700, 800 ) ,   
  ( 900, 1000, 700, 800 ) ,   ( 1000, 1200, 700, 900 ) ,   ( 1200, 1300, 900, 1000 ) ,   
  ( 1300, 1400, 900, 1000 ) ,   ( 1400, 1500, 400, 500 ) ,   ( 1500, 1600, 1000, 1100 ) ,   
  ( 1600, 1700, 1300, 1400 ) ,   ( 1700, 1800, 1200, 1300 ) ,   ( 1800, 1900, 1100, 1200 ) ,   
  ( 1900, 2000, 1400, 1500 ) ,   ( 2300, 2500, 1500, 1700 ) ]
```

**结果分析**：

- 找到17个匹配片段，显示出复杂的重排模式
- 查询序列的后段（2000-2300）在参考序列中没有找到匹配
- 存在多个短片段的重叠和跳跃匹配

## 6. 可视化结果

实验生成了SVG格式的可视化图像：

- `pair1_matches_visualization.svg`：展示数据集1的匹配模式
- `pair2_matches_visualization.svg`：展示数据集2的匹配模式

可视化图像包含：

- 正向匹配（蓝色线段和点）
- 反向互补匹配（红色线段和点）
- 最优路径高亮显示
- 网格坐标系便于定位

## 7. 算法优化策略

### 7.1 并行化处理

- 使用多进程并行处理查询块匹配
- 自动回退到单线程处理以保证兼容性

### 7.2 内存优化

- 使用生成器和迭代器减少内存占用
- 及时释放不再需要的中间结果

### 7.3 剪枝策略

- 设置K-mer阈值过滤低质量匹配
- 使用惩罚机制避免不合理的路径跳跃

## 8. 实验总结

本实验成功实现了一个高效的DNA序列匹配算法，主要贡献包括：

1. **算法设计**：结合了局部匹配和全局优化的动态规划方法
2. **性能优化**：通过并行处理和合理的数据结构选择提高了效率
3. **实用性**：支持正向和反向互补匹配，适用于真实的生物序列分析场景
4. **可视化**：提供直观的匹配结果展示

实验结果表明，该算法能够有效处理不同规模的DNA序列匹配问题，在保证准确性的同时具有良好的计算效率。对于大规模序列（如29K bp），算法能在合理时间内完成计算；对于小规模序列，响应时间非常快速。
