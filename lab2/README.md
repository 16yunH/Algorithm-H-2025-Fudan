# DNA序列比对算法 - Python实现

高效的生物信息学序列比对工具，提供了模块化的DNA序列比对解决方案。

## 功能特点

- **模块化设计**: 清晰的类结构和接口
- **两种实现**: 内存高效版本和预计算版本
- **完整功能**: 涵盖序列比对的所有关键算法
- **易于使用**: 简单的API和详细的文档

## 项目结构

```
python_implementation/
├── __init__.py              # 包初始化
├── common.py                # 公共数据结构和工具
├── memory_efficient.py     # 内存高效实现
├── precomputed.py          # 预计算实现
├── example.py              # 使用示例和测试
└── README.md               # 说明文档
```

## 核心类

### 1. Segment
数据类，表示匹配片段：
- `x0, y0`: 起始坐标
- `length`: 片段长度
- `direction`: 方向（1=主对角线，0=反对角线）
- `distance`: 错误数

### 2. MemoryEfficientAligner
内存高效的比对器：
- 动态计算网格值
- 适用于大序列
- 内存使用量低

### 3. PrecomputedAligner
预计算的比对器：
- 预建完整点阵图
- 使用累积表快速查询
- 适用于中等序列，查询速度快

## 使用方法

### 基本使用

```python
from python_implementation import MemoryEfficientAligner, PrecomputedAligner

# 准备序列
query = "ATCGATCGATCG"
reference = "ATCGATCGATCG"

# 使用内存高效算法
aligner1 = MemoryEfficientAligner()
results1 = aligner1.align(query, reference)

# 使用预计算算法
aligner2 = PrecomputedAligner()
results2 = aligner2.align(query, reference)

# 查看结果
for segment in results1:
    print(f"匹配片段: {segment}")
```

### 高级使用

```python
from python_implementation import SequenceUtils

# 计算得分
score = SequenceUtils.calculate_score(results1, query, reference)
print(f"总得分: {score}")

# 格式化输出
output = SequenceUtils.format_output(results1)
print(output)
```

## 算法流程

两种实现都遵循相同的核心流程：

1. **对角线片段提取**: 寻找连续的匹配或互补匹配
2. **片段合并**: 合并相近的片段
3. **区间覆盖**: 选择最少片段覆盖序列
4. **空白填充**: 在间隙中寻找新片段
5. **片段扩展**: 向前后扩展片段边界
6. **长度过滤**: 只保留足够长的片段

## 性能对比

| 特性     | MemoryEfficientAligner | PrecomputedAligner |
| -------- | ---------------------- | ------------------ |
| 内存使用 | O(1)                   | O(M×N)             |
| 查询速度 | 慢                     | 快                 |
| 适用场景 | 大序列                 | 中等序列           |
| 构建时间 | 无                     | O(M×N)             |

## 运行示例

```bash
cd python_implementation
python example.py
```

这将运行完整的测试示例，包括：
- 基本功能测试
- 性能比较
- 详细结果分析

## 依赖要求

- Python 3.7+
- NumPy
- 标准库模块

## 安装

```bash
pip install numpy
```

## 设计特点

该实现具有以下特点：

1. **清晰可读**: 优雅的代码结构和命名规范
2. **模块化设计**: 功能分离，便于维护和扩展
3. **类型安全**: 完整的类型注解支持
4. **鲁棒性**: 完善的边界检查和错误处理
5. **易于测试**: 内置示例和性能测试工具

## 算法参数

关键参数可以通过修改各方法的默认值来调整：

- `error_rate`: 错误率阈值 (0.06-0.1)
- `min_length`: 最小片段长度 (5-30)
- `max_gap`: 最大间隙 (1-5)

## 扩展性

该实现易于扩展：
- 可以添加新的匹配策略
- 支持不同的评分函数
- 可以集成到更大的生物信息学流水线中

## 注意事项

1. 对于非常大的序列，建议使用 MemoryEfficientAligner
2. 对于需要频繁查询的应用，PrecomputedAligner 更合适
3. 参数调优可能需要根据具体应用场景进行
