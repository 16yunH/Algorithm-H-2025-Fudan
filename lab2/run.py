import os
import sys
import time
from pathlib import Path

# 导入main.py中的主要函数
from main import find_and_visualize_matches

def read_sequence_file(filepath):
    """
    读取序列文件
    
    Args:
        filepath (str): 序列文件路径
        
    Returns:
        str: 序列内容（去除换行符和空格）
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            sequence = f.read().replace('\n', '').replace(' ', '').strip().upper()
        return sequence
    except FileNotFoundError:
        print(f"错误: 文件 {filepath} 不存在")
        return None
    except Exception as e:
        print(f"读取文件 {filepath} 时出错: {e}")
        return None

def run_dna_analysis(query_file, reference_file, output_name, **params):
    """
    运行DNA序列分析
    
    Args:
        query_file (str): 查询序列文件路径
        reference_file (str): 参考序列文件路径
        output_name (str): 输出文件名前缀
        **params: 传递给find_and_visualize_matches的参数
    """
    print(f"\n{'='*60}")
    print(f"开始分析: {output_name}")
    print(f"查询序列: {query_file}")
    print(f"参考序列: {reference_file}")
    print(f"{'='*60}")
    
    # 检查文件是否存在
    if not os.path.exists(query_file):
        print(f"错误: 查询序列文件 {query_file} 不存在")
        return False
        
    if not os.path.exists(reference_file):
        print(f"错误: 参考序列文件 {reference_file} 不存在")
        return False
    
    # 获取脚本目录，生成输出文件路径
    script_dir = os.path.dirname(os.path.abspath(__file__))
    result_dir = os.path.join(script_dir, "result")
    svg_output = os.path.join(result_dir, f"{output_name}_matches_visualization.svg")
    
    # 确保结果目录存在
    os.makedirs(result_dir, exist_ok=True)
    
    # 记录开始时间
    start_time = time.time()
    
    try:
        # 调用主分析函数
        find_and_visualize_matches(
            query_f=query_file,
            ref_f=reference_file,
            svg_out_f=svg_output,
            **params
        )
        
        # 计算运行时间
        end_time = time.time()
        runtime = end_time - start_time
        
        print(f"\n分析完成!")
        print(f"运行时间: {runtime:.2f} 秒")
        print(f"SVG可视化结果已保存到: {svg_output}")
        return True
        
    except Exception as e:
        print(f"分析过程中出错: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """
    主函数：运行所有数据对的分析
    """
    print("DNA序列匹配分析系统")
    print("=" * 60)
    
    # 获取脚本所在目录，确保路径正确
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, "data")
    
    print(f"当前工作目录: {os.getcwd()}")
    print(f"脚本目录: {script_dir}")
    print(f"数据目录: {data_dir}")
    
    # 检查数据目录是否存在
    if not os.path.exists(data_dir):
        print(f"错误: 数据目录 {data_dir} 不存在")
        return    
    # 定义要分析的数据对
    data_pairs = [
        {
            "query": os.path.join(data_dir, "query1.txt"),
            "reference": os.path.join(data_dir, "reference1.txt"),
            "name": "pair1"
        },
        {
            "query": os.path.join(data_dir, "query2.txt"),
            "reference": os.path.join(data_dir, "reference2.txt"),
            "name": "pair2"
        }
    ]
    
    # 检查所有数据文件是否存在
    missing_files = []
    for pair in data_pairs:
        if not os.path.exists(pair["query"]):
            missing_files.append(pair["query"])
        if not os.path.exists(pair["reference"]):
            missing_files.append(pair["reference"])
    
    if missing_files:
        print("错误: 以下文件不存在:")
        for file in missing_files:
            print(f"  {file}")
        print("\n请检查data目录中的文件名是否正确")
        return
    
    # 定义分析参数 - 可以根据需要调整
    analysis_params = {
        "k": 6,                          # K-mer 大小
        "chunk_sz": 25,                  # 比较的序列块大小
        "kmer_thresh": 22,               # K-mer 匹配阈值
        "switching_penalty": 10,         # 切换匹配类型的惩罚
        "collinear_jump_penalty": 9,     # 共线跳跃的惩罚
        "allowed_chunk_gap": 18          # 允许合并的共线段之间的最大查询块间隔
    }
    
    print("分析参数:")
    for param, value in analysis_params.items():
        print(f"  {param}: {value}")
    
    # 运行分析
    success_count = 0
    total_start_time = time.time()
    
    for i, pair in enumerate(data_pairs, 1):
        print(f"\n处理数据对 {i}/{len(data_pairs)}: {pair['name']}")
        
        success = run_dna_analysis(
            query_file=pair["query"],
            reference_file=pair["reference"],
            output_name=pair["name"],
            **analysis_params
        )
        
        if success:
            success_count += 1
    
    # 总结
    total_end_time = time.time()
    total_runtime = total_end_time - total_start_time
    
    print(f"\n{'='*60}")
    print("分析总结:")
    print(f"成功分析: {success_count}/{len(data_pairs)} 个数据对")
    print(f"总运行时间: {total_runtime:.2f} 秒")
    print(f"结果保存在 result/ 目录下")
    print(f"{'='*60}")

def run_custom_analysis():
    """
    运行自定义参数的分析
    """
    print("\n自定义参数分析模式")
    print("-" * 40)
    
    # 用户可以在这里自定义参数
    custom_params = {
        "k": 7,                          # 可以尝试不同的k值
        "chunk_sz": 50,                  # 可以尝试不同的块大小
        "kmer_thresh": 45,               # 可以调整阈值
        "switching_penalty": 1000,       # 可以调整惩罚参数
        "collinear_jump_penalty": 10,
        "allowed_chunk_gap": 0
    }
    
    print("自定义参数:")
    for param, value in custom_params.items():
        print(f"  {param}: {value}")
    
    # 获取脚本目录，构建文件路径
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, "data")
    
    # 只分析第一对数据作为示例
    query_file = os.path.join(data_dir, "query1.txt")
    reference_file = os.path.join(data_dir, "reference1.txt")
    
    success = run_dna_analysis(
        query_file=query_file,
        reference_file=reference_file,
        output_name="custom_analysis",
        **custom_params
    )
    
    return success

if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) > 1:
        if sys.argv[1] == "--custom":
            # 运行自定义分析
            run_custom_analysis()
        elif sys.argv[1] == "--help":
            print("用法:")
            print("  python run.py          # 运行标准分析")
            print("  python run.py --custom # 运行自定义参数分析")
            print("  python run.py --help   # 显示帮助信息")
        else:
            print(f"未知参数: {sys.argv[1]}")
            print("使用 --help 查看用法")
    else:
        # 运行标准分析
        main()
