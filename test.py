class Graph:
    def __init__(self):
        self.edges = {}
    def insert_edge(self, u, v, cost, info):
        self.edges.setdefault(u, {})[v] = {'cost':cost, 'info':info}
    def have_edge(self, u, v):
        return u in self.edges and v in self.edges[u]
    def all_edges(self, u):
        if u not in self.edges:
            return []
        return self.edges[u].keys()
    def get_edge_info(self, u, v):
        return self.edges[u][v]['info']
    def get_edge_cost(self, u, v):
        return self.edges[u][v]['cost']

class Heap:
    def __init__(self, cmp=lambda x, y: x < y):
        self.data = []
        self.size = 0
        self.cmp = cmp
    
    def heapify(self, i):
        while i > 0:
            j = (i - 1) // 2
            if self.cmp(self.data[i][0], self.data[j][0]):
                self.data[i], self.data[j] = self.data[j], self.data[i]
                i = j
            else:
                break
        
        while i < self.size:
            j = 2 * i + 1
            if j >= self.size:
                break
            if j + 1 < self.size and self.cmp(self.data[j + 1][0], self.data[j][0]):
                j += 1
            if self.cmp(self.data[j][0], self.data[i][0]):
                self.data[i], self.data[j] = self.data[j], self.data[i]
                i = j
            else:
                break
    
    def insert(self, value, item):
        self.data.append((value, item))
        self.size += 1
        self.heapify(self.size - 1)
    
    def top(self):
        if self.size == 0:
            return None
        return self.data[0]
    
    def pop(self):
        if self.size == 0:
            return None
        
        ret = self.data[0]
        self.data[0] = self.data[self.size - 1]
        self.data.pop()
        self.size -= 1
        
        if self.size > 0:
            self.heapify(0)
            
        return ret
    
    def get_size(self):  
        return self.size
    
    def is_empty(self):
        return self.size == 0
    
def reverse_complement(dna):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))


def find_occurrences_with_rolling_hash(query, reference, min_length=3, max_length=100):
    """使用滚动哈希查找在reference中出现过（原串或互补串）的query子串
    
    Args:
        query: 查询序列
        reference: 参考序列
        min_length: 最小子串长度
        max_length: 最大子串长度
        
    Returns:
        包含所有出现信息的列表，每项包含(ref_start, qry_start, length, count, is_reverse)
    """
    n = len(query)
    m = len(reference)
    
    results = []
    rev_reference = reverse_complement(reference)
    
    # 限制最大长度，避免过度计算
    max_length = min(max_length, m, n)
    
    # 定义哈希函数参数
    base = 4  # DNA只有4种碱基
    mod = 10**9 + 7  # 大质数，避免哈希碰撞
    
    # DNA碱基到数字的映射
    base_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    # 对参考序列中的子串进行哈希
    for length in range(min_length, max_length + 1):
        # 预计算base^(length-1) % mod
        highest_power = pow(base, length-1, mod)
        
        # 参考序列哈希表
        ref_hash_table = {}
        
        # 计算参考序列第一个窗口的哈希值
        ref_hash = 0
        rev_hash = 0
        for i in range(length):
            if i < m and reference[i] in base_to_num:
                ref_hash = (ref_hash * base + base_to_num[reference[i]]) % mod
            if i < m and rev_reference[i] in base_to_num:
                rev_hash = (rev_hash * base + base_to_num[rev_reference[i]]) % mod
        
        # 存储第一个窗口的哈希值
        ref_hash_table[ref_hash] = [(0, reference[0:length], False)]
        ref_hash_table.setdefault(rev_hash, []).append((0, rev_reference[0:length], True))
        
        # 滚动计算参考序列其余窗口的哈希值
        for i in range(1, m - length + 1):
            # 移除前一个字符
            if reference[i-1] in base_to_num:
                ref_hash = (ref_hash - base_to_num[reference[i-1]] * highest_power) % mod
            # 添加新字符
            if reference[i+length-1] in base_to_num:
                ref_hash = (ref_hash * base + base_to_num[reference[i+length-1]]) % mod
            # 确保哈希值为正
            ref_hash = (ref_hash + mod) % mod
                
            # 同样处理反向互补序列
            if rev_reference[i-1] in base_to_num:
                rev_hash = (rev_hash - base_to_num[rev_reference[i-1]] * highest_power) % mod
            if rev_reference[i+length-1] in base_to_num:
                rev_hash = (rev_hash * base + base_to_num[rev_reference[i+length-1]]) % mod
            rev_hash = (rev_hash + mod) % mod
            
            # 存储哈希值对应的位置和子串
            ref_segment = reference[i:i+length]
            rev_segment = rev_reference[i:i+length]
            ref_hash_table.setdefault(ref_hash, []).append((i, ref_segment, False))
            ref_hash_table.setdefault(rev_hash, []).append((i, rev_segment, True))
        
        # 计算查询序列的哈希值并查找匹配
        for qry_start in range(n - length + 1):
            # 为每个起始位置单独计算哈希值，确保准确性
            qry_hash = 0
            qry_segment = query[qry_start:qry_start + length]
            
            for i in range(length):
                if qry_segment[i] in base_to_num:
                    qry_hash = (qry_hash * base + base_to_num[qry_segment[i]]) % mod
            
            # 检查是否有匹配
            if qry_hash in ref_hash_table:
                for ref_start, ref_segment, is_reverse in ref_hash_table[qry_hash]:
                    if qry_segment == ref_segment:  # 哈希值匹配后进行确认比较
                        count = 1
                        next_pos = qry_start + length
                        while next_pos + length <= n and query[next_pos:next_pos + length] == qry_segment:
                            count += 1
                            next_pos += length
                        results.append((ref_start, qry_start, length, count, is_reverse))
    
    return results
       

def find_occurrences(query, reference):
    """原始方法查找在reference中出现过的query子串（保留以便比较）"""
    n = len(query)
    m = len(reference)
    
    results = []  
    ref_segments = {}  
    rev_reference = reverse_complement(reference)

    # 限制序列长度，避免超时
    max_length = min(100, m//10, n//10)
    min_length = 3
    
    for length in range(min_length, max_length + 1):  
        for ref_start in range(m - length + 1):
            ref_segment = reference[ref_start:ref_start + length]
            rev_segment = rev_reference[m - (ref_start + length):m - ref_start]
            ref_segments.setdefault(ref_segment, []).append((ref_start, length, False))
            ref_segments.setdefault(rev_segment, []).append((ref_start, length, True))

    for length in range(min_length, max_length + 1):  
        for qry_start in range(n - length + 1):
            qry_segment = query[qry_start:qry_start + length]
            if qry_segment in ref_segments:
                for (ref_start, ref_len, is_reverse) in ref_segments[qry_segment]:
                    count = 1
                    next_pos = qry_start + length
                    while next_pos + length <= n and query[next_pos:next_pos + length] == qry_segment:
                        count += 1
                        next_pos += length
                    results.append((ref_start, qry_start, length, count, is_reverse))
    
    return results

def build_occur_graph(repeat_table):
    graph = Graph()
    for ref_start, qry_start, length, count, is_reverse in repeat_table:
        qry_end = qry_start + length * count
        graph.insert_edge(qry_start, qry_end, count, {'ref_start':ref_start, 'length':length, 'count':count, 'is_reverse':is_reverse})
    return graph

def get_current_move(query, reference, graph, repeat_table, ref_pos, qry_pos):
    all_repeat_edges = graph.all_edges(qry_pos)
    if all_repeat_edges == []:
        if qry_pos < len(query) and ref_pos < len(reference) and query[qry_pos] == reference[ref_pos]:
            return [(ref_pos + 1, qry_pos + 1, 0, 0)]
        else: 
            print(f"Failed on {qry_pos}, {ref_pos}")
            return None # Fail to complete the alignment
    else:
        result = []
        if qry_pos < len(query) and ref_pos < len(reference) and query[qry_pos] == reference[ref_pos]:
            result.append((ref_pos + 1, qry_pos + 1, 0, 0))
        for repeat_edge in all_repeat_edges:
            cost = graph.get_edge_cost(qry_pos, repeat_edge)
            info = graph.get_edge_info(qry_pos, repeat_edge)
            curr_ref_pos = ref_pos
            if info['ref_start'] == curr_ref_pos:
                curr_ref_pos += info['length']
                cost -= 1   # 第一次相同的情况, 重复次数减一
            reverse = info['is_reverse']
            result.append((curr_ref_pos, qry_pos + info['length'] * info['count'], cost, reverse))
    return result

def get_best_path(query, reference):
    # 使用滚动哈希算法查找重复序列
    repeat_results = find_occurrences_with_rolling_hash(query, reference, min_length=3, max_length=100)
    
    repeat_table = {}
    for repeat_result in repeat_results:
        ref_start, qry_start, length, count, is_reverse = repeat_result
        repeat_table.setdefault(qry_start, {}).setdefault(ref_start, []).append({'ref_start':ref_start, 'length':length, 'count':count, 'is_reverse':is_reverse})
    graph = build_occur_graph(repeat_results)
    n = len(query)
    minHeap = Heap()
    dist = [float('inf') for i in range(n + 1)]
    dist[0] = 0
    visited = set()
    all_moves = get_current_move(query, reference, graph, repeat_table, 0, 0)
    pre_state = {}
    minHeap.insert(0, (0, 0, 0, 0))
    while not minHeap.is_empty():
        top = minHeap.pop()
        cost = top[0]
        item = top[1]
        ref_pos, qry_pos, cost, reverse = item
        if (ref_pos, qry_pos) in visited:
            continue
        visited.add((ref_pos, qry_pos))
        all_moves = get_current_move(query, reference, graph, repeat_table, ref_pos, qry_pos)
        if all_moves is None:
            continue
        for move in all_moves:
            next_ref_pos, next_qry_pos, next_cost, next_reverse = move
            if next_qry_pos >= n:
                print(f"最优路径代价:{dist[qry_pos]}")
                path = []
                curr_state = (ref_pos, qry_pos)
                while curr_state != (0,0):
                    path.append(curr_state)
                    curr_state = pre_state[curr_state]
                path.append((0,0))  # 确保起点也包含在路径中
                print(f"路径长度: {len(path)}")
                print_ans(path[::-1], query, reference, repeat_table, graph)
                return path
            if dist[next_qry_pos] > dist[qry_pos] + next_cost:
                dist[next_qry_pos] = dist[qry_pos] + next_cost
                minHeap.insert(dist[next_qry_pos], (next_ref_pos, next_qry_pos, next_cost, next_reverse))
                pre_state[(next_ref_pos, next_qry_pos)] = (ref_pos, qry_pos)
    print("Fail to complete the alignment")


def print_repeat_table(repeat_table):
    for qry_start in repeat_table:
        print(f'qry_start:{qry_start}')
        for ref_start in repeat_table[qry_start]:
            print(f'    ref_start:{ref_start}')
            for item in repeat_table[qry_start][ref_start]:
                print(f'        length:{item["length"]}\tcount:{item["count"]}\tis_reverse:{item["is_reverse"]}')
            
def print_results(results):
    if not results:
        print("未检测到重复变异")
        return
        
    print("序号\tref位置\tqry位置\t长度\t次数\t反向重复")
    for idx, (ref_pos, qry_pos, length, count, is_reverse) in enumerate(results, 1):
        print(f"{idx}\t{ref_pos}\t{qry_pos}\t{length}\t{count}\t{is_reverse}")

def print_ans(path, query, reference, repeat_table, graph):
    """
    从最优路径中提取并打印重复变异的详细信息
    
    参数:
        path: 最优路径的节点列表，格式为 [(ref_pos, qry_pos), ...]
        query: 查询序列
        reference: 参考序列
        repeat_table: 包含重复信息的表格
        graph: 有向图对象
    """
    print("\n最优路径中的重复变异:")
    print("序号\tref位置\tqry位置\t长度\t次数\t重复序列\t反向互补")
    
    repeat_count = 0
    for i in range(len(path) - 1):
        curr_ref, curr_qry = path[i]
        next_ref, next_qry = path[i + 1]
        
        if next_qry - curr_qry > 1:
            if graph.have_edge(curr_qry, next_qry):
                repeat_count += 1
                
                info = graph.get_edge_info(curr_qry, next_qry)
                ref_start = info['ref_start']
                length = info['length']
                count = info['count']
                is_reverse = info['is_reverse']
                
                repeat_seq = query[curr_qry:curr_qry + length]
                
                actual_count = (next_qry - curr_qry) // length
                
                if ref_start == curr_ref:
                    actual_count -= 1  # 减去对齐部分
                
                if actual_count > 0:
                    print(f"{repeat_count}\t{ref_start}\t{curr_qry}\t{length}\t{actual_count}\t{repeat_seq}\t{is_reverse}")
            else:
                print(f"从位置 {curr_qry} 到 {next_qry} 的对齐")
        else:
            pass
    
    return

def read_sequence(filename):
    """从文件中读取序列数据"""
    with open(filename, 'r') as f:
        return f.read().strip()

def main():
    # 可以在此处选择使用测试数据还是实际数据
    use_test_data = False
    
    if use_test_data:
        query = read_sequence("query_test.txt")
        reference = read_sequence("ref_test.txt")
    else:
        query = read_sequence("data/query1.txt")
        reference = read_sequence("data/reference1.txt")
        
    print(f"查询序列长度: {len(query)}")
    print(f"参考序列长度: {len(reference)}")
    
    get_best_path(query, reference)

if __name__ == "__main__":
    main()