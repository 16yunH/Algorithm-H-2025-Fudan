import heapq
from utils import reverse_complement
from graph import Graph

def search_occurrences(query_seq, ref_seq):
    q_len = len(query_seq)
    r_len = len(ref_seq)
    
    occurrences = []  
    seg_dict = {}  
    rev_ref = reverse_complement(ref_seq)

    for seg_length in range(1, r_len + 1):
        for r_start in range(r_len - seg_length + 1):
            seg = ref_seq[r_start:r_start + seg_length]
            rev_seg = rev_ref[r_len - (r_start + seg_length):r_len - r_start]
            seg_dict.setdefault(seg, []).append((r_start, seg_length, False))
            seg_dict.setdefault(rev_seg, []).append((r_start, seg_length, True))

    for seg_length in range(1, q_len + 1):
        for q_start in range(q_len - seg_length + 1):
            sub_query = query_seq[q_start:q_start + seg_length]
            if sub_query in seg_dict:
                for (r_start, r_seg_length, is_rev) in seg_dict[sub_query]:
                    count = 1
                    pos = q_start + seg_length
                    while pos + seg_length <= q_len and query_seq[pos:pos + seg_length] == sub_query:
                        count += 1
                        pos += seg_length
                    occurrences.append((r_start, q_start, seg_length, count, is_rev))
    
    return occurrences

def build_graph_from_occurrences(occurrence_list):
    graph = Graph()
    for r_start, q_start, seg_length, count, is_rev in occurrence_list:
        q_end = q_start + seg_length * count
        graph.add_edge(q_start, q_end, count, {'ref_start': r_start, 'seg_length': seg_length, 'count': count, 'is_rev': is_rev})
    return graph

def get_next_shifts(query_seq, ref_seq, graph, curr_r, curr_q):
    neighbors = graph.neighbors(curr_q)
    shifts = []
    if not neighbors:
        if curr_q < len(query_seq) and curr_r < len(ref_seq) and query_seq[curr_q] == ref_seq[curr_r]:
            shifts.append((curr_r + 1, curr_q + 1, 0, 0))
        else: 
            print(f"Alignment failed at query position {curr_q}, reference position {curr_r}.")
            return None
    else:
        if curr_q < len(query_seq) and curr_r < len(ref_seq) and query_seq[curr_q] == ref_seq[curr_r]:
            shifts.append((curr_r + 1, curr_q + 1, 0, 0))
        for dest in neighbors:
            weight = graph.get_weight(curr_q, dest)
            meta = graph.get_meta(curr_q, dest)
            new_r = curr_r
            if meta['ref_start'] == curr_r:
                new_r += meta['seg_length']
                weight -= 1
            shifts.append((new_r, dest, weight, meta['is_rev']))
    return shifts

def show_alignment_details(path, query_seq, graph):
    print("\nRepeated mutations in the optimal path:")
    header = f"{'Index':<6} {'RefStart':<10} {'QueryStart':<12} {'SegLength':<10} {'Repeats':<8} {'Segment':<15} {'Reverse':<10}"
    print(header)
    rep_index = 0
    for i in range(len(path) - 1):
        curr_r, curr_q = path[i]
        next_r, next_q = path[i + 1]
        
        if next_q - curr_q > 1:
            if graph.edge_exists(curr_q, next_q):
                rep_index += 1
                meta = graph.get_meta(curr_q, next_q)
                seg = query_seq[curr_q:curr_q + meta['seg_length']]
                actual_repeat = (next_q - curr_q) // meta['seg_length']
                if meta['ref_start'] == curr_r:
                    actual_repeat -= 1
                if actual_repeat > 0:
                    print(f"{rep_index:<6} {meta['ref_start']:<10} {curr_q:<12} {meta['seg_length']:<10} {actual_repeat:<8} {seg:<15} {meta['is_rev']!s:<10}")
            else:
                print(f"Alignment from query pos {curr_q} to {next_q}.")
    return

def print_occurrence_list(occurrences):
    if not occurrences:
        print("No repeat mutations detected")
        return
    print(f"{'Index':<6} {'RefStart':<10} {'QueryStart':<12} {'SegLength':<10} {'RepeatCount':<12} {'Reverse':<10}")
    for i, (r_start, q_start, seg_length, count, is_rev) in enumerate(occurrences, 1):
        print(f"{i:<6} {r_start:<10} {q_start:<12} {seg_length:<10} {count:<12} {is_rev!s:<10}")

def optimal_alignment(query_seq, ref_seq):
    occ_list = search_occurrences(query_seq, ref_seq)
    occ_table = {}
    for occ in occ_list:
        r_start, q_start, seg_length, count, is_rev = occ
        occ_table.setdefault(q_start, {}).setdefault(r_start, []).append({
            'ref_start': r_start,
            'seg_length': seg_length,
            'count': count,
            'is_rev': is_rev
        })
    
    graph = build_graph_from_occurrences(occ_list)
    q_len = len(query_seq)
    heap = []
    cost_arr = [float('inf')] * (q_len + 1)
    cost_arr[0] = 0
    visited = set()
    prev_state = {}
    heapq.heappush(heap, (0, (0, 0, 0, 0)))
    
    while heap:
        curr_total_cost, state = heapq.heappop(heap)
        curr_r, curr_q, local_cost, rev_flag = state
        if (curr_r, curr_q) in visited:
            continue
        visited.add((curr_r, curr_q))
        
        next_steps = get_next_shifts(query_seq, ref_seq, graph, curr_r, curr_q)
        if next_steps is None:
            continue
        for next_r, next_q, step_cost, next_rev in next_steps:
            if next_q >= q_len:
                route = []
                curr_state = (curr_r, curr_q)
                while curr_state != (0, 0):
                    route.append(curr_state)
                    curr_state = prev_state.get(curr_state, (0,0))
                route.reverse()
                show_alignment_details(route, query_seq, graph)
                return route
            if cost_arr[next_q] > cost_arr[curr_q] + step_cost:
                cost_arr[next_q] = cost_arr[curr_q] + step_cost
                heapq.heappush(heap, (cost_arr[next_q], (next_r, next_q, step_cost, next_rev)))
                prev_state[(next_r, next_q)] = (curr_r, curr_q)
    
    print("Failed to complete the alignment")