"""
SVG可视化模块
生成匹配结果的交互式可视化图像
"""

def generate_svg_grid(query_s_len, ref_s_len, 
                      fwd_segs, fwd_pts,
                      rc_segs, rc_pts,
                      chunk_sz,
                      num_q_chunks, num_r_chunks,
                      output_filename="matches_grid_visualization.svg",
                      svg_padding=70, cell_size=10, point_radius=3,
                      long_path=None):
    """
    生成匹配结果的二维网格 SVG 可视化图像。
    包含线段、孤立点和最长路径高亮。
    """
    grid_width = num_q_chunks * cell_size
    grid_height = num_r_chunks * cell_size    # 处理序列或区块数量为零的情况
    if num_q_chunks == 0 or num_r_chunks == 0:
        grid_width = max(grid_width, cell_size) # 至少保证一个单元格的宽度
        grid_height = max(grid_height, cell_size) # 至少保证一个单元格的高度

    svg_width = grid_width + 2 * svg_padding
    svg_height = grid_height + 2 * svg_padding

    svg_content = f'<svg width="{svg_width}" height="{svg_height}" xmlns="http://www.w3.org/2000/svg" font-family="Arial, sans-serif">\n'
    svg_content += '  <style>\n'
    svg_content += '    .grid_line { stroke: #E0E0E0; stroke-width: 0.5; }\n'
    svg_content += '    .fwd_match_node { fill: blue; fill-opacity: 0.7; }\n' # 正向匹配点样式
    svg_content += '    .rc_match_node { fill: red; fill-opacity: 0.7; }\n'   # 反向互补匹配点样式
    svg_content += '    .fwd_match_segment { stroke: blue; stroke-width: 2; stroke-opacity: 0.7; }\n' # 正向匹配线段样式
    svg_content += '    .rc_match_segment { stroke: red; stroke-width: 2; stroke-opacity: 0.7; }\n'   # 反向互补匹配线段样式
    svg_content += '    .longest_path_line { stroke: purple; stroke-width: 2.5; stroke-opacity: 0.9; }\n' # 最长路径连线样式
    svg_content += '    .longest_path_node_highlight { stroke: gold; stroke-width: 2px; fill-opacity: 1 !important; }\n' # 最长路径节点高亮样式
    svg_content += '    .axis_label { font-size: 10px; text-anchor: middle; }\n' # 坐标轴标签样式
    svg_content += '    .title_label { font-size: 14px; text-anchor: middle; font-weight: bold; }\n' # 标题样式
    svg_content += '    .no_match_text { font-size: 12px; text-anchor: middle; fill: #555; }\n' # 无匹配时提示文本样式
    svg_content += '    circle:hover { stroke: black; stroke-width: 1px; }\n' # 鼠标悬停效果
    svg_content += '  </style>\n'

    # 图表标题
    svg_content += f'  <text x="{svg_width / 2}" y="{svg_padding / 2.5}" class="title_label">查询块 vs 参考块 匹配网格</text>\n'    # 绘制网格背景和边框
    svg_content += f'  <rect x="{svg_padding}" y="{svg_padding}" width="{grid_width}" height="{grid_height}" fill="none" stroke="#333" stroke-width="1"/>\n'
    
    if num_q_chunks == 0 or num_r_chunks == 0:
        svg_content += f'  <text x="{svg_padding + grid_width / 2}" y="{svg_padding + grid_height / 2}" class="no_match_text">无数据显示 (查询或参考块数量为零)。</text>\n'
    else:
        # 绘制网格线
        for i in range(1, num_q_chunks): # 内部垂直线
            x = svg_padding + i * cell_size
            svg_content += f'  <line x1="{x}" y1="{svg_padding}" x2="{x}" y2="{svg_padding + grid_height}" class="grid_line" />\n'
        for j in range(1, num_r_chunks): # 内部水平线
            y = svg_padding + j * cell_size
            svg_content += f'  <line x1="{svg_padding}" y1="{y}" x2="{svg_padding + grid_width}" y2="{y}" class="grid_line" />\n'
        
        # 绘制轴标签
        # X 轴标签 (查询块)
        label_step_x = max(1, num_q_chunks // 20) # 调整标签密度避免重叠
        for i in range(num_q_chunks):
            if i % label_step_x == 0:
                x = svg_padding + (i + 0.5) * cell_size
                svg_content += f'  <text x="{x}" y="{svg_padding + grid_height + 15}" class="axis_label">{i}</text>\n'
        svg_content += f'  <text x="{svg_padding + grid_width / 2}" y="{svg_padding + grid_height + 35}" class="axis_label">查询块索引 (块大小: {chunk_sz}bp, 总块数: {num_q_chunks})</text>\n'
        
        # Y 轴标签 (参考块)
        label_step_y = max(1, num_r_chunks // 20)
        for j in range(num_r_chunks):
            if j % label_step_y == 0:
                y_pos = svg_padding + (j + 0.5) * cell_size
                x_pos = svg_padding - 15
                svg_content += f'  <text x="{x_pos}" y="{y_pos}" class="axis_label" transform="rotate(-90 {x_pos},{y_pos})">{j}</text>\n'
        svg_content += f'  <text x="{svg_padding - 45}" y="{svg_padding + grid_height / 2}" class="axis_label" transform="rotate(-90 {svg_padding - 45},{svg_padding + grid_height / 2})">参考块索引 (总块数: {num_r_chunks})</text>\n'        # 绘制匹配节点和线段
        # 正向匹配 (蓝色)
        for q_s, r_s, q_e, r_e in fwd_segs: # 绘制线段
            x1 = svg_padding + (q_s + 0.5) * cell_size
            y1 = svg_padding + (r_s + 0.5) * cell_size
            x2 = svg_padding + (q_e + 0.5) * cell_size
            y2 = svg_padding + (r_e + 0.5) * cell_size
            svg_content += f'  <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" class="fwd_match_segment">\n'
            svg_content += f'    <title>正向线段\n查询块: {q_s}-{q_e}\n参考块: {r_s}-{r_e}</title>\n'
            svg_content += f'  </line>\n'
        for q_idx, r_idx, score in fwd_pts: # 绘制孤立点
            center_x = svg_padding + (q_idx + 0.5) * cell_size
            center_y = svg_padding + (r_idx + 0.5) * cell_size
            svg_content += f'  <circle cx="{center_x}" cy="{center_y}" r="{point_radius}" class="fwd_match_node">\n'
            svg_content += f'    <title>正向匹配 (孤立点)\n查询块: {q_idx}\n参考块: {r_idx}\n分数: {score}</title>\n'
            svg_content += f'  </circle>\n'

        # 反向互补匹配 (红色)
        for q_s, r_s, q_e, r_e in rc_segs: # 绘制线段
            x1 = svg_padding + (q_s + 0.5) * cell_size
            y1 = svg_padding + (r_s + 0.5) * cell_size
            x2 = svg_padding + (q_e + 0.5) * cell_size
            y2 = svg_padding + (r_e + 0.5) * cell_size
            svg_content += f'  <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" class="rc_match_segment">\n'
            svg_content += f'    <title>反向互补线段\n查询块: {q_s}-{q_e}\n参考块: {r_s}-{r_e}</title>\n'
            svg_content += f'  </line>\n'
        for q_idx, r_idx, score in rc_pts: # 绘制孤立点
            center_x = svg_padding + (q_idx + 0.5) * cell_size
            center_y = svg_padding + (r_idx + 0.5) * cell_size
            svg_content += f'  <circle cx="{center_x}" cy="{center_y}" r="{point_radius}" class="rc_match_node">\n'
            svg_content += f'    <title>反向互补匹配 (孤立点)\n查询块: {q_idx}\n参考块: {r_idx}\n分数: {score}</title>\n'
            svg_content += f'  </circle>\n'

        # 高亮显示最长路径（如果提供）
        if long_path and len(long_path) > 0:
            if len(long_path) > 1: # 绘制连接线
                for p_idx in range(len(long_path) - 1):
                    curr_n = long_path[p_idx]
                    next_n = long_path[p_idx+1]
                    x1 = svg_padding + (curr_n['q'] + 0.5) * cell_size
                    y1 = svg_padding + (curr_n['r'] + 0.5) * cell_size
                    x2 = svg_padding + (next_n['q'] + 0.5) * cell_size
                    y2 = svg_padding + (next_n['r'] + 0.5) * cell_size
                    svg_content += f'  <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" class="longest_path_line" />\n'
            for n_data in long_path: # 高亮路径上的节点
                center_x = svg_padding + (n_data['q'] + 0.5) * cell_size
                center_y = svg_padding + (n_data['r'] + 0.5) * cell_size
                base_class = "fwd_match_node" if n_data['type'] == 'forward' else "rc_match_node"
                svg_content += f'  <circle cx="{center_x}" cy="{center_y}" r="{point_radius + 1}" class="{base_class} longest_path_node_highlight">\n' # 高亮节点半径稍大
                svg_content += f'    <title>{n_data["type"].capitalize()} 匹配 (最长路径)\n查询块: {n_data["q"]}\n参考块: {n_data["r"]}\n原始分数: {n_data["original_score"]}\n路径分数: {n_data.get("path_score_at_node", "N/A"):.2f}</title>\n'
                svg_content += f'  </circle>\n'
    svg_content += '</svg>\n'
    
    with open(output_filename, 'w', encoding='utf-8') as f:
        f.write(svg_content)
    print(f"SVG 网格可视化已保存到 {output_filename}")
