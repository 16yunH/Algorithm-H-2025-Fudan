class Graph:
    def __init__(self):
        self.adjacency = {}
        
    def add_edge(self, src, dest, weight, meta):
        self.adjacency.setdefault(src, {})[dest] = {'weight': weight, 'meta': meta}
        
    def edge_exists(self, src, dest):
        return src in self.adjacency and dest in self.adjacency[src]
    
    def neighbors(self, src):
        if src not in self.adjacency:
            return []
        return list(self.adjacency[src].keys())
    
    def get_meta(self, src, dest):
        return self.adjacency[src][dest]['meta']
    
    def get_weight(self, src, dest):
        return self.adjacency[src][dest]['weight']