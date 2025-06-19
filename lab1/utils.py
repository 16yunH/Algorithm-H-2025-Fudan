def reverse_complement(dna_seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(nuc, nuc) for nuc in reversed(dna_seq))

def read_seq(filename):
    with open(filename, 'r') as f:
        return f.read().strip()