from utils import read_seq
from alignment import optimal_alignment

def main():
    query = read_seq("data/query2.txt")
    reference = read_seq("data/reference2.txt")
    optimal_alignment(query, reference)

if __name__ == "__main__":
    main()