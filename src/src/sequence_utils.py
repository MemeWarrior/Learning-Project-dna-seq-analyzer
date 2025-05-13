from typing import List, Tuple, Union
from Bio.Seq import Seq

# === Functions ===

SeqType = Union[str, Seq] #A custom type which for readability encapsulates the str and Seq types

def base_counts(seq: SeqType) -> dict[str, int]:
    return {nt: seq.count(nt) for nt in "ACGT"}

def gc_contents(seq: SeqType) -> float:

    gc = 0
    cc = 0

    for base in seq:
        if base == "G":
            gc+=1
        elif base == "C":
            cc+=1

    return 100 * ((gc+cc)/len(seq))

def sliding_gc(seq: SeqType, window_size: int, step: int) -> Tuple[List[int], List[float]]:

    positions = []
    gc_values = []

    for i in range(0, len(seq)-window_size+1,step):
        positions.append(i)
        gc_values.append(gc_contents(seq[i:i+window_size]))
        
    return positions, gc_values


