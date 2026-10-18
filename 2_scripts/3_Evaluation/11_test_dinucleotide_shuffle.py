from collections import Counter
from pathlib import Path
import random

DNA = set("ACGT")

def dinucleotide_counts(sequence):
    return Counter(
        sequence[i:i+2]
        for i in range(len(sequence) - 1)
    )

def shuffle_dinucleotides(sequence, rng):
    """Randomize a sequence while preserving its dinucleotide counts."""
    sequence = sequence.upper()

    if len(sequence) < 2 or not set(sequence) <= DNA:
        return None

    # Build a directed multigraph from adjacent bases.
    edges = {}
    for a, b in zip(sequence[:-1], sequence[1:]):
        edges.setdefault(a, []).append(b)

    # Randomize edge order at every base.
    for base in edges:
        rng.shuffle(edges[base])

    # Reconstruct an Eulerian trail using Hierholzer's algorithm.
    stack = [sequence[0]]
    trail = []

    while stack:
        base = stack[-1]
        if edges.get(base):
            stack.append(edges[base].pop())
        else:
            trail.append(stack.pop())

    shuffled = "".join(reversed(trail))
    return shuffled if len(shuffled) == len(sequence) else None

rng = random.Random(42)

examples = [
    "ACGTACGTACGT",
    "AAAACCCCGGGGTTTT",
    "ATATCGCGATATCGCG",
    "ACGTTGCAACGTAGCT",
]

for original in examples:
    shuffled = shuffle_dinucleotides(original, rng)

    if shuffled is None:
        raise RuntimeError("Shuffling failed.")

    original_counts = dinucleotide_counts(original)
    shuffled_counts = dinucleotide_counts(shuffled)

    print("Original: ", original)
    print("Shuffled: ", shuffled)
    print("Same dinucleotide counts:", original_counts == shuffled_counts)

    assert original_counts == shuffled_counts
    assert Counter(original) == Counter(shuffled)
    print()

print("PASS: all test sequences preserve mono- and dinucleotide counts.")
