import itertools
import sys
from typing import List, Set, Tuple
from collections import defaultdict

def kmer_possibilities(sequence: str, chunks: int) -> List[Set[str]]:
    """
    Partition a sequence in almost equal sized chunks. Return all possibilities.

    Example sequence ABCDEFGH with 3 chunks. Possibilities:
    ["ABC", "DEF", "GH"]; ["ABC", "DE", "FGH"]; ["AB", "CDE", "FGH"]

    :param sequence: The sequence to b
    :param chunks:
    :return: A list of lists with all kmers
    """
    chunk_size = len(sequence) // (chunks)
    remainder = len(sequence) % (chunks)
    chunk_sizes: List[int] = remainder * [chunk_size + 1] + (chunks - remainder) * [
        chunk_size
    ]
    possible_orderings = set(itertools.permutations(chunk_sizes))
    kmer_sets = []
    for chunk_list in possible_orderings:
        offset = 0
        chunk_set = set()
        for size in chunk_list:
            chunk_set.add(sequence[offset : offset + size])
            offset += size
        kmer_sets.append(chunk_set)
    return kmer_sets


# A SearchSet is an offset combined with a list of possible kmer sets which
# should appear after this offset
SearchSet = Tuple[int, List[Set[str]]]


def find_optimal_kmers(search_sets: List[SearchSet]) -> List[Tuple[str, int]]:
    minimal_score = sys.maxsize
    best_combination = None
    offsets = [offset for offset, kmer_set_list in search_sets]
    kmer_set_lists = [kmer_set_list for offset, kmer_set_list in search_sets]
    for kmer_sets in itertools.product(*kmer_set_lists):
        check_set = set()
        for s in kmer_sets:
            check_set |= s
        if len(check_set) < minimal_score:
            best_combination = kmer_sets
            minimal_score = len(check_set)
    kmer_and_offsets_dict = defaultdict(list)
    for offset, kmer_set in zip(offsets, best_combination):
        for kmer in kmer_set:
            kmer_and_offsets_dict[kmer].append(offset)
    kmers_and_offsets: List[Tuple[str, int]] = []
    for kmer, offsets in kmer_and_offsets_dict.items():
        if len(offsets) == 1:
            offset = offsets[0]
        elif 0 in offsets:
            offset = 0
        # If all offsets have the same sign then choose the offset closest to
        # the start.
        elif all(offset < 0 for offset in offsets) or all(
            offset > 0 for offset in offsets
        ):
            offset = min(offsets)
        else:  # Mixed positive and negative: search the entire sequence
            offset = 0
        kmers_and_offsets.append((kmer, offset))
    return kmers_and_offsets


def create_kmers_and_offsets(
    adapter: str, min_overlap: int, error_rate: float
) -> List[Tuple[str, int]]:
    adapter_length = len(adapter)
    max_errors = int(adapter_length * error_rate)
    error_lengths = []
    max_error = 1
    search_sets: List[SearchSet] = []
    for i in range(adapter_length + 1):
        if i * error_rate >= max_error:
            error_lengths.append(i)
            max_error += 1

    # Add a couple of directly matching 1, 2, 3 and 4-mer searches.
    # The probability of a false positive is just to high when for example
    # a 3-mer is evaluated in more than one position.
    min_overlap_kmer_length = 5
    if min_overlap < min_overlap_kmer_length:
        for i in range(min_overlap, min_overlap_kmer_length):
            search_sets.append(
                (
                    -i,
                    [
                        {
                            adapter[:i],
                        }
                    ],
                )
            )
        min_overlap = min_overlap_kmer_length
    # Build up the array with chunks which should occur at the tail end
    # if the adapter overlaps with the end.
    min_overlap_kmer = adapter[:min_overlap]
    min_overlap_kmer_offset = -(error_lengths[0] - 1)
    search_sets.append(
        (
            min_overlap_kmer_offset,
            [
                {
                    min_overlap_kmer,
                }
            ],
        )
    )
    for i, error_length in enumerate(error_lengths):
        if (i + 1) < len(error_lengths):
            next_length = error_lengths[i + 1]
        else:
            next_length = adapter_length
        offset = -(next_length - 1)
        number_of_errors = i + 1
        kmer_sets = kmer_possibilities(adapter[:error_length], number_of_errors + 1)
        search_sets.append((offset, kmer_sets))

    # Create kmers at least one of which should be in the read when there is a
    # an adapter
    kmer_sets = kmer_possibilities(adapter, max_errors + 1)
    search_sets.append((0, kmer_sets))
    return find_optimal_kmers(search_sets)


def kmer_probability_analysis(
    kmers_and_offsets: List[Tuple[str, int]], default_length: int = 150
):
    print("kmer\toffset\tconsidered sites\thit chance by random sequence (%)")
    accumulated_not_hit_chance = 1.0
    for kmer, offset in kmers_and_offsets:
        kmer_length = len(kmer)
        check_length = -offset if offset < 0 else default_length - offset
        considered_sites = check_length - kmer_length + 1
        single_kmer_hit_chance = 1 / 4**kmer_length
        not_hit_chance = (1 - single_kmer_hit_chance) ** considered_sites
        accumulated_not_hit_chance *= not_hit_chance
        print(
            f"{kmer:10}\t{offset}\t{considered_sites}\t{(1 - not_hit_chance) * 100:.2f}"
        )
    print(
        f"Chance for profile hit by random sequence: {(1 - accumulated_not_hit_chance) * 100:.2f}%"
    )
