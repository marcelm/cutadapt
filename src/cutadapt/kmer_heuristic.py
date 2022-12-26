import itertools
import sys
from typing import List, Optional, Set, Tuple
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


# A SearchSet is a start and stop combined with a list of possible kmer sets
# which should appear between this start and stop. Start and stop follow python
# indexing rules. (Negative start is a position relative to the end. None end
# is to the end of the sequence)
SearchSet = Tuple[int, Optional[int], List[Set[str]]]


def find_optimal_kmers(search_sets: List[SearchSet]) -> List[Tuple[str, int]]:
    minimal_score = sys.maxsize
    best_combination = None
    positions = [(start, stop) for start, stop, kmer_set_list in search_sets]
    kmer_set_lists = [kmer_set_list for start, stop, kmer_set_list in search_sets]
    for kmer_sets in itertools.product(*kmer_set_lists):
        check_set = set()
        for s in kmer_sets:
            check_set |= s
        if len(check_set) < minimal_score:
            best_combination = kmer_sets
            minimal_score = len(check_set)
    kmer_and_offsets_dict = defaultdict(list)
    for position, kmer_set in zip(positions, best_combination):  # type: ignore
        for kmer in kmer_set:
            kmer_and_offsets_dict[kmer].append(position)
    kmers_and_positions: List[Tuple[str, int, int]] = []
    for kmer, postions in kmer_and_offsets_dict.items():
        if len(positions) == 1:
            kmers_and_positions.append((kmer, *positions[0]))
            continue

        starts = [start for start, stop in postions]
        stops = [stop for start, stop in postions]
        if 0 in starts:
            start = 0
        # If all offsets have the same sign then choose the offset closest to
        # the start.
        elif all(start < 0 for start in starts) or all(start > 0 for start in starts):
            start = min(starts)
        else:  # Mixed positive and negative: search the entire sequence
            start = 0
        if None in stops:
            stop = None
        elif all(stop < 0 for stop in stops) or all(stop > 0 for stop in stops):
            stop = max(stops)
        else:
            stop = None
        kmers_and_positions.append((kmer, start, stop))
    return kmers_and_positions


def create_back_overlap_searchsets(
    adapter: str, min_overlap: int, error_rate: float
) -> List[SearchSet]:
    adapter_length = len(adapter)
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
                    None,
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
    min_overlap_kmer_start = (
        -(error_lengths[0] - 1) if error_lengths else -(adapter_length - 1)
    )
    search_sets.append(
        (
            min_overlap_kmer_start,
            None,
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
        start = -(next_length - 1)
        number_of_errors = i + 1
        kmer_sets = kmer_possibilities(adapter[:error_length], number_of_errors + 1)
        search_sets.append((start, None, kmer_sets))
    return search_sets


def create_kmers_and_offsets(
    adapter: str, min_overlap: int, error_rate: float
) -> List[Tuple[str, int]]:
    max_errors = int(len(adapter) * error_rate)
    search_sets = create_back_overlap_searchsets(adapter, min_overlap, error_rate)
    kmer_sets = kmer_possibilities(adapter, max_errors + 1)
    search_sets.append((0, None, kmer_sets))
    return find_optimal_kmers(search_sets)


def kmer_probability_analysis(
    kmers_and_offsets: List[Tuple[str, int, Optional[int]]], default_length: int = 150
):
    print("kmer\tstart\tstop\tconsidered sites\thit chance by random sequence (%)")
    accumulated_not_hit_chance = 1.0
    for kmer, start, stop in kmers_and_offsets:
        kmer_length = len(kmer)
        if stop is None:
            check_length = -start if start < 0 else default_length - start
        else:
            start = default_length - start if start < 0 else start
            check_length = max(stop - start, 0)
        considered_sites = check_length - kmer_length + 1
        single_kmer_hit_chance = 1 / 4**kmer_length
        not_hit_chance = (1 - single_kmer_hit_chance) ** considered_sites
        accumulated_not_hit_chance *= not_hit_chance
        print(
            f"{kmer:10}\t{start}\t{stop}\t{considered_sites}\t{(1 - not_hit_chance) * 100:.2f}"
        )
    print(
        f"Chance for profile hit by random sequence: {(1 - accumulated_not_hit_chance) * 100:.2f}%"
    )
