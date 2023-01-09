import io
import itertools
import sys
from typing import List, Optional, Set, Tuple
from collections import defaultdict


def kmer_possibilities(sequence: str, chunks: int) -> List[Set[str]]:
    """
    Partition a sequence in almost equal sized chunks. Return all possibilities.

    Example sequence ABCDEFGH with 3 chunks. Possibilities:
    ["ABC", "DEF", "GH"]; ["ABC", "DE", "FGH"]; ["AB", "CDE", "FGH"]
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


def minimize_kmer_search_list(
    kmer_search_list: List[Tuple[str, int, Optional[int]]]
) -> List[Tuple[str, int, Optional[int]]]:
    kmer_and_offsets_dict = defaultdict(list)
    for kmer, start, stop in kmer_search_list:  # type: ignore
        kmer_and_offsets_dict[kmer].append((start, stop))
    kmers_and_positions: List[Tuple[str, int, Optional[int]]] = []
    for kmer, positions in kmer_and_offsets_dict.items():
        if len(positions) == 1:
            start, stop = positions[0]
            kmers_and_positions.append((kmer, start, stop))
            continue
        if (0, None) in positions:
            kmers_and_positions.append((kmer, 0, None))
            continue
        front_searches = [(start, stop) for start, stop in positions if start == 0]
        back_searches = [(start, stop) for start, stop in positions if stop is None]
        middle_searches = [
            (start, stop)
            for start, stop in positions
            if start != 0 and stop is not None
        ]
        if middle_searches:
            raise NotImplementedError(
                "Situations with searches starting in the middle have not been considered."
            )
        if front_searches:
            # (0, None) condition is already caught, so stop is never None.
            kmers_and_positions.append(
                (kmer, 0, max(stop for start, stop in front_searches))  # type: ignore
            )
        if back_searches:
            kmers_and_positions.append(
                (kmer, min(start for start, stop in back_searches), None)
            )
    return kmers_and_positions


def find_optimal_kmers(
    search_sets: List[SearchSet],
) -> List[Tuple[int, Optional[int], List[str]]]:
    minimal_score = sys.maxsize
    best_combination = []
    positions = [(start, stop) for start, stop, kmer_set_list in search_sets]
    kmer_set_lists = [kmer_set_list for start, stop, kmer_set_list in search_sets]
    for kmer_sets in itertools.product(*kmer_set_lists):
        kmer_search_list = []
        for kmer_set, (start, stop) in zip(kmer_sets, positions):
            for kmer in kmer_set:
                kmer_search_list.append((kmer, start, stop))
        minimized_search_list = minimize_kmer_search_list(kmer_search_list)
        if len(minimized_search_list) < minimal_score:
            best_combination = minimized_search_list
            minimal_score = len(minimized_search_list)
    result_dict = defaultdict(list)
    for kmer, start, stop in best_combination:
        result_dict[(start, stop)].append(kmer)
    return [(start, stop, kmers) for (start, stop), kmers in result_dict.items()]


def create_back_overlap_searchsets(
    adapter: str, min_overlap: int, error_rate: float
) -> List[SearchSet]:
    adapter_length = len(adapter)
    error_lengths = []
    max_error = 0
    search_sets: List[SearchSet] = []
    for i in range(adapter_length + 1):
        if int(i * error_rate) > max_error:
            error_lengths.append((max_error, i - 1))
            max_error += 1
    error_lengths.append((max_error, adapter_length))

    minimum_length = min_overlap
    for max_errors, length in error_lengths:
        if minimum_length > length:
            continue
        if max_errors == 0:
            # Add a couple of directly matching 1, 2, 3 and 4-mer searches.
            # The probability of a false positive is just to high when for
            # example a 3-mer is evaluated in more than one position.
            min_overlap_kmer_length = 5
            if minimum_length < min_overlap_kmer_length:
                for i in range(minimum_length, min_overlap_kmer_length):
                    search_set = (-i, None, [{adapter[:i]}])
                    search_sets.append(search_set)
                minimum_length = min_overlap_kmer_length
        kmer_sets = kmer_possibilities(adapter[:minimum_length], max_errors + 1)
        search_sets.append((-length, None, kmer_sets))
        minimum_length = length + 1
    return search_sets


def create_positions_and_kmers(
    adapter: str,
    min_overlap: int,
    error_rate: float,
    back_adapter: bool,
    front_adapter: bool,
    internal: bool = True,
) -> List[Tuple[int, Optional[int], List[str]]]:
    """
    Create a set of position and words combinations where at least one of the
    words needs to occur at its specified position. If not an alignment
    algorithm will not be able to find a solution. This can be checked very
    quickly and allows for skipping alignment in cases where the adapter would
    not align anyway.

    Example: looking for AAAAATTTTT with at most one error. This means either
    AAAAA or TTTTT (or both) must be present, otherwise alignment will not
    succeed.

    This function returns the positions and the accompanying words while also
    taking into account partial overlap for back and front adapters.
    """
    max_errors = int(len(adapter) * error_rate)
    search_sets = []
    if back_adapter:
        search_sets.extend(
            create_back_overlap_searchsets(adapter, min_overlap, error_rate)
        )
    if front_adapter:
        # To create a front adapter the code is practically the same except
        # with some parameters set differently. Reversing the adapter, running
        # the back adapter code and reversing all the kmers and positions has
        # the same effect without needing to duplicate the code.
        reversed_back_search_sets = create_back_overlap_searchsets(
            adapter[::-1], min_overlap, error_rate
        )
        front_search_sets = []
        for start, stop, kmer_sets in reversed_back_search_sets:
            new_kmer_sets = [
                {kmer[::-1] for kmer in kmer_set} for kmer_set in kmer_sets
            ]
            front_search_sets.append((0, -start, new_kmer_sets))
        search_sets.extend(front_search_sets)
    if internal:
        kmer_sets = kmer_possibilities(adapter, max_errors + 1)
        search_sets.append((0, None, kmer_sets))
    return find_optimal_kmers(search_sets)


def kmer_probability_analysis(
    kmers_and_offsets: List[Tuple[int, Optional[int], List[str]]],
    default_length: int = 150,
) -> str:  # pragma: no cover  # only for debugging use
    """
    Returns a tab separated table with for each kmer a start, stop, the number
    of considered sites and the hit chance on a randomly generated sequence
    containing only A, C, G and T. Assumes kmers only consist of A, C, G and T
    too.

    Useful for investigating whether the create_positions_and_kmers function
    creates a useful runtime heuristic.
    """
    out = io.StringIO()
    out.write(
        "kmer\tstart\tstop\tconsidered sites\thit chance by random sequence (%)\n"
    )
    accumulated_not_hit_chance = 1.0
    for start, stop, kmers in kmers_and_offsets:
        if stop is None:
            check_length = -start if start < 0 else default_length - start
        else:
            start = default_length - start if start < 0 else start
            check_length = max(stop - start, 0)
        for kmer in kmers:
            kmer_length = len(kmer)
            considered_sites = check_length - kmer_length + 1
            single_kmer_hit_chance = 1 / 4**kmer_length
            not_hit_chance = (1 - single_kmer_hit_chance) ** considered_sites
            accumulated_not_hit_chance *= not_hit_chance
            out.write(
                f"{kmer:10}\t{start}\t{stop}\t{considered_sites}\t{(1 - not_hit_chance) * 100:.2f}\n"
            )
    out.write(
        f"Chance for profile hit by random sequence: {(1 - accumulated_not_hit_chance) * 100:.2f}%\n"
    )
    return out.getvalue()
