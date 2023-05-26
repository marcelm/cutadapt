from typing import List, Optional, Tuple

MAXIMUM_WORD_SIZE: int

class KmerFinder:
    def __init__(
        self,
        positions_and_kmers: List[Tuple[int, Optional[int], List[str]]],
        ref_wildcards: bool = False,
        query_wildcards: bool = False,
    ): ...
    def kmers_present(self, __sequence: str) -> bool: ...
