from typing import List, Tuple

class KmerFinder:
    def __init__(
        self,
        kmers_and_offsets: List[Tuple[str, int]],
        ref_wildcards: bool,
        query_wildcards: bool,
    ): ...
    def kmers_present(self, __sequence: str) -> bool: ...
