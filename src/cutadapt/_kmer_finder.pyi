from typing import List, Optional, Tuple

class KmerFinder:
    def __init__(
        self,
        kmers_and_positions: List[Tuple[int, Optional[int], List[str]]],
        ref_wildcards: bool = False,
        query_wildcards: bool = False,
    ): ...
    def kmers_present(self, __sequence: str) -> bool: ...
