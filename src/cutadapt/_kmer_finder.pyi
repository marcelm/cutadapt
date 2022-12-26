from typing import List, Optional, Tuple

class KmerFinder:
    def __init__(
        self,
        kmers_and_positions: List[Tuple[str, int, Optional[int]]],
        ref_wildcards: bool,
        query_wildcards: bool,
    ): ...
    def kmers_present(self, __sequence: str) -> bool: ...
