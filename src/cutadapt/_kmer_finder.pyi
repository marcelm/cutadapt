MAXIMUM_WORD_SIZE: int

class KmerFinder:
    def __init__(
        self,
        positions_and_kmers: list[tuple[int, int | None, list[str]]],
        ref_wildcards: bool = False,
        query_wildcards: bool = False,
    ): ...
    def kmers_present(self, /, sequence: str) -> bool: ...
