from dnaio import SequenceRecord
from cutadapt.adapters import Match

class ModificationInfo:
    matches: list[Match] = ...
    original_read: SequenceRecord
    cut_prefix: str
    cut_suffix: str
    is_rc: bool | None
    def __init__(self, read: SequenceRecord) -> None: ...
    def __repr__(self) -> str: ...
