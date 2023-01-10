from dnaio import SequenceRecord

TABLE_SIZE: int
NUMBER_OF_PHREDS: int
NUMBER_OF_NUCS: int

class QCMetrics:
    number_of_reads: int
    def __init__(self): ...
    def add_read(self, __read: SequenceRecord): ...
    def count_table_view(self) -> memoryview: ...
