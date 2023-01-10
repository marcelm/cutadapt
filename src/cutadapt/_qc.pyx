from cpython.memoryview cimport PyMemoryView_FromMemory
from cpython.buffer cimport PyBUF_READ
from cpython.mem cimport PyMem_Realloc, PyMem_Free
from cpython.unicode cimport PyUnicode_GET_LENGTH
from cpython.object cimport PyObject_GetAttr, PyTypeObject, Py_TYPE
from cpython.type cimport PyType_CheckExact
from libc.stdint cimport uint8_t

cdef extern from "Python.h":
    void *PyUnicode_DATA(object)

from dnaio import SequenceRecord

if not PyType_CheckExact(SequenceRecord):
    raise RuntimeError("SequenceRecord is not a type class")

cdef PyTypeObject *sequence_record_class = <PyTypeObject *>SequenceRecord

DEF PHRED_MAX=93
DEF NUC_MAX=5
ctypedef size_t counter_t

cdef inline uint8_t nucleotide_index_from_char(uint8_t c):
    # Nice trick from fastp: A,C, G, T, N all have different last three
    # bits. So this requires only 8 entries per count array.
    # We can take this further. Note that by design, upper and lowercase
    # ASCII alphabetic characters share the last 5 bits.
    # char     &0b1111  >> 1   &0b11
    # A, a     0001     000    00
    # C, c     0011     001    01
    # G, g     0111     011    11
    # T, t     0100     010    10
    # N, n     1110     111    11
    # N >> 3 & 0b1 == 1. For A,C,G,T this is 0.
    # The following calculation is branchless and guaranteed to be
    # 0 (A), 1 (C), 2 (T), 3 (G) or 4 (N). Works also with lowercase.
    # Other characters don't throw errors. It is garbage in, garbage out.
    # The guarantees mean no bounds check is required on the result.
    return ((c & 0b111) >> 1) + ((c >> 3) & 0b1)

cdef class QCMetrics:
    cdef:
        object seq_name
        object qual_name
        uint8_t phred_offset
        counter_t[PHRED_MAX][NUC_MAX] *count_tables
        Py_ssize_t max_length

    def __cinit__(self):
        self.seq_name = "sequence"
        self.qual_name = "qualities"
        self.phred_offset = 33
        self.max_length = 0
        self.count_tables = NULL

    def __dealloc__(self):
        PyMem_Free(self.count_tables)

    def add_read(self, read):
        if Py_TYPE(read) != sequence_record_class:
            raise TypeError(f"read should be a dnaio.SequenceRecord object, "
                            f"got {type(read)}")
        cdef:
            # Use PyObject_GetAttr rather than relying on Cython or using
            # PyObject_GetAttrString, that will create a new unicode object
            # With PyUnicode_FromString *every single time*.
            object sequence_obj = PyObject_GetAttr(read, self.seq_name)
            object qualities_obj = PyObject_GetAttr(read, self.qual_name)
            # Thanks to guarantees in dnaio, we know all strings are ascii
            const uint8_t *sequence = <uint8_t *>PyUnicode_DATA(sequence_obj)
            const uint8_t *qualities = <uint8_t *>PyUnicode_DATA(qualities_obj)
            # dnaio guarantees sequence and qualities have the same length.
            Py_ssize_t sequence_length = PyUnicode_GET_LENGTH(sequence_obj)
            size_t i
            uint8_t c, q
            uint8_t c_index
        if sequence_length > self.max_length:
            self.count_tables =<counter_t (*)[PHRED_MAX][NUC_MAX]>PyMem_Realloc(self.count_tables, sequence_length * sizeof(counter_t) * NUC_MAX * PHRED_MAX)
            self.max_length = sequence_length

        for i in range(<size_t>sequence_length):
            c=sequence[i]
            q=qualities[i] - self.phred_offset
            if q > PHRED_MAX:
                raise ValueError(f"Not a valid phred character: {<char>qualities[i]}")
            c_index = nucleotide_index_from_char(c)
            self.count_tables[i][c_index][q] += 1

def count_table_view(self):
        return PyMemoryView_FromMemory(
            <char *>self.count_tables,
            self.max_length * sizeof(counter_t) * NUC_MAX * PHRED_MAX,
            PyBUF_READ)
