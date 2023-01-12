from cpython.memoryview cimport PyMemoryView_FromMemory
from cpython.buffer cimport PyBUF_READ
from cpython.mem cimport PyMem_Realloc, PyMem_Free
from cpython.unicode cimport PyUnicode_GET_LENGTH
from cpython.object cimport PyObject_GetAttr, PyTypeObject, Py_TYPE
from cpython.type cimport PyType_CheckExact
from libc.stdint cimport uint8_t
from libc.string cimport memset

cdef extern from "Python.h":
    void *PyUnicode_DATA(object)

from dnaio import SequenceRecord

if not PyType_CheckExact(SequenceRecord):
    raise RuntimeError("SequenceRecord is not a type class")

cdef PyTypeObject *sequence_record_class = <PyTypeObject *>SequenceRecord

# Nice trick from fastp: A,C, G, T, N all have different last three
# bits. So this requires only 8 entries per count array. Fastp performs
# a bitwise and of 0b111 on every character.
# This can be taken further by using  a lookup table. A=1, C=2, G=3, T=4.
# Lowercase a,c,g and t are supported. All other characters are index 0 and
# are considered N. This way we can make a very dense table, and don't have to
# check every nucleotide if it is within bounds. Furthermore, odd characters
# such as IUPAC K will map to N, unlike the fastp method where K will map to C.
cdef extern from *:
    """
    #include <stdint.h>
    const uint8_t NUCLEOTIDE_TO_INDEX[128] = {
        // Control characters
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // Interpunction numbers etc
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // A, B, C, D, E, F, G, H, I, J, K, L, M, N, O,
        0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
     // P, Q, R, S, T, U, V, W, X, Y, Z,  
        0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        // a, b, c, d, e, f, g, h, i, j, k, l, m, n, o,
        0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
     // p, q, r, s, t, u, v, w, x, y, z, 
        0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    };
    """
    const uint8_t NUCLEOTIDE_TO_INDEX[128]

ctypedef size_t counter_t
DEF PHRED_MAX=93
# Any phred value above the limit is capped at the limit. This way the
# table is significantly reduced in size. With a phred limit of 50, there
# are 51 possible threads. With 5 nucleotides and a 8 byte counter_t this
# means the table is slightly smaller than 2k. Which means we can fit 2 in a
# memory page. This improves performance and memory usage. Phreds above 50
# are rare (non existent?) and for the average higher phred scores do not
# contribute much compared to lower phred scores. So capping at 50 will have
# no meaningful effect on the resulting stats.
DEF PHRED_LIMIT=50
DEF PHRED_TABLE_SIZE=PHRED_LIMIT + 1
# Nice trick from fastp: A,C, G, T, N all have different last three
# bits. So this requires only 8 entries per count array.
DEF NUC_TABLE_SIZE=5
NUMBER_OF_NUCS = NUC_TABLE_SIZE
NUMBER_OF_PHREDS = PHRED_TABLE_SIZE
TABLE_SIZE = PHRED_TABLE_SIZE * NUC_TABLE_SIZE

ctypedef counter_t[PHRED_TABLE_SIZE][NUC_TABLE_SIZE] counttable_t
cdef class QCMetrics:
    cdef:
        object seq_name
        object qual_name
        uint8_t phred_offset
        # Illumina reads often use a limited range of qualities across the
        # sprectrum. Putting phreds first therefore gives nice dense spots
        # in the table that are accessed often leading to better cache locality.
        counttable_t *count_tables
        readonly Py_ssize_t max_length
        readonly size_t number_of_reads
    def __cinit__(self):
        self.seq_name = "sequence"
        self.qual_name = "qualities"
        self.phred_offset = 33
        self.max_length = 0
        self.count_tables = NULL

    def __dealloc__(self):
        PyMem_Free(self.count_tables)

    cdef _resize(self, Py_ssize_t new_size):
        cdef:
            counttable_t *tmp
        tmp =<counttable_t *>PyMem_Realloc(self.count_tables, new_size * sizeof(counttable_t))
        if tmp == NULL:
            raise MemoryError()
        self.count_tables = tmp
        memset(self.count_tables + self.max_length, 0, (new_size - self.max_length) * sizeof(counttable_t))
        self.max_length = new_size
        self.number_of_reads = 0
        
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
            self._resize(sequence_length)

        self.number_of_reads += 1
        for i in range(<size_t>sequence_length):
            c=sequence[i]
            q=qualities[i] - self.phred_offset
            if q > PHRED_MAX:
                raise ValueError(f"Not a valid phred character: {<char>qualities[i]}")
            if q > PHRED_LIMIT:
                q = PHRED_LIMIT
            c_index = NUCLEOTIDE_TO_INDEX[c]
            self.count_tables[i][q][c_index] += 1

    def count_table_view(self):
            return PyMemoryView_FromMemory(
                <char *>self.count_tables,
                self.max_length * sizeof(counttable_t),
                PyBUF_READ)
