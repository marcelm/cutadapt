from cpython.object cimport PyObject_GetAttr
from cpython.unicode cimport PyUnicode_GET_LENGTH
from cpython.object cimport PyTypeObject, Py_TYPE
from cpython.type cimport PyType_CheckExact
from libc.stdint cimport uint8_t

cdef extern from "Python.h":
    void *PyUnicode_DATA(object)

from dnaio import SequenceRecord

if not PyType_CheckExact(SequenceRecord):
    raise RuntimeError("SequenceRecord is not a type class")

cdef PyTypeObject *sequence_record_class = <PyTypeObject *>SequenceRecord

DEF PHRED_MAX=93

cdef class QCMetrics:
    cdef:
        object seq_name
        object qual_name
        uint8_t phred_offset

    def __cinit__(self):
        self.seq_name = "sequence"
        self.qual_name = "qualities"
        self.phred_offset = 33


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
        for i in range(<size_t>sequence_length):
            c=sequence[i]
            q=qualities[i] - self.phred_offset
            if q > PHRED_MAX:
                raise ValueError(f"Not a valid phred character: {<char>qualities[i]}")
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
            c_index = ((c & 0b111) >> 1) + ((c >> 3) & 0b1)

