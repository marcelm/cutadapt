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

cdef class QCMetrics:
    cdef:
        object seq_name
        object qual_name

    def __cinit__(self):
        self.seq_name = "sequence"
        self.qual_name = "qualities"


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
            Py_ssize_t sequence_length = PyUnicode_GET_LENGTH(sequence_obj)
            Py_ssize_t qualities_length = PyUnicode_GET_LENGTH(qualities_obj)
