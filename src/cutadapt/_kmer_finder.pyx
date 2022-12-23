# cython: language_level=3

from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.string cimport memset, strlen
from cpython.unicode cimport PyUnicode_CheckExact, PyUnicode_GET_LENGTH
from libc.stdint cimport uint8_t

# Dnaio conveniently ensures that all sequences are ASCII only.
DEF ASCII_CHAR_COUNT = 128

cdef extern from "Python.h":
    void *PyUnicode_DATA(object o)
    bint PyUnicode_IS_COMPACT_ASCII(object o)

ctypedef struct KmerEntry:
    size_t kmer_length
    size_t mask_offset
    ssize_t search_offset


cdef class KmerFinder:
    """
    Find kmers in strings. To replace the following code:

        kmers_and_offsets = [("AGA", -10), ("AGCATGA", 0)]
        for kmer, offset in kmers_and_offsets:
            sequence.find(kmer, offset)

    This has a lot of python overhead. The following code is equivalent:

        kmers_and_offsets = [("AGA", -10), ("AGCATGA", 0)]
        kmer_finder = KmerFinder(kmers_and_offsets)
        kmer_finder.kmers_present(sequence)

    This is more efficient as the kmers_present method can be applied to a lot
    of sequences and all the necessary unpacking for each kmer into C variables
    happens only once.
    """
    cdef:
        KmerEntry *kmer_entries
        size_t *kmer_masks
        size_t number_of_kmers
        object kmers_and_offsets
        readonly bint ref_wildcards
        readonly bint query_wildcards

    def __cinit__(self, kmers_and_offsets, ref_wildcards=False, query_wildcards=False):
        self.kmer_masks = NULL
        self.kmer_entries = NULL
        self.number_of_kmers = 0
        self.ref_wildcards = ref_wildcards
        self.query_wildcards = query_wildcards
        kmers = [kmer for kmer, _ in kmers_and_offsets]
        kmer_total_length = sum(len(kmer) for kmer in kmers)
        number_of_entries = len(kmers_and_offsets)
        self.kmer_entries = <KmerEntry *>PyMem_Malloc(number_of_entries * sizeof(KmerEntry))
        # for the kmers the NULL bytes also need space.
        self.kmer_masks = <size_t *>PyMem_Malloc(number_of_entries * sizeof(size_t) * ASCII_CHAR_COUNT)
        self.number_of_kmers = number_of_entries
        cdef size_t mask_offset = 0
        cdef char *kmer_ptr
        cdef Py_ssize_t kmer_length
        for i, (kmer, offset) in enumerate(kmers_and_offsets):
            if not PyUnicode_CheckExact(kmer):
                raise TypeError(f"Kmer should be a string not {type(kmer)}")
            if not PyUnicode_IS_COMPACT_ASCII(kmer):
                raise ValueError("Only ASCII strings are supported")
            self.kmer_entries[i].search_offset  = offset
            self.kmer_entries[i].mask_offset = mask_offset
            kmer_length = PyUnicode_GET_LENGTH(kmer)
            self.kmer_entries[i].kmer_length = kmer_length
            kmer_ptr = <char *>PyUnicode_DATA(kmer)
            populate_needle_mask(self.kmer_masks + mask_offset, kmer_ptr, kmer_length,
                                 self.ref_wildcards, self.query_wildcards)
            mask_offset += ASCII_CHAR_COUNT
        self.kmers_and_offsets = kmers_and_offsets

    def __reduce__(self):
        return KmerFinder, (self.kmers_and_offsets,)

    def kmers_present(self, str sequence):
        cdef:
            KmerEntry entry
            size_t i
            size_t kmer_offset
            size_t kmer_length
            ssize_t search_offset
            size_t *mask_ptr
            char *search_ptr
            char *search_result
            size_t search_length
        if not PyUnicode_IS_COMPACT_ASCII(sequence):
            raise ValueError("Only ASCII strings are supported")
        cdef char *seq = <char *>PyUnicode_DATA(sequence)
        cdef Py_ssize_t seq_length = PyUnicode_GET_LENGTH(sequence)
        for i in range(self.number_of_kmers):
            entry = self.kmer_entries[i]
            search_offset = entry.search_offset
            if search_offset < 0:
                search_offset = seq_length + search_offset
                if search_offset < 0:
                    search_offset = 0
            if search_offset > seq_length:
                continue
            kmer_length = entry.kmer_length
            mask_ptr = self.kmer_masks + entry.mask_offset
            search_ptr = seq + search_offset
            search_length = seq_length - (search_ptr - seq)
            search_result = shift_and_search(search_ptr, search_length,
                                             mask_ptr, kmer_length)
            if search_result:
                return True
        return False

    def __dealloc__(self):
        PyMem_Free(self.kmer_masks)
        PyMem_Free(self.kmer_entries)


cdef void set_masks(size_t *needle_mask, size_t pos, char *chars):
    cdef char c
    cdef size_t i
    for i in range(strlen(chars)):
        needle_mask[<uint8_t>chars[i]] &= ~(1UL << pos)

cdef populate_needle_mask(size_t *needle_mask, char *needle, size_t needle_length,
                          bint ref_wildcards, query_wildcards):
    cdef size_t i
    cdef char c
    if needle_length > (sizeof(size_t) * 8 - 1):
        raise ValueError("The pattern is too long!")
    memset(needle_mask, 0xff, sizeof(size_t) * ASCII_CHAR_COUNT)
    for i in range(needle_length):
        c = needle[i]
        if c == b"A" or c == b"a":
            set_masks(needle_mask, i, "Aa")
            if query_wildcards:
                set_masks(needle_mask, i, "RWMDHVNrwmdhvn")
        elif c == b"C" or c == b"c":
            set_masks(needle_mask, i, "Cc")
            if query_wildcards:
                set_masks(needle_mask, i, "YSMBHVNysmbhvn")
        elif c == b"G" or c == b"g":
            set_masks(needle_mask, i, "Gg")
            if query_wildcards:
                set_masks(needle_mask, i, "RSKBDVNrskbdvn")
        elif c == b"T" or c == b"t":
            set_masks(needle_mask, i, "Tt")
            if query_wildcards:
                set_masks(needle_mask, i, "YWKBDHNywkbdhn")
        elif (c == b"U" or c == b"u") and ref_wildcards:
            set_masks(needle_mask, i, "TtUu")
            if query_wildcards:
                set_masks(needle_mask, i, "YWKBDHNywkbdhn")
        elif (c == b"R" or c == b"r") and ref_wildcards:
            set_masks(needle_mask, i, "AaGg")
            if query_wildcards:
                set_masks(needle_mask, i, "RDVNrdvn")
        elif (c == b"Y" or c == b"y") and ref_wildcards:
            set_masks(needle_mask, i, "CcTt")
            if query_wildcards:
                set_masks(needle_mask, i, "YBHNybhn")
        else:
            needle_mask[<uint8_t>c] &= ~(1UL << i)


cdef char *shift_and_search(char *haystack, size_t haystack_length,
                            size_t *needle_mask, size_t needle_length):
    cdef:
        size_t R = ~1
        size_t i

    if needle_length == 0:
        return haystack

    for i in range(haystack_length):
        # Update the bit array
        R |= needle_mask[<uint8_t>haystack[i]]
        R <<= 1
        if (0 == (R & (1UL << needle_length))):
            return (haystack + i - needle_length) + 1

    return NULL
