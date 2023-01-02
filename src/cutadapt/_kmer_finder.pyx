# cython: profile=False, emit_code_comments=False, language_level=3

from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.string cimport memset, strlen
from cpython.unicode cimport PyUnicode_CheckExact, PyUnicode_GET_LENGTH
from libc.stdint cimport uint8_t

# Dnaio conveniently ensures that all sequences are ASCII only.
DEF BITMASK_INDEX_SIZE = 128

# Make bitmask type definable. size_t is the largest unsigned integer available
# to the machine.
ctypedef size_t bitmask_t

cdef extern from "Python.h":
    void *PyUnicode_DATA(object o)
    bint PyUnicode_IS_COMPACT_ASCII(object o)

ctypedef struct KmerEntry:
    size_t kmer_length
    size_t mask_offset
    ssize_t search_start
    ssize_t search_stop


cdef class KmerFinder:
    """
    Find kmers in strings. To replace the following code:

        kmers_and_positions = [("AGA", -10, None), ("AGCATGA", 0, None)]
        for sequence in sequences:
            for kmer, start, stop in kmers_and_positions:
                if sequence.find(kmer, start, stop) != -1:
                    # do something
                    pass

    This has a lot of python overhead. The following code is equivalent:

        kmers_and_positions = [("AGA", -10, None), ("AGCATGA", 0, None)]
        kmer_finder = KmerFinder(kmers_and_positions)
        for sequence in sequences:
            if kmer_finder.kmers_present(sequence):
                # do something
                pass

    This is more efficient as the kmers_present method can be applied to a lot
    of sequences and all the necessary unpacking for each kmer into C variables
    happens only once.
    """
    cdef:
        KmerEntry *kmer_entries
        bitmask_t *kmer_masks
        size_t number_of_kmers
        readonly object kmers_and_positions
        readonly bint ref_wildcards
        readonly bint query_wildcards

    def __cinit__(self, kmers_and_positions, ref_wildcards=False, query_wildcards=False):
        self.kmer_masks = NULL
        self.kmer_entries = NULL
        self.number_of_kmers = 0
        self.ref_wildcards = ref_wildcards
        self.query_wildcards = query_wildcards
        number_of_entries = len(kmers_and_positions)
        self.kmer_entries = <KmerEntry *>PyMem_Malloc(number_of_entries * sizeof(KmerEntry))
        # for the kmers the NULL bytes also need space.
        self.kmer_masks = <bitmask_t *>PyMem_Malloc(number_of_entries * sizeof(bitmask_t) * BITMASK_INDEX_SIZE)
        self.number_of_kmers = number_of_entries
        cdef size_t mask_offset = 0
        cdef char *kmer_ptr
        cdef Py_ssize_t kmer_length
        for i, (kmer, start, stop) in enumerate(kmers_and_positions):
            if not PyUnicode_CheckExact(kmer):
                raise TypeError(f"Kmer should be a string not {type(kmer)}")
            if not PyUnicode_IS_COMPACT_ASCII(kmer):
                raise ValueError("Only ASCII strings are supported")
            self.kmer_entries[i].search_start  = start
            if stop is None:
                stop = 0
            self.kmer_entries[i].search_stop = stop
            self.kmer_entries[i].mask_offset = mask_offset
            kmer_length = PyUnicode_GET_LENGTH(kmer)
            self.kmer_entries[i].kmer_length = kmer_length
            kmer_ptr = <char *>PyUnicode_DATA(kmer)
            populate_needle_mask(self.kmer_masks + mask_offset, kmer_ptr, kmer_length,
                                 self.ref_wildcards, self.query_wildcards)
            mask_offset += BITMASK_INDEX_SIZE
        self.kmers_and_positions = kmers_and_positions

    def __reduce__(self):
        return KmerFinder, (self.kmers_and_positions, self.ref_wildcards, self.query_wildcards)

    def kmers_present(self, str sequence):
        cdef:
            KmerEntry entry
            size_t i
            size_t kmer_offset
            size_t kmer_length
            ssize_t start, stop
            const bitmask_t *mask_ptr
            const char *search_ptr
            const char *search_result
            ssize_t search_length
        if not PyUnicode_IS_COMPACT_ASCII(sequence):
            raise ValueError("Only ASCII strings are supported")
        cdef const char *seq = <char *>PyUnicode_DATA(sequence)
        cdef Py_ssize_t seq_length = PyUnicode_GET_LENGTH(sequence)
        for i in range(self.number_of_kmers):
            entry = self.kmer_entries[i]
            start = entry.search_start 
            stop = entry.search_stop
            if start < 0:
                start = seq_length + start
                if start < 0:
                    start = 0
            elif start > seq_length:
                continue
            if stop < 0:
                step = seq_length + stop
                if stop <= 0:  # No need to search
                    continue
            elif stop == 0:
                stop = seq_length
            search_length = stop - start
            if search_length <= 0:
                continue
            search_ptr = seq + start
            kmer_length = entry.kmer_length
            mask_ptr = self.kmer_masks + entry.mask_offset
            search_result = shift_or_search(search_ptr, search_length,
                                            mask_ptr, kmer_length)
            if search_result:
                return True
        return False

    def __dealloc__(self):
        PyMem_Free(self.kmer_masks)
        PyMem_Free(self.kmer_entries)


cdef void set_masks(bitmask_t *needle_mask, size_t pos, const char *chars):
    cdef char c
    cdef size_t i
    for i in range(strlen(chars)):
        needle_mask[<uint8_t>chars[i]] &= ~(<bitmask_t>1ULL << pos)

cdef populate_needle_mask(bitmask_t *needle_mask, const char *needle, size_t needle_length,
                          bint ref_wildcards, bint query_wildcards):
    cdef size_t i
    cdef char c
    cdef uint8_t j
    if needle_length > (sizeof(bitmask_t) * 8 - 1):
        raise ValueError("The pattern is too long!")
    memset(needle_mask, 0xff, sizeof(bitmask_t) * BITMASK_INDEX_SIZE)
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
                set_masks(needle_mask, i, "RSWKMBDHVNrswkmvdhvn")
        elif (c == b"Y" or c == b"y") and ref_wildcards:
            set_masks(needle_mask, i, "CcTtUu")
            if query_wildcards:
                set_masks(needle_mask, i, "YSWKMBDHVNyswkmbdhvn")
        elif (c == b"S" or c == b"s") and ref_wildcards:
            set_masks(needle_mask, i, "GgCc")
            if query_wildcards:
                set_masks(needle_mask, i, "YRSKMBDHVNyrskmbdhvn")
        elif (c == b"W" or c == b"w") and ref_wildcards:
            set_masks(needle_mask, i, "AaTtUu")
            if query_wildcards:
                set_masks(needle_mask, i, "YRWKMBDHVNyrskmbdhvn")
        elif (c == b"K" or c == b"k") and ref_wildcards:
            set_masks(needle_mask, i, "GgTtUu")
            if query_wildcards:
                set_masks(needle_mask, i, "YRWSKBDHVNyrskmbdhvn")
        elif (c == b"M" or c == b"m") and ref_wildcards:
            set_masks(needle_mask, i, "AaCc")
            if query_wildcards:
                set_masks(needle_mask, i, "YRWSMBDHVNyrskmbdhvn")
        elif (c == b"B" or  c == b"b") and ref_wildcards:
            set_masks(needle_mask, i, "CcGgTtUu")
            if query_wildcards:
                set_masks(needle_mask, i, "RYSWKMBDHVNryswkmbdhvn")
        elif (c == b"D" or c == b"d") and ref_wildcards:
            set_masks(needle_mask, i, "AaGgTtUu")
            if query_wildcards:
                set_masks(needle_mask, i, "RYSWKMBDHVNryswkmbdhvn")
        elif (c == b"H" or c == b"h") and ref_wildcards:
            set_masks(needle_mask, i, "AaCcTtUu")
            if query_wildcards:
                set_masks(needle_mask, i, "RYSWKMBDHVNryswkmbdhvn")
        elif (c == b"V" or c == b"v") and ref_wildcards:
            set_masks(needle_mask, i, "AaCcGg")
            if query_wildcards:
                set_masks(needle_mask, i, "RYSWKMBDHVNryswkmbdhvn")
        elif (c == b"N" or c == b"n") and ref_wildcards:
            if query_wildcards:  # Proper IUPAC matching
                set_masks(needle_mask, i, "ACGTURYSWKMBDHVNacgturyswkmbdhvn")
            else:  # N matches literally everything except \00
                for j in range(1,128):
                    needle_mask[j] &= ~(<bitmask_t>1ULL << i)
        else:
            needle_mask[<uint8_t>c] &= ~(<bitmask_t>1ULL << i)


cdef const char *shift_or_search(const char *haystack, size_t haystack_length,
                           const bitmask_t *needle_mask, size_t needle_length):
    cdef:
        bitmask_t R = ~1
        size_t i

    if needle_length == 0:
        return haystack

    for i in range(haystack_length):
        # Update the bit array
        R |= needle_mask[<uint8_t>haystack[i]]
        R <<= 1
        if (0 == (R & (<bitmask_t>1ULL << needle_length))):
            return (haystack + i - needle_length) + 1

    return NULL
