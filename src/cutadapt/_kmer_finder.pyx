# cython: profile=False, emit_code_comments=False, language_level=3

from cpython.mem cimport PyMem_Realloc, PyMem_Free
from libc.string cimport memcpy, memset, strlen

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

ctypedef struct KmerSearchEntry:
    size_t mask_offset
    ssize_t search_start
    ssize_t search_stop
    bitmask_t zero_mask
    bitmask_t found_mask


cdef class KmerFinder:
    """
    Find kmers in strings. Allows case-independent and IUPAC matching.
    ``ref_wildcards=True`` allows IUPAC characters in the kmer sequences.
    ``query_wildcards=True`` allows IUPAC characters in the sequences fed to
    the ``kmers_present`` method.

    Replaces the following code:

        kmers_and_positions = [("AGA", -10, None), ("AGCATGA", 0, None)]
        for sequence in sequences:
            for kmer, start, stop in kmers_and_positions:
                if sequence.find(kmer, start, stop) != -1:
                    # do something
                    pass

    This has a lot of python overhead. The following code accomplishes the
    same task and allows for case-independent matching:

        positions_and_kmers = [(-10, None, ["AGA"]), (0, None, ["AGCATGA"])]
        kmer_finder = KmerFinder(positions_and_kmers)
        for sequence in sequences:
            if kmer_finder.kmers_present(sequence):
                # do something
                continue

    This is more efficient as the kmers_present method can be applied to a lot
    of sequences and all the necessary unpacking for each kmer into C variables
    happens only once.
    Note that multiple kmers can be given per position. Kmerfinder finds
    all of these simultaneously using a multiple pattern matching algorithm.
    """
    cdef:
        KmerSearchEntry *search_entries
        bitmask_t *search_masks
        size_t number_of_searches
        readonly object positions_and_kmers
        readonly bint ref_wildcards
        readonly bint query_wildcards

    def __cinit__(self, positions_and_kmers, ref_wildcards=False, query_wildcards=False):
        cdef char[64] search_word
        self.search_masks = NULL
        self.search_entries = NULL
        self.number_of_searches = 0
        self.ref_wildcards = ref_wildcards
        self.query_wildcards = query_wildcards
        self.number_of_searches = 0
        cdef size_t mask_offset = 0
        cdef char *kmer_ptr
        cdef size_t offset
        cdef bitmask_t zero_mask, found_mask
        cdef Py_ssize_t kmer_length
        # The maximum length of words + null bytes that we can store in
        # the bitmask array
        cdef size_t max_total_length = sizeof(bitmask_t) * 8
        # The maximum length of a word. Since word_length + 1 bits are needed to search.
        cdef ssize_t max_word_length = max_total_length - 1

        for (start, stop, kmers) in positions_and_kmers:
            index = 0 
            while index < len(kmers):
                memset(search_word, 0, 64)
                offset = 0
                zero_mask = ~0
                found_mask = 0
                # Run an inner loop in case all words combined are larger than
                # the maximum bitmask size. In that case we create multiple
                # bitmasks to hold all the words.
                while index < len(kmers):
                    kmer = kmers[index]
                    if not PyUnicode_CheckExact(kmer):
                        raise TypeError(f"Kmer should be a string not {type(kmer)}")
                    if not PyUnicode_IS_COMPACT_ASCII(kmer):
                        raise ValueError("Only ASCII strings are supported")
                    kmer_length = PyUnicode_GET_LENGTH(kmer)
                    if kmer_length > max_word_length:
                        raise ValueError(f"{kmer} of length {kmer_length} is longer "
                                         f"than the maximum of {max_word_length}.")
                    if (offset + kmer_length + 1) > max_total_length:
                        break
                    zero_mask ^= <bitmask_t>1ULL << offset
                    kmer_ptr = <char *> PyUnicode_DATA(kmer)
                    memcpy(search_word + offset, kmer_ptr, kmer_length)
                    search_word[offset + kmer_length] = 0
                    found_mask |= <bitmask_t>1ULL << (offset + kmer_length)
                    offset = offset + kmer_length + 1
                    index += 1
                i = self.number_of_searches  # Save the index position for the mask and entry
                self.number_of_searches += 1
                self.search_entries = <KmerSearchEntry *>PyMem_Realloc(self.search_entries, self.number_of_searches * sizeof(KmerSearchEntry))
                self.search_masks = <bitmask_t *>PyMem_Realloc(self.search_masks, self.number_of_searches * sizeof(bitmask_t) * BITMASK_INDEX_SIZE)
                mask_offset = i * BITMASK_INDEX_SIZE
                self.search_entries[i].search_start  = start
                if stop is None:
                    stop = 0
                self.search_entries[i].search_stop = stop
                self.search_entries[i].mask_offset = mask_offset
                self.search_entries[i].zero_mask = zero_mask
                self.search_entries[i].found_mask = found_mask
                # Offset -1 because we don't count the last NULL byte
                populate_needle_mask(self.search_masks + mask_offset, search_word, offset - 1,
                                     self.ref_wildcards, self.query_wildcards)
        self.positions_and_kmers = positions_and_kmers

    def __reduce__(self):
        return KmerFinder, (self.positions_and_kmers, self.ref_wildcards, self.query_wildcards)

    def kmers_present(self, str sequence):
        cdef:
            KmerSearchEntry entry
            size_t i
            size_t kmer_offset
            bitmask_t zero_mask
            bitmask_t found_mask
            ssize_t start, stop
            const bitmask_t *mask_ptr
            const char *search_ptr
            bint search_result
            ssize_t search_length
        if not PyUnicode_IS_COMPACT_ASCII(sequence):
            raise ValueError("Only ASCII strings are supported")
        cdef const char *seq = <char *>PyUnicode_DATA(sequence)
        cdef Py_ssize_t seq_length = PyUnicode_GET_LENGTH(sequence)
        for i in range(self.number_of_searches):
            entry = self.search_entries[i]
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
            zero_mask = entry.zero_mask
            found_mask = entry.found_mask
            mask_ptr = self.search_masks + entry.mask_offset
            search_result = shift_or_multiple_is_present(
                search_ptr, search_length, mask_ptr, zero_mask, found_mask)
            if search_result:
                return True
        return False

    def __dealloc__(self):
        PyMem_Free(self.search_masks)
        PyMem_Free(self.search_entries)


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
        if c == 0:
            continue
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
            if ref_wildcards:
                set_masks(needle_mask, i, "Uu")
            if query_wildcards:
                set_masks(needle_mask, i, "UYWKBDHNuywkbdhn")
        elif (c == b"U" or c == b"u") and (ref_wildcards or query_wildcards):
            set_masks(needle_mask, i, "TtUu")
            if query_wildcards:
                set_masks(needle_mask, i, "YWKBDHNywkbdhn")
        elif (c == b"R" or c == b"r") and ref_wildcards:
            set_masks(needle_mask, i, "AaGg")
            if query_wildcards:
                set_masks(needle_mask, i, "RSWKMBDHVNrswkmbdhvn")
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
                set_masks(needle_mask, i, "YRWKMBDHVNyrwkmbdhvn")
        elif (c == b"K" or c == b"k") and ref_wildcards:
            set_masks(needle_mask, i, "GgTtUu")
            if query_wildcards:
                set_masks(needle_mask, i, "YRWSKBDHVNyrwskbdhvn")
        elif (c == b"M" or c == b"m") and ref_wildcards:
            set_masks(needle_mask, i, "AaCc")
            if query_wildcards:
                set_masks(needle_mask, i, "YRWSMBDHVNyrwsmbdhvn")
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
        elif query_wildcards and not ref_wildcards:
            # All non-AGCT characters match to N
            set_masks(needle_mask, i, "Nn")
        elif ref_wildcards:
            # ref and query wildcards are True. Perform proper IUPAC matching.
            # All unknown characters do not match.
            pass
        else:
            if chr(c).isalpha():
                bothcase = chr(c).lower() + chr(c).upper()
                set_masks(needle_mask, i, bothcase.encode("ascii"))
            else:
                needle_mask[<uint8_t>c] &= ~(<bitmask_t>1ULL << i)


cdef bint shift_or_multiple_is_present(
    const char *haystack,
    size_t haystack_length,
    const bitmask_t *needle_mask,
    bitmask_t zero_mask,
    bitmask_t found_mask):
    cdef:
        bitmask_t R = zero_mask
        size_t i

    for i in range(haystack_length):
        # Update the bit array
        R |= needle_mask[<uint8_t>haystack[i]]
        R <<= 1
        R &= zero_mask
        if ((R & found_mask) != found_mask):
            return True

    return False
