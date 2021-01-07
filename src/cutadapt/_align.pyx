# cython: profile=False, emit_code_comments=False, language_level=3
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc

DEF MATCH_SCORE = 16384

# structure for a DP matrix entry
ctypedef struct _Entry:
    int cost
    int score    # score for this alignment (mostly keeps track of matches)
    int origin   # where the alignment originated: negative for positions within seq1, positive for pos. within seq2


ctypedef struct _Match:
    int origin
    int cost
    int matches
    int score
    int ref_stop
    int query_stop


def _acgt_table():
    """
    Return a translation table that maps A, C, G, T characters to the lower
    four bits of a byte. Other characters (including possibly IUPAC characters)
    are mapped to the most significant bit (0x80).

    Lowercase versions are also translated, and U is treated the same as T.
    """
    d = dict(A=1, C=2, G=4, T=8, U=8)
    t = bytearray([0x80]) * 256
    for c, v in d.items():
        t[ord(c)] = v
        t[ord(c.lower())] = v
    return bytes(t)


def _iupac_table():
    """
    Return a translation table for IUPAC characters.

    The table maps ASCII-encoded IUPAC nucleotide characters to bytes in which
    the four least significant bits are used to represent one nucleotide each.

    For the "N" wildcard, additionally the most significant bit is set (0x80),
    which allows it to match characters that are not A, C, G or T if _acgt_table
    was used to encode them.

    Whether two encoded characters x and y match can then be checked with the
    expression "x & y != 0".
    """
    A = 1
    C = 2
    G = 4
    T = 8
    iupac = dict(
        X=0,
        A=A,
        C=C,
        G=G,
        T=T,
        U=T,
        R=A|G,
        Y=C|T,
        S=G|C,
        W=A|T,
        K=G|T,
        M=A|C,
        B=C|G|T,
        D=A|G|T,
        H=A|C|T,
        V=A|C|G,
        N=A|C|G|T + 0x80,
    )
    t = bytearray(b'\0') * 256
    for c, v in iupac.items():
        t[ord(c)] = v
        t[ord(c.lower())] = v
    return bytes(t)


cdef:
    bytes ACGT_TABLE = _acgt_table()
    bytes IUPAC_TABLE = _iupac_table()


class DPMatrix:
    """
    Representation of the dynamic-programming matrix.

    This is used only when debugging is enabled in the Aligner class since the
    matrix is normally not stored in full.

    Entries in the matrix may be None, in which case that value was not
    computed.
    """
    def __init__(self, reference, query):
        m = len(reference)
        n = len(query)
        self._rows = [ [None] * (n+1) for _ in range(m + 1) ]
        self.reference = reference
        self.query = query

    def set_entry(self, int i, int j, cost):
        """
        Set an entry in the dynamic programming matrix.
        """
        self._rows[i][j] = cost

    def __str__(self):
        """
        Return a representation of the matrix as a string.
        """
        rows = ['     ' + ' '.join(c.rjust(2) for c in self.query)]
        for c, row in zip(' ' + self.reference, self._rows):
            r = c + ' ' + ' '.join('  ' if v is None else '{:2d}'.format(v) for v in row)
            rows.append(r)
        return '\n'.join(rows)


cdef int score_to_matches(int score) nogil:
    return (score + MATCH_SCORE - 1) // MATCH_SCORE


cdef class Aligner:
    """
    Find a full or partial occurrence of a query string in a reference string
    allowing errors (mismatches, insertions, deletions).

    By default, unit costs are used, meaning that mismatches, insertions and
    deletions are counted as one error (edit distance).

    Semi-global alignments allow skipping a suffix and/or prefix of query or
    reference at no cost. Combining semi-global alignment with edit distance is
    a bit unusual because the trivial “optimal” solution at edit distance 0
    would be to skip all of the reference and all of the query, like this:

        REFERENCE-----
        ---------QUERY

    Conceptually, the algorithm used here instead tests all possible overlaps
    between the two sequences and chooses the overlap which maximizes the
    number of matches in the overlapping part, while the error rate must not
    go above a threshold.

    TODO working here

    To allow skipping of a prefix of string1 at no cost, set the
    START_IN_REFERENCE flag.
    To allow skipping of a prefix of string2 at no cost, set the
    START_IN_QUERY flag.
    If both are set, a prefix of string1 or of string1 is skipped,
    never both.
    Similarly, set STOP_IN_REFERENCE and STOP_IN_QUERY to
    allow skipping of suffixes of string1 or string2. Again, when both
    flags are set, never suffixes in both strings are skipped.
    If all flags are set, this results in standard semiglobal alignment.

    The skipped parts are described with two intervals (start1, stop1),
    (start2, stop2).

    For example, an optimal semiglobal alignment of SISSI and MISSISSIPPI looks like this:

    ---SISSI---
    MISSISSIPPI

    start1, stop1 = 0, 5
    start2, stop2 = 3, 8
    (with zero errors)

    The aligned parts are string1[start1:stop1] and string2[start2:stop2].

    The error rate is: errors / length where length is (stop1 - start1).

    An optimal alignment fulfills all of these criteria:

    - its error_rate is at most max_error_rate
    - Among those alignments with error_rate <= max_error_rate, the alignment contains
      a maximal number of matches (there is no alignment with more matches).
    - If there are multiple alignments with the same no. of matches, then one that
      has minimal no. of errors is chosen.
    - If there are still multiple candidates, choose the alignment that starts at the
      leftmost position within the read.

    The alignment itself is not returned, only the tuple
    (start1, stop1, start2, stop2, matches, errors), where the first four fields have the
    meaning as described, matches is the number of matches and errors is the number of
    errors in the alignment.

    It is always the case that at least one of start1 and start2 is zero.

    IUPAC wildcard characters can be allowed in the reference and the query
    by setting the appropriate flags.

    If neither flag is set, the full ASCII alphabet is used for comparison.
    If any of the flags is set, all non-IUPAC characters in the sequences
    compare as 'not equal'.
    """
    cdef:
        int m
        _Entry* column  # one column of the DP matrix
        double max_error_rate
        bint start_in_reference
        bint start_in_query
        bint stop_in_reference
        bint stop_in_query
        int _insertion_cost
        int _deletion_cost
        int _min_overlap
        bint wildcard_ref
        bint wildcard_query
        bint debug
        object _dpmatrix
        str reference  # reference as set by the user (as str)
        bytes _reference  # internal, bytes version of reference (possibly translated to a non-ASCII representation)
        readonly int effective_length
        int* n_counts  # n_counts[i] == number of N characters in reference[:i]

    def __cinit__(
        self,
        str reference,
        double max_error_rate,
        int flags=15,
        bint wildcard_ref=False,
        bint wildcard_query=False,
        int indel_cost=1,
        int min_overlap=1,
    ):
        self.max_error_rate = max_error_rate
        self.start_in_reference = flags & 1
        self.start_in_query = flags & 2
        self.stop_in_reference = flags & 4
        self.stop_in_query = flags & 8
        self.wildcard_ref = wildcard_ref
        self.wildcard_query = wildcard_query
        self._set_reference(reference)
        self._min_overlap = min_overlap
        self.debug = False
        self._dpmatrix = None
        if indel_cost < 1:
            raise ValueError('indel_cost must be at least 1')
        self._insertion_cost = indel_cost
        self._deletion_cost = indel_cost

    def _compute_flags(self):
        cdef int flags = 0
        if self.start_in_reference:
            flags |= 1
        if self.start_in_query:
            flags |= 2
        if self.stop_in_reference:
            flags |= 4
        if self.stop_in_query:
            flags |= 8
        return flags

    def __reduce__(self):
        return (Aligner, (self.reference, self.max_error_rate, self._compute_flags(), self.wildcard_ref, self.wildcard_query, self._insertion_cost, self._min_overlap))

    def __repr__(self):
        return (
            f"Aligner(reference='{self.reference}', max_error_rate={self.max_error_rate}, "
            f"flags={self._compute_flags()}, wildcard_ref={self.wildcard_ref}, "
            f"wildcard_query={self.wildcard_query}, indel_cost={self._insertion_cost}, "
            f"min_overlap={self._min_overlap})"
        )

    def _set_reference(self, str reference):
        mem = <_Entry*> PyMem_Realloc(self.column, (len(reference) + 1) * sizeof(_Entry))
        if not mem:
            raise MemoryError()
        mem_nc = <int*> PyMem_Realloc(self.n_counts, (len(reference) + 1) * sizeof(int))
        if not mem_nc:
            raise MemoryError()
        self.column = mem
        self.n_counts = mem_nc
        self._reference = reference.encode('ascii')
        self.m = len(reference)
        self.effective_length = self.m
        n_count = 0
        for i in range(self.m):
            self.n_counts[i] = n_count
            if reference[i] == 'n' or reference[i] == 'N':
                n_count += 1
        self.n_counts[self.m] = n_count
        assert self.n_counts[self.m] == reference.count('N') + reference.count('n')
        if self.wildcard_ref:
            self.effective_length = self.m - self.n_counts[self.m]
            if self.effective_length == 0:
                raise ValueError("Cannot have only N wildcards in the sequence")
            self._reference = self._reference.translate(IUPAC_TABLE)
        elif self.wildcard_query:
            self._reference = self._reference.translate(ACGT_TABLE)
        self.reference = reference

    property dpmatrix:
        """
        The dynamic programming matrix as a DPMatrix object. This attribute is
        usually None, unless debugging has been enabled with enable_debug().
        """
        def __get__(self):
            return self._dpmatrix

    def enable_debug(self):
        """
        Store the dynamic programming matrix while running the locate() method
        and make it available in the .dpmatrix attribute.
        """
        self.debug = True

    def locate(self, str query):
        """
        locate(query) -> (refstart, refstop, querystart, querystop, matches, errors)

        Find the query within the reference associated with this aligner. The
        intervals (querystart, querystop) and (refstart, refstop) give the
        location of the match.

        That is, the substrings query[querystart:querystop] and
        self.reference[refstart:refstop] were found to align best to each other,
        with the given number of matches and the given number of errors.

        The alignment itself is not returned.
        """
        cdef:
            char* s1 = self._reference
            bytes query_bytes = query.encode('ascii')
            char* s2
            int m = self.m
            int n = len(query)
            _Entry* column = self.column  # Current column of the DP matrix
            double max_error_rate = self.max_error_rate
            bint stop_in_query = self.stop_in_query
            bint compare_ascii = False

        if self.wildcard_query:
            query_bytes = query_bytes.translate(IUPAC_TABLE)
        elif self.wildcard_ref:
            query_bytes = query_bytes.translate(ACGT_TABLE)
        else:
            # TODO Adding the .upper() increases overall runtime slightly even
            # when I remove the .upper() from Adapter.match_to().
            query_bytes = query_bytes.upper()
            compare_ascii = True
        s2 = query_bytes
        """
        DP Matrix:
                   query (j)
                 ----------> n
                |
        ref (i) |
                |
                V
               m
        """
        cdef int i, j

        # maximum no. of errors
        cdef int k = <int> (max_error_rate * m)

        # Determine largest and smallest column we need to compute
        cdef int max_n = n
        cdef int min_n = 0
        if not self.start_in_query:
            # costs can only get worse after column m
            max_n = min(n, m + k)
        if not self.stop_in_query:
            min_n = max(0, n - m - k)

        # Fill column min_n.
        #
        # Four cases:
        # not startin1, not startin2: c(i,j) = max(i,j); origin(i, j) = 0
        #     startin1, not startin2: c(i,j) = j       ; origin(i, j) = min(0, j - i)
        # not startin1,     startin2: c(i,j) = i       ; origin(i, j) =
        #     startin1,     startin2: c(i,j) = min(i,j)

        # TODO (later)
        # fill out columns only until 'last'
        if not self.start_in_reference and not self.start_in_query:
            for i in range(m + 1):
                column[i].score = 0
                column[i].cost = max(i, min_n) * self._insertion_cost
                column[i].origin = 0
        elif self.start_in_reference and not self.start_in_query:
            for i in range(m + 1):
                column[i].score = 0
                column[i].cost = min_n * self._insertion_cost
                column[i].origin = min(0, min_n - i)
        elif not self.start_in_reference and self.start_in_query:
            for i in range(m + 1):
                column[i].score = 0
                column[i].cost = i * self._insertion_cost
                column[i].origin = max(0, min_n - i)
        else:
            for i in range(m + 1):
                column[i].score = 0
                column[i].cost = min(i, min_n) * self._insertion_cost
                column[i].origin = min_n - i

        if self.debug:
            self._dpmatrix = DPMatrix(self.reference, query)
            for i in range(m + 1):
                self._dpmatrix.set_entry(i, min_n, column[i].cost)
        cdef _Match best
        best.ref_stop = m
        best.query_stop = n
        best.cost = m + n + 1
        best.origin = 0
        best.score = 0
        best.matches = 0

        # Ukkonen's trick: index of the last cell that is at most k
        cdef int last = min(m, k + 1)
        if self.start_in_reference:
            last = m

        cdef:
            int cost_diag
            int cost_deletion
            int cost_insertion
            int origin, cost, matches, score
            int length
            int ref_start
            int cur_effective_length
            bint characters_equal
            bint is_acceptable
            # We keep only a single column of the DP matrix in memory.
            # To access the diagonal cell to the upper left,
            # we store it here before overwriting it.
            _Entry diag_entry

        with nogil:
            # iterate over columns
            for j in range(min_n + 1, max_n + 1):
                # remember first entry before overwriting
                diag_entry = column[0]

                # fill in first entry in this column
                if self.start_in_query:
                    column[0].origin = j
                else:
                    column[0].cost = j * self._insertion_cost
                for i in range(1, last + 1):
                    if compare_ascii:
                        characters_equal = (s1[i-1] == s2[j-1])
                    else:
                        characters_equal = (s1[i-1] & s2[j-1]) != 0
                    if characters_equal:
                        # If the characters match, we can skip computing costs for
                        # insertion and deletion as they are at least as high.
                        cost = diag_entry.cost
                        origin = diag_entry.origin
                        # Among the optimal alignments whose edit distance is within the
                        # maximum allowed error rate, we prefer those that have the
                        # maximum number of matches. And among those, we want to
                        # prefer those with mismatches over those with indels.
                        # To achieve this, scores are as follows
                        # - match: MATCH_SCORE
                        # - mismatch: 0
                        # - indel: -1
                        score = diag_entry.score + MATCH_SCORE
                    else:
                        # Characters do not match.
                        cost_diag = diag_entry.cost + 1
                        cost_deletion = column[i].cost + self._deletion_cost
                        cost_insertion = column[i-1].cost + self._insertion_cost

                        if cost_diag <= cost_deletion and cost_diag <= cost_insertion:
                            # MISMATCH
                            cost = cost_diag
                            origin = diag_entry.origin
                            score = diag_entry.score
                        elif cost_insertion <= cost_deletion:
                            # INSERTION
                            cost = cost_insertion
                            origin = column[i-1].origin
                            # penalize insertions slightly
                            score = column[i-1].score - 1
                        else:
                            # DELETION
                            cost = cost_deletion
                            origin = column[i].origin
                            # penalize deletions slightly
                            score = column[i].score - 1

                    # Remember the current cell for next iteration
                    diag_entry = column[i]

                    column[i].cost = cost
                    column[i].origin = origin
                    column[i].score = score
                if self.debug:
                    with gil:
                        for i in range(last + 1):
                            self._dpmatrix.set_entry(i, j, column[i].cost)
                while last >= 0 and column[last].cost > k:
                    last -= 1
                # last can be -1 here, but will be incremented next.
                # TODO if last is -1, can we stop searching?
                if last < m:
                    last += 1
                elif stop_in_query:
                    # Found a match. If requested, find best match in last row.
                    # length of the aligned part of the reference
                    length = m + min(column[m].origin, 0)
                    cur_effective_length = length
                    if self.wildcard_ref:
                        if length < m:
                            # Recompute effective length so that it only takes into
                            # account the matching suffix of the reference
                            cur_effective_length = length - self.n_counts[length]
                        else:
                            cur_effective_length = self.effective_length
                    cost = column[m].cost
                    score = column[m].score
                    matches = score_to_matches(score)
                    origin = column[m].origin
                    is_acceptable = (
                        length >= self._min_overlap
                        and cost <= cur_effective_length * max_error_rate
                    )
                    if is_acceptable and (
                        # no best match recorded so far
                        (best.cost == m + n + 1)
                        # same start position, use score to judge whether this is better
                        or (origin == best.origin and score > best.score)
                        # different start position, only matches count
                        or (matches > best.matches or (matches == best.matches and cost < best.cost))
                    ):
                        # The case "matches == best.matches and cost < best.cost" applies
                        # when the query contains the reference twice, and where the right
                        # occurrence contains fewer errors
                        best.score = score
                        best.matches = matches
                        best.cost = cost
                        best.origin = origin
                        best.ref_stop = m
                        best.query_stop = j
                        if cost == 0 and matches == m:
                            # exact match, stop early
                            break
                # column finished

        if max_n == n:
            first_i = 0 if self.stop_in_reference else m
            # search in last column # TODO last?
            for i in range(first_i, m+1):
                length = i + min(column[i].origin, 0)
                cost = column[i].cost
                score = column[i].score
                matches = score_to_matches(score)
                if self.wildcard_ref:
                    if length < m:
                        # Recompute effective length so that it only takes into
                        # account the matching part of the reference
                        ref_start = -min(column[i].origin, 0)
                        assert 0 <= ref_start <= m
                        cur_effective_length = length - (self.n_counts[i] - self.n_counts[ref_start])
                    else:
                        cur_effective_length = self.effective_length
                else:
                    cur_effective_length = length
                assert 0 <= cur_effective_length and cur_effective_length <= length
                assert cur_effective_length <= self.effective_length

                is_acceptable = (
                    length >= self._min_overlap
                    and cost <= cur_effective_length * max_error_rate
                )
                if is_acceptable and (
                    # no best match recorded so far
                    (best.cost == m + n + 1)
                    # same start position, use score to judge whether this is better
                    or (origin == best.origin and score > best.score)
                    # different start position, only matches count
                    or (matches > best.matches or (matches == best.matches and cost < best.cost))
                ):
                    best.matches = matches
                    best.score = score
                    best.cost = cost
                    best.origin = column[i].origin
                    best.ref_stop = i
                    best.query_stop = n
        if best.cost == m + n + 1:
            # best.cost was initialized with this value.
            # If it is unchanged, no alignment was found that has
            # an error rate within the allowed range.
            return None

        cdef int start1, start2
        if best.origin >= 0:
            start1 = 0
            start2 = best.origin
        else:
            start1 = -best.origin
            start2 = 0

        return (start1, best.ref_stop, start2, best.query_stop, best.matches, best.cost)

    def __dealloc__(self):
        PyMem_Free(self.column)
        PyMem_Free(self.n_counts)


cdef class PrefixComparer:
    """
    A version of the Aligner that is specialized in the following way:

    - it does not allow indels
    - it allows only 5' anchored adapters

    This is a separate class, not simply a function, in order to be able
    to cache the reference (avoiding to convert it from str to bytes on
    every invocation)
    """
    cdef:
        bytes reference
        bint wildcard_ref
        bint wildcard_query
        int m
        int max_k  # max. number of errors
        readonly int effective_length
        int min_overlap

    # __init__ instead of __cinit__ because we need to override this in SuffixComparer
    def __init__(
        self,
        str reference,
        double max_error_rate,
        bint wildcard_ref=False,
        bint wildcard_query=False,
        int min_overlap=1,
    ):
        self.wildcard_ref = wildcard_ref
        self.wildcard_query = wildcard_query
        self.m = len(reference)
        self.effective_length = self.m
        if self.wildcard_ref:
            self.effective_length -= reference.count('N') - reference.count('n')
            if self.effective_length == 0:
                raise ValueError("Cannot have only N wildcards in the sequence")
        if not (0 <= max_error_rate <= 1.):
            raise ValueError("max_error_rate must be between 0 and 1")
        self.max_k = int(max_error_rate * self.effective_length)
        self.reference = reference.encode('ascii').upper()
        if min_overlap < 1:
            raise ValueError("min_overlap must be at least 1")
        self.min_overlap = min_overlap
        if self.wildcard_ref:
            self.reference = self.reference.translate(IUPAC_TABLE)
        elif self.wildcard_query:
            self.reference = self.reference.translate(ACGT_TABLE)

    def __repr__(self):
        return "{}(reference={!r}, max_k={}, wildcard_ref={}, "\
            "wildcard_query={})".format(
                self.__class__.__name__,
                self.reference, self.max_k, self.wildcard_ref,
                self.wildcard_query)

    def locate(self, str query):
        """
        Find out whether one string is the prefix of the other one, allowing
        IUPAC wildcards in ref and/or query if the appropriate flag is set.

        This is used to find an anchored 5' adapter (type 'FRONT') in the 'no indels' mode.
        This is very simple as only the number of errors needs to be counted.

        This function returns a tuple compatible with what Aligner.locate outputs.
        """
        cdef:
            bytes query_bytes = query.encode('ascii')
            char* r_ptr = self.reference
            char* q_ptr
            int i, matches = 0
            int n = len(query_bytes)
            int length = min(self.m, n)
            bint compare_ascii = False
            int errors

        if self.wildcard_query:
            query_bytes = query_bytes.translate(IUPAC_TABLE)
        elif self.wildcard_ref:
            query_bytes = query_bytes.translate(ACGT_TABLE)
        else:
            query_bytes = query_bytes.upper()
            compare_ascii = True
        q_ptr = query_bytes

        if compare_ascii:
            for i in range(length):
                if r_ptr[i] == q_ptr[i]:
                    matches += 1
        else:
            for i in range(length):
                if (r_ptr[i] & q_ptr[i]) != 0:
                    matches += 1

        errors = length - matches
        if errors > self.max_k or length < self.min_overlap:
            return None
        return (0, length, 0, length, matches, length - matches)


cdef class SuffixComparer(PrefixComparer):

    def __init__(
        self,
        str reference,
        double max_error_rate,
        bint wildcard_ref=False,
        bint wildcard_query=False,
        int min_overlap=1,
    ):
        super().__init__(reference[::-1], max_error_rate, wildcard_ref, wildcard_query, min_overlap)

    def locate(self, str query):
        cdef int n = len(query)
        result = super().locate(query[::-1])
        if result is None:
            return None
        _, length, _, _, matches, errors = result
        return (self.m - length, self.m, n - length, n, matches, errors)
