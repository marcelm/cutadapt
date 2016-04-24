# cython: profile=False, emit_code_comments=False
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc

DEF START_WITHIN_SEQ1 = 1
DEF START_WITHIN_SEQ2 = 2
DEF STOP_WITHIN_SEQ1 = 4
DEF STOP_WITHIN_SEQ2 = 8
DEF SEMIGLOBAL = 15

# structure for a DP matrix entry
ctypedef struct _Entry:
	int cost
	int matches  # no. of matches in this alignment
	int origin   # where the alignment originated: negative for positions within seq1, positive for pos. within seq2


ctypedef struct _Match:
	int origin
	int cost
	int matches
	int ref_stop
	int query_stop


def _acgt_table():
	"""
	Return a translation table that maps A, C, G, T characters to the lower
	four bits of a byte. Other characters (including possibly IUPAC characters)
	are mapped to zero.

	Lowercase versions are also translated, and U is treated the same as T.
	"""
	d = dict(A=1, C=2, G=4, T=8, U=8)
	t = bytearray(b'\0') * 256
	for c, v in d.items():
		t[ord(c)] = v
		t[ord(c.lower())] = v
	return bytes(t)


def _iupac_table():
	"""
	Return a translation table for IUPAC characters.

	The table maps ASCII-encoded IUPAC nucleotide characters to bytes in which
	the four least significant bits are used to represent one nucleotide each.

	Whether two characters x and y match can then be checked with the
	expression "x & y != 0".
	"""
	A = 1
	C = 2
	G = 4
	T = 8
	d = dict(
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
		N=A|C|G|T
	)
	t = bytearray(b'\0') * 256
	for c, v in d.items():
		t[ord(c)] = v
		t[ord(c.lower())] = v
	return bytes(t)


cdef bytes ACGT_TABLE = _acgt_table()
cdef bytes IUPAC_TABLE = _iupac_table()


class DPMatrix:
	"""
	Representation of the dynamic-programming matrix.

	This used only when debugging is enabled in the Aligner class since the
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
			r = c + ' ' + ' '.join('  ' if v is None else '{0:2d}'.format(v) for v in row)
			rows.append(r)
		return '\n'.join(rows)


cdef class Aligner:
	"""
	TODO documentation still uses s1 (reference) and s2 (query).

	Locate one string within another by computing an optimal semiglobal
	alignment between string1 and string2.

	The alignment uses unit costs, which means that mismatches, insertions and deletions are
	counted as one error.

	flags is a bitwise 'or' of the allowed flags.
	To allow skipping of a prefix of string1 at no cost, set the
	START_WITHIN_SEQ1 flag.
	To allow skipping of a prefix of string2 at no cost, set the
	START_WITHIN_SEQ2 flag.
	If both are set, a prefix of string1 or of string1 is skipped,
	never both.
	Similarly, set STOP_WITHIN_SEQ1 and STOP_WITHIN_SEQ2 to
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
	cdef int m
	cdef _Entry* column  # one column of the DP matrix
	cdef double max_error_rate
	cdef int flags
	cdef int _insertion_cost
	cdef int _deletion_cost
	cdef int _min_overlap
	cdef bint wildcard_ref
	cdef bint wildcard_query
	cdef bint debug
	cdef object _dpmatrix
	cdef bytes _reference  # TODO rename to translated_reference or so
	cdef str str_reference

	def __cinit__(self, str reference, double max_error_rate, int flags=SEMIGLOBAL, bint wildcard_ref=False, bint wildcard_query=False):
		self.max_error_rate = max_error_rate
		self.flags = flags
		self.wildcard_ref = wildcard_ref
		self.wildcard_query = wildcard_query
		self.str_reference = reference
		self.reference = reference
		self._min_overlap = 1
		self.debug = False
		self._dpmatrix = None
		self._insertion_cost = 1
		self._deletion_cost = 1

	property min_overlap:
		def __get__(self):
			return self._min_overlap

		def __set__(self, int value):
			if value < 1:
				raise ValueError('Minimum overlap must be at least 1')
			self._min_overlap = value

	property indel_cost:
		"""
		Matches cost 0, mismatches cost 1. Only insertion/deletion costs can be
		changed.
		"""
		def __set__(self, value):
			if value < 1:
				raise ValueError('Insertion/deletion cost must be at leat 1')
			self._insertion_cost = value
			self._deletion_cost = value

	property reference:
		def __get__(self):
			return self._reference

		def __set__(self, str reference):
			mem = <_Entry*> PyMem_Realloc(self.column, (len(reference) + 1) * sizeof(_Entry))
			if not mem:
				raise MemoryError()
			self.column = mem
			self._reference = reference.encode('ascii')
			self.m = len(reference)
			if self.wildcard_ref:
				self._reference = self._reference.translate(IUPAC_TABLE)
			elif self.wildcard_query:
				self._reference = self._reference.translate(ACGT_TABLE)
			self.str_reference = reference

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
		cdef char* s1 = self._reference
		cdef bytes query_bytes = query.encode('ascii')
		cdef char* s2 = query_bytes
		cdef int m = self.m
		cdef int n = len(query)
		cdef _Entry* column = self.column
		cdef double max_error_rate = self.max_error_rate
		cdef bint start_in_ref = self.flags & START_WITHIN_SEQ1
		cdef bint start_in_query = self.flags & START_WITHIN_SEQ2
		cdef bint stop_in_ref = self.flags & STOP_WITHIN_SEQ1
		cdef bint stop_in_query = self.flags & STOP_WITHIN_SEQ2

		if self.wildcard_query:
			query_bytes = query_bytes.translate(IUPAC_TABLE)
			s2 = query_bytes
		elif self.wildcard_ref:
			query_bytes = query_bytes.translate(ACGT_TABLE)
			s2 = query_bytes
		cdef bint compare_ascii = not (self.wildcard_query or self.wildcard_ref)
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
		if not start_in_query:
			# costs can only get worse after column m
			max_n = min(n, m + k)
		if not stop_in_query:
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
		if not start_in_ref and not start_in_query:
			for i in range(m + 1):
				column[i].matches = 0
				column[i].cost = max(i, min_n) * self._insertion_cost
				column[i].origin = 0
		elif start_in_ref and not start_in_query:
			for i in range(m + 1):
				column[i].matches = 0
				column[i].cost = min_n * self._insertion_cost
				column[i].origin = min(0, min_n - i)
		elif not start_in_ref and start_in_query:
			for i in range(m + 1):
				column[i].matches = 0
				column[i].cost = i * self._insertion_cost
				column[i].origin = max(0, min_n - i)
		else:
			for i in range(m + 1):
				column[i].matches = 0
				column[i].cost = min(i, min_n) * self._insertion_cost
				column[i].origin = min_n - i

		if self.debug:
			self._dpmatrix = DPMatrix(self.str_reference, query)
			for i in range(m + 1):
				self._dpmatrix.set_entry(i, min_n, column[i].cost)
		cdef _Match best
		best.ref_stop = m
		best.query_stop = n
		best.cost = m + n
		best.origin = 0
		best.matches = 0

		# Ukkonen's trick: index of the last cell that is less than k.
		cdef int last = min(m, k + 1)
		if start_in_ref:
			last = m

		cdef int cost_diag
		cdef int cost_deletion
		cdef int cost_insertion
		cdef int origin, cost, matches
		cdef int length
		cdef bint characters_equal
		cdef _Entry tmp_entry

		with nogil:
			# iterate over columns
			for j in range(min_n + 1, max_n + 1):
				# remember first entry
				tmp_entry = column[0]

				# fill in first entry in this column
				if start_in_query:
					column[0].origin = j
				else:
					column[0].cost = j * self._insertion_cost
				for i in range(1, last + 1):
					if compare_ascii:
						characters_equal = (s1[i-1] == s2[j-1])
					else:
						characters_equal = (s1[i-1] & s2[j-1]) != 0
					if characters_equal:
						# Characters match: This cannot be an indel.
						cost = tmp_entry.cost
						origin = tmp_entry.origin
						matches = tmp_entry.matches + 1
					else:
						# Characters do not match.
						cost_diag = tmp_entry.cost + 1
						cost_deletion = column[i].cost + self._deletion_cost
						cost_insertion = column[i-1].cost + self._insertion_cost

						if cost_diag <= cost_deletion and cost_diag <= cost_insertion:
							# MISMATCH
							cost = cost_diag
							origin = tmp_entry.origin
							matches = tmp_entry.matches
						elif cost_insertion <= cost_deletion:
							# INSERTION
							cost = cost_insertion
							origin = column[i-1].origin
							matches = column[i-1].matches
						else:
							# DELETION
							cost = cost_deletion
							origin = column[i].origin
							matches = column[i].matches

					# remember current cell for next iteration
					tmp_entry = column[i]

					column[i].cost = cost
					column[i].origin = origin
					column[i].matches = matches
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
					cost = column[m].cost
					matches = column[m].matches
					if length >= self._min_overlap and cost <= length * max_error_rate and (matches > best.matches or (matches == best.matches and cost < best.cost)):
						# update
						best.matches = matches
						best.cost = cost
						best.origin = column[m].origin
						best.ref_stop = m
						best.query_stop = j
						if cost == 0 and matches == m:
							# exact match, stop early
							break
				# column finished

		if max_n == n:
			first_i = 0 if stop_in_ref else m
			# search in last column # TODO last?
			for i in range(first_i, m+1):
				length = i + min(column[i].origin, 0)
				cost = column[i].cost
				matches = column[i].matches
				if length >= self._min_overlap and cost <= length * max_error_rate and (matches > best.matches or (matches == best.matches and cost < best.cost)):
					# update best
					best.matches = matches
					best.cost = cost
					best.origin = column[i].origin
					best.ref_stop = i
					best.query_stop = n
		if best.cost == m + n:
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

		assert best.ref_stop - start1 > 0  # Do not return empty alignments.
		return (start1, best.ref_stop, start2, best.query_stop, best.matches, best.cost)

	def __dealloc__(self):
		PyMem_Free(self.column)


def locate(str reference, str query, double max_error_rate, int flags=SEMIGLOBAL, bint wildcard_ref=False, bint wildcard_query=False, int min_overlap=1):
	aligner = Aligner(reference, max_error_rate, flags, wildcard_ref, wildcard_query)
	aligner.min_overlap = min_overlap
	return aligner.locate(query)


def compare_prefixes(str ref, str query, bint wildcard_ref=False, bint wildcard_query=False):
	"""
	Find out whether one string is the prefix of the other one, allowing
	IUPAC wildcards in ref and/or query if the appropriate flag is set.

	This is used to find an anchored 5' adapter (type 'FRONT') in the 'no indels' mode.
	This is very simple as only the number of errors needs to be counted.

	This function returns a tuple compatible with what Aligner.locate outputs.
	"""
	cdef int m = len(ref)
	cdef int n = len(query)
	cdef bytes query_bytes = query.encode('ascii')
	cdef bytes ref_bytes = ref.encode('ascii')
	cdef char* r_ptr
	cdef char* q_ptr
	cdef int length = min(m, n)
	cdef int i, matches = 0
	cdef bint compare_ascii = False

	if wildcard_ref:
		ref_bytes = ref_bytes.translate(IUPAC_TABLE)
	elif wildcard_query:
		ref_bytes = ref_bytes.translate(ACGT_TABLE)
	else:
		compare_ascii = True
	if wildcard_query:
		query_bytes = query_bytes.translate(IUPAC_TABLE)
	elif wildcard_ref:
		query_bytes = query_bytes.translate(ACGT_TABLE)

	if compare_ascii:
		for i in range(length):
			if ref[i] == query[i]:
				matches += 1
	else:
		r_ptr = ref_bytes
		q_ptr = query_bytes
		for i in range(length):
			if (r_ptr[i] & q_ptr[i]) != 0:
				matches += 1

	# length - matches = no. of errors
	return (0, length, 0, length, matches, length - matches)
