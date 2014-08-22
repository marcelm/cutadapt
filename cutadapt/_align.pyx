from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc
DEF START_WITHIN_SEQ1 = 1
DEF START_WITHIN_SEQ2 = 2
DEF STOP_WITHIN_SEQ1 = 4
DEF STOP_WITHIN_SEQ2 = 8
DEF SEMIGLOBAL = 15
DEF ALLOW_WILDCARD_SEQ1 = 1
DEF ALLOW_WILDCARD_SEQ2 = 2

DEF INSERTION_COST = 1
DEF DELETION_COST = 1
DEF MATCH_COST = 0
DEF MISMATCH_COST = 1
DEF WILDCARD_CHAR = b'N'

# structure for a DP matrix entry
ctypedef struct Entry:
	int cost
	int matches  # no. of matches in this alignment
	int origin   # where the alignment originated: negative for positions within seq1, positive for pos. within seq2


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

	The alignment itself is not returned, only the tuple
	(start1, stop1, start2, stop2, matches, errors), where the first four fields have the
	meaning as described, matches is the number of matches and errors is the number of
	errors in the alignment.

	It is always the case that at least one of start1 and start2 is zero.
	"""
	cdef int m
	cdef Entry* column  # one column of the DP matrix
	cdef double max_error_rate
	cdef int flags
	cdef int degenerate
	cdef bytes reference

	def __cinit__(self, char* reference, double max_error_rate, int flags = SEMIGLOBAL, int degenerate = 0):
		self.m = len(reference)
		self.column = <Entry*> PyMem_Malloc((self.m + 1) * sizeof(Entry))
		if not self.column:
			raise MemoryError()
		self.reference = reference
		self.degenerate = degenerate
		self.flags = flags
		self.max_error_rate = max_error_rate

	property reference:
		def __get__(self):
			return self.reference

		def __set__(self, char* reference):
			mem = <Entry*> PyMem_Realloc(self.column, (len(reference) + 1) * sizeof(Entry))
			if not mem:
				raise MemoryError()
			self.column = mem
			self.reference = reference
			self.m = len(reference)

	def locate(self, char* query):
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
		cdef char* s1 = self.reference
		cdef char* s2 = query
		cdef char* ref = self.reference  # s1.   s2 = query.
		cdef int n = len(query)
		cdef int m = self.m
		cdef Entry* column = self.column
		cdef double max_error_rate = self.max_error_rate
		cdef bint start_in_ref = self.flags & START_WITHIN_SEQ1
		cdef bint start_in_query = self.flags & START_WITHIN_SEQ2
		cdef bint stop_in_ref = self.flags & STOP_WITHIN_SEQ1
		cdef bint stop_in_query = self.flags & STOP_WITHIN_SEQ2
		cdef bint wildcard1 = self.degenerate & ALLOW_WILDCARD_SEQ1
		cdef bint wildcard2 = self.degenerate & ALLOW_WILDCARD_SEQ2
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
			for i in range(0, m + 1):
				column[i].matches = 0
				column[i].cost = max(i, min_n)
				column[i].origin = 0
		elif start_in_ref and not start_in_query:
			for i in range(0, m + 1):
				column[i].matches = 0
				column[i].cost = min_n
				column[i].origin = min(0, min_n - i)
		elif not start_in_ref and start_in_query:
			for i in range(0, m + 1):
				column[i].matches = 0
				column[i].cost = i
				column[i].origin = max(0, min_n - i)
		else:
			for i in range(0, m + 1):
				column[i].matches = 0
				column[i].cost = min(i, min_n)
				column[i].origin = min_n - i

		cdef int best_i = m
		cdef int best_j = 0
		cdef int best_cost = column[m].cost
		cdef int best_matches = 0
		cdef int best_origin = column[m].origin

		# Ukkonen's trick: index of the last cell that is less than k.
		cdef int last = k + 1
		if start_in_ref:
			last = m

		cdef int match
		cdef int cost_diag
		cdef int cost_deletion
		cdef int cost_insertion
		cdef int origin, cost, matches
		cdef int length
		cdef Entry tmp_entry

		# iterate over columns
		for j in range(min_n + 1, max_n + 1):
			# remember first entry
			tmp_entry = column[0]

			# fill in first entry in this column TODO move out of loop
			assert column[0].matches == 0
			if start_in_query:
				assert column[0].cost == 0
				assert column[0].origin == j-1
				column[0].cost = 0
				column[0].origin = j
				column[0].matches = 0
			else:
				assert column[0].cost == j-1
				assert column[0].origin == 0
				column[0].cost = j * INSERTION_COST
				column[0].origin = 0
				column[0].matches = 0
			for i in range(1, last + 1):
				if s1[i-1] == s2[j-1] or \
						(wildcard1 and s1[i-1] == WILDCARD_CHAR) or \
						(wildcard2 and s2[j-1] == WILDCARD_CHAR):
					match = 1
				else:
					match = 0

				cost_diag = tmp_entry.cost + (MATCH_COST if match else MISMATCH_COST)
				cost_deletion = column[i].cost + DELETION_COST
				cost_insertion = column[i-1].cost + INSERTION_COST

				if cost_diag <= cost_deletion and cost_diag <= cost_insertion:
					# MATCH or MISMATCH
					cost = cost_diag
					origin = tmp_entry.origin
					matches = tmp_entry.matches + match
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
			while last >= 0 and column[last].cost > k:
				last -= 1
			# last can be -1 here, but will be incremented next.
			# TODO if last is -1, can we stop searching?
			if last < m:
				last += 1
			else:
				# Found. If requested, find best match in last row
				if stop_in_query:
					# length of the aligned part of string1
					length = m + min(column[m].origin, 0)
					cost = column[m].cost
					matches = column[m].matches
					if cost <= length * max_error_rate and (matches > best_matches or (matches == best_matches and cost < best_cost)):
						# update
						best_matches = matches
						best_cost = cost
						best_origin = column[m].origin
						best_i = m
						best_j = j
			# column finished

		if max_n == n and stop_in_ref:
			# search in last column # TODO last?
			for i in range(0, m+1):
				length = i + min(column[i].origin, 0)
				cost = column[i].cost
				matches = column[i].matches
				if cost <= length * max_error_rate and (matches > best_matches or (matches == best_matches and cost < best_cost)):
					# update best
					best_matches = matches
					best_cost = cost
					best_origin = column[i].origin
					best_i = i
					best_j = n

		cdef int start1, start2
		if best_origin >= 0:
			start1 = 0
			start2 = best_origin
		else:
			start1 = -best_origin
			start2 = 0

		return (start1, best_i, start2, best_j, best_matches, best_cost)

	def __dealloc__(self):
		PyMem_Free(self.column)


def locate(char* reference, char* query, double max_error_rate, int flags = SEMIGLOBAL, int degenerate = 0):
	aligner = Aligner(reference, max_error_rate, flags, degenerate)
	return aligner.locate(query)


def compare_prefixes(char* s1, char* s2, int degenerate = 0):
	"""
	Find out whether two strings have a common prefix, allowing mismatches.

	This is used to find an anchored 5' adapter (type 'FRONT') in the 'no indels' mode.
	This is very simple as only the number of errors needs to be counted.

	This function returns a tuple compatible with what globalalign_locate outputs.
	"""
	cdef int m = len(s1)
	cdef int n = len(s2)
	cdef int wildcard1 = degenerate & ALLOW_WILDCARD_SEQ1
	cdef int wildcard2 = degenerate & ALLOW_WILDCARD_SEQ2
	cdef int length = min(m, n)
	cdef int i, matches = 0
	for i in range(length):
		if s1[i] == s2[i] or \
				(wildcard1 and s1[i] == WILDCARD_CHAR) or \
				(wildcard2 and s2[i] == WILDCARD_CHAR):
			matches += 1
	# length - matches = no. of errors
	return (0, length, 0, length, matches, length - matches)
