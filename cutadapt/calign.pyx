# kate: syntax Python;
from cpython.mem cimport PyMem_Malloc, PyMem_Free

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
DEF WILDCARD_CHAR = 'N'

# structure for a DP matrix entry
ctypedef struct Entry:
	int cost
	int matches  # no. of matches in this alignment
	int origin   # where the alignment originated: negative for positions within seq1, positive for pos. within seq2


def globalalign_locate(char* s1, char* s2, double max_error_rate, int flags = SEMIGLOBAL, int degenerate = 0):
	"""
	globalalign_locate(string1, string2, max_error_rate, flags=SEMIGLOBAL)
		-> (start1, stop1, start2, stop2, matches, errors)

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
	cdef int m = len(s1)
	cdef int n = len(s2)

	"""
	DP Matrix:
	           s2 (j)
	        ----------> n
	       |
	s1 (i) |
	       |
	       V
	       m
	"""
	# only a single column of the DP matrix is stored
	cdef Entry* column = <Entry*> PyMem_Malloc((m+1) * sizeof(Entry))
	#TODO if column == NULL:
	#	return NULL;
	cdef int i, j, best_i, best_j, best_cost, best_matches, best_origin

	# initialize first column
	for i in range(0, m+1):
		column[i].matches = 0
		column[i].cost = 0 if (flags & START_WITHIN_SEQ1) else (i * DELETION_COST)
		column[i].origin = -i if (flags & START_WITHIN_SEQ1) else 0

	best_i = m
	best_j = 0
	best_cost = column[m].cost
	best_matches = 0
	best_origin = column[m].origin

	# maximum no. of errors
	cdef int k = <int> (max_error_rate * m)
	cdef int last = k + 1

	cdef int match
	cdef int cost_diag
	cdef int cost_deletion
	cdef int cost_insertion
	cdef int origin, cost, matches
	cdef int length
	cdef int wildcard1 = degenerate & ALLOW_WILDCARD_SEQ1
	cdef int wildcard2 = degenerate & ALLOW_WILDCARD_SEQ2
	cdef int tmp, tmp2
	cdef Entry tmp_entry

	if flags & START_WITHIN_SEQ1:
		last = m

	# determine largest column we need to compute
	cdef int max_n = n
	if not (flags & START_WITHIN_SEQ2):
		# costs can only get worse after column m
		max_n = min(max_n, m+k)
	# iterate over columns
	for j in range(1, max_n+1):
		# remember first entry
		tmp_entry = column[0]

		# fill in first entry in this column TODO move out of loop
		if flags & START_WITHIN_SEQ2:
			column[0].cost = 0
			column[0].origin = j
			column[0].matches = 0
		else:
			column[0].cost = j * INSERTION_COST
			column[0].origin = 0
			column[0].matches = 0
		for i in range(1, last + 1):
			# the following code is ugly due to a bug in Cython up to at least 0.17.2
			match = (s1[i-1] == s2[j-1])
			tmp = wildcard1
			tmp2 = (s1[i-1] == WILDCARD_CHAR)
			tmp = tmp and tmp2
			match = match or tmp
			tmp = wildcard2
			tmp2 = (s2[j-1] == WILDCARD_CHAR)
			tmp = tmp and tmp2
			match = match or tmp

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

		while column[last].cost > k:
			last -= 1
		if last < m:
			last += 1
		else:
			# found
			# if requested, find best match in last row
			if flags & STOP_WITHIN_SEQ2:
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

	if max_n == n and flags & STOP_WITHIN_SEQ1:
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

	PyMem_Free(column)
	cdef int start1, start2
	if best_origin >= 0:
		start1 = 0
		start2 = best_origin
	else:
		start1 = -best_origin
		start2 = 0

	return (start1, best_i, start2, best_j, best_matches, best_cost)
