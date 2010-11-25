# kate: syntax python;

from array import array

# direction constants for backtrace table
cdef enum DIRECTION:
	LEFT = 1
	TOP = 2
	DIAG = 3

ctypedef signed short score_t

def semiglobalalign(char* s1, char* s2, print_table=False):
	"""
Computes an end-gap free alignment (also called free-shift alignment) of
strings s1 and s2, using unit costs.

Return a tuple (r1, r2, start1, stop1, start2, stop2, errors)
where r1 and r2 are sequences of the same length containing the alignment
(an INDEL is marked by '-').

start1 is the position within r1 at which the part of s1, that is aligned, starts.
start2 is the position within r1 at which the part of s1, that is aligned, ends.
The same holds for start2, stop2.

It is always the case that at least one of start1 and start2 is zero.

It is always the case that either stop1==len(r1) or stop2==len(r2) or both
(note that len(r1)==len(r2)). This is a property of semiglobal alignments.

errors is the number of errors in the alignment.

For example, pysemiglobalalign("SISSI", "MISSISSIPPI") returns:

r1 = [ '-', '-', '-', 'S', 'I', 'S', 'S', 'I', '-', '-', '-']
r2 = [ 'M', 'I', 'S', 'S', 'I', 'S', 'S', 'I', 'P', 'P', 'I']
start1, stop1 = 0, 5
start2, stop2 = 3, 8
errors = 0

This corresponds to the following alignment (formatted with print_alignment):
   SISSI
   |||||
MISSISSIPPI
"""

	#             s2 (n, j)
	#      --------------->
	#     |
	#  s1 | (m, i)
	#     |
	#     V
	cdef int m = len(s1)
	cdef int n = len(s2)
	
	cdef score_t score, tmp
	cdef DIRECTION bt

	# the DP and backtrace table are both stored column-wise
	# It's much faster to use two tables instead
	# of one with tuples. (Tried both.)
	columns = [ array('h', (m+1)*[0]) for x in xrange(n+1) ]
	backtrace = [ array('h', (m+1)*[0]) for x in xrange(n+1) ]

	# calculate alignment (using unit costs)
	# outer loop goes over columns
	for j in range(1, n+1):
		cur_column = columns[j]
		prev_column = columns[j-1]
		for i in range(1, m+1):
			bt = DIAG
			score = prev_column[i-1] + (1 if (s1[i-1] == s2[j-1]) else -1)
			tmp = cur_column[i-1] - 1
			if tmp > score:
				bt = LEFT
				score = tmp
			tmp = prev_column[i] - 1
			if tmp > score:
				bt = TOP
				score = tmp
			cur_column[i] = score
			backtrace[j][i] = bt

	# find position with highest score in last column or last row
	last_row = [ (columns[j][m], j) for j in xrange(n+1) ]
	row_best, best_j = max(last_row)

	last_column = [ (columns[n][i], i) for i in xrange(m+1) ]
	col_best, best_i = max(last_column)

	if row_best > col_best:
		best_i = m
		best = row_best
	else:
		best_j = n
		best = col_best

	# now trace back (from "lower right" corner)
	r1 = []
	r2 = []
	i = m
	j = n

	# first, walk from the lower right corner to the
	# position where we found the maximum score
	if i == best_i: # we are in the last row
		while j > best_j:
			r1.append('-')
			r2.append(s2[j-1])
			j -= 1
	else: # we are in the last column
		while i > best_i:
			r1.append(s1[i-1])
			r2.append('-')
			i -= 1

	rlen = len(r1)

	errors = 0
	# now, the actual backtracing
	# we build reverse sequences while backtracing and
	# reverse them afterwards.
	while i > 0 and j > 0: #columns[j][i] > 0:
		score = columns[j][i]
		direction = backtrace[j][i]
		if direction == DIAG:
			if s1[i-1] != s2[j-1]:
				errors += 1
			r1.append(s1[i-1])
			r2.append(s2[j-1])
			i -= 1
			j -= 1
		elif direction == TOP:
			r1.append('-')
			r2.append(s2[j-1])
			errors += 1
			j -= 1
		elif direction == LEFT:
			r1.append(s1[i-1])
			r2.append('-')
			errors += 1
			i -= 1

	# compute the length of the actual alignment (ignoring ends)
	length = len(r1) - rlen

	#begin = max(i, j)
	start1 = i
	start2 = j

	# continue from where we are to the upper left corner
	#while i > 0 and j > 0:
		#r1.append(s1[i-1])
		#r2.append(s2[j-1])
		#i -= 1
		#j -= 1
	while j > 0:
		r1.append('-')
		r2.append(s2[j-1])
		j -= 1
	while i > 0:
		r1.append(s1[i-1])
		r2.append('-')
		i -= 1
	assert i == 0 and j == 0

	#print "best score:", columns[best_j][best_i]
	#print "length:", length
	#print "errors:", errors
	assert columns[best_j][best_i] == length - 2*errors
	# reverse result
	r1.reverse()
	r2.reverse()

	stop1 = best_i
	stop2 = best_j
	return (r1, r2, start1, stop1, start2, stop2, errors)

