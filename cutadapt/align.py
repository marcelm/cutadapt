"""
Alignment module.
"""
from __future__ import print_function
from array import array
import sys

# flags for global alignment

# The interpretation of the first flag is:
# An initial portion of seq1 may be skipped at no cost.
# This is equivalent to saying that in the alignment,
# gaps in the beginning of seq2 are free.
#
# The other flags have an equivalent meaning.
START_WITHIN_SEQ1 = 1
START_WITHIN_SEQ2 = 2
STOP_WITHIN_SEQ1 = 4
STOP_WITHIN_SEQ2 = 8
ALLOW_WILDCARD_SEQ1 = 1
ALLOW_WILDCARD_SEQ2 = 2

GLOBAL = 0

# Use this to get regular semiglobal alignment
# (all gaps in the beginning or end are free)
SEMIGLOBAL = START_WITHIN_SEQ1 | START_WITHIN_SEQ2 | STOP_WITHIN_SEQ1 | STOP_WITHIN_SEQ2

SCORE_MATCH = 1
SCORE_MISMATCH = -1
SCORE_DELETION = -2
SCORE_INSERTION = SCORE_DELETION


def _ansired(s):
	return "\x1b[1;31m" + s + "\x1b[00m"


# Pure Python implementation, fallback for when the C module is not available.
# Also useful for testing.
def pysemiglobalalign(s1, s2, print_table=False):
	"""
	Compute a semiglobal alignment of strings s1 and s2.

	Return a tuple (r1, r2, start1, stop1, start2, stop2, errors)
	where r1 and r2 are sequences of the same length containing the alignment
	(an INDEL is marked by '-').

	start1 is the position within r1 at which the part of s1, that is aligned, starts.
	stop2 is the position within r1 at which the part of s1, that is aligned, ends.
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
	m = len(s1)
	n = len(s2)

	# the DP and backtrace table are both stored column-wise
	# It's much faster to use two tables instead
	# of one with tuples. (Tried both.)
	columns = [ array('h', (m+1)*[0]) for x in xrange(n+1) ]
	backtrace = [ array('h', (m+1)*[0]) for x in xrange(n+1) ]

	# direction constants for backtrace table
	LEFT = 1
	TOP = 2
	DIAG = 3

	# fill DP table (using unit costs)
	# outer loop goes over columns
	for j in range(1, n+1):
		cur_column = columns[j]
		prev_column = columns[j-1]
		for i in range(1, m+1):
			bt = DIAG
			score = prev_column[i-1] + (SCORE_MATCH if (s1[i-1] == s2[j-1]) else SCORE_MISMATCH)
			tmp = cur_column[i-1] + SCORE_INSERTION
			if tmp > score:
				bt = LEFT
				score = tmp
			tmp = prev_column[i] + SCORE_DELETION
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

	if print_table:
		print("s1:", s1)
		print("s2:", s2)
		print("best i,j", best_i, best_j)
		print("      ", end="")
		for j in xrange(n):
			print("%3c" % s2[j], end="")
		print()
		for i in xrange(m+1):
			if i > 0:
				print("%3c" % s1[i-1], end="")
			else:
				print("   ", end="")
			for j in xrange(n+1):
				if (i,j) == (best_i, best_j):
					st = _ansired("%3d")
				else:
					st = "%3d"
				sys.stdout.write(st % columns[j][i])
			print()
		print()

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
	# now, the actual backtracing.
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

	start1 = i
	start2 = j

	# continue from where we are to the upper left corner
	while j > 0:
		r1.append('-')
		r2.append(s2[j-1])
		j -= 1
	while i > 0:
		r1.append(s1[i-1])
		r2.append('-')
		i -= 1
	assert i == 0 and j == 0

	#TODO assert columns[best_j][best_i] == length - 2*errors
	r1.reverse()
	r2.reverse()

	stop1 = best_i
	stop2 = best_j
	return (r1, r2, start1, stop1, start2, stop2, errors)


from .calign import globalalign_locate
