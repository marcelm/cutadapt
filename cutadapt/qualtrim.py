"""
Quality trimming.
"""
import sys
if sys.version_info[0] >= 3:
	xrange = range


def quality_trim_index(qualities, cutoff, base=33):
	"""
	Find the position at which to trim a low-quality end from a nucleotide sequence.

	Qualities are assumed to be ASCII-encoded as chr(qual + base).

	The algorithm is the same as the one used by BWA within the function
	'bwa_trim_read':
	- Subtract the cutoff value from all qualities.
	- Compute partial sums from all indices to the end of the sequence.
	- Trim sequence at the index at which the sum is minimal.
	"""
	s = 0
	max_qual = 0
	max_i = len(qualities)
	for i in reversed(xrange(len(qualities))):
		q = ord(qualities[i:i+1]) - base
		s += cutoff - q
		if s < 0:
			break
		if s > max_qual:
			max_qual = s
			max_i = i
	return max_i
