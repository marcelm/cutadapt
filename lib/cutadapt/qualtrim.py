"""
Quality trimming.
"""

def quality_trim_index(qualities, cutoff):
	"""
	Find the position at which to trim a low-quality end from a nucleotide sequence.

	Qualities are assumed to be ASCII-encoded as chr(qual + 33).

	>>> trim("", "", 10) TODO
	with qualities based on bwa_trim_read

	The algorithm is the same as the one used by BWA within the function
	'bwa_trim_read':
	- Subtract the cutoff value from all qualities.
	- Compute partial sums from all indices to the end of the sequence.
	- Trim sequence at the index at which the sum is minimal.
	"""
	s = 0
	max_qual = 0
	max_i = len(qualities)
	for i in xrange(len(qualities)-1, -1, -1):
		q = ord(qualities[i]) - 33
		s += cutoff - q
		if s < 0:
			break
		if s > max_qual:
			max_qual = s
			max_i = i
	return max_i
