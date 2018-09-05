import multiprocessing
import re


def available_cpu_count():
	"""
	Return the number of available virtual or physical CPUs on this system.
	The number of available CPUs can be smaller than the total number of CPUs
	when the cpuset(7) mechanism is in use, as is the case on some cluster
	systems.

	Adapted from http://stackoverflow.com/a/1006301/715090
	"""
	try:
		with open('/proc/self/status') as f:
			status = f.read()
		m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', status)
		if m:
			res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
			if res > 0:
				return min(res, multiprocessing.cpu_count())
	except IOError:
		pass

	return multiprocessing.cpu_count()
