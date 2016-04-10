"""
Classes for managing output files.
"""
import os
from . import seqio

import sys
UNCLOSEABLE = (sys.stdout, sys.stderr)


class SingleEndWriter(object):
	"""
	Wrapper around a single output file. Can safely wrap
	sys.stdout and sys.stderr.
	"""
	def __init__(self, handle):
		self.handle = handle
		self.written = 0
		self.read1_bp = 0
		self.read2_bp = 0
	
	def write(self, read1, read2=None):
		assert read2 is None
		self.handle.write(read1)
		self.written += 1
		self.read1_bp += len(read1)
	
	@property
	def written_bp(self):
		"""
		Returns the total number of read1 and read2 bases
		written to the output file.
		"""
		return (self.read1_bp, self.read2_bp)
	
	def close(self):
		if self.handle not in UNCLOSEABLE:
			self.handle.close()


class PairedEndWriter(SingleEndWriter):
	"""
	Wrapper around a pair of output files.
	"""
	def write(self, read1, read2):
		self.handle.write(read1, read2)
		self.written += 1
		self.read1_bp += len(read1)
		self.read2_bp += len(read2)


class Writers(object):
	"""
	Manage output files. Each output file (single-end) or pair of files
	(paired-end) is associated with a filter type. Multiplexing is also
	supported, in which case a filter type of `None` indicates that the
	name of the output file to be used is determined by replacing '{name}'
	in `name_pattern` with the adapter name. Files are opened lazily so 
	that a `Writers` object can be constructed in the main thread and 
	then be safely passed to a separate writer thread.
	"""
	def __init__(self, multiplexed=False, name_pattern=None, **seqio_open_args):
		"""
		Create a Writers object.
		
		multiplexed -- whether the default writer should be muliplexed,
		  i.e. a separate file will be opened for each adapter name. If True, 
		  `name_pattern` must contain '{name}', which will be replaced with the
		  adapter name.
		name_pattern -- The file name pattern for multiplexed files; ignored 
		  unless `multiplexed` = True
		seqio_open_args -- arguments to pass to seqio.open
		"""
		if multiplexed:
			assert '{name}' in name_pattern
		else:
			output = None
		self.multiplexed = multiplexed
		self.name_pattern = name_pattern
		self.seqio_open_args = seqio_open_args
		self._writers = {}
		self._seqfile_paths = {}
		self._force_create = {}
		self.discarded = 0
	
	def add_writer(self, filter_type, file1, file2=None, force_create=False):
		"""
		Add an output file (or pair of output files) for a specific
		filter type.
		
		filter_type -- class of filter; reads that trigger this filter type
		  will be written to the specified output file(s).
		file1 -- read1 output file
		file2 -- read2 output file (if any)
		force_create -- whether the output file(s) should always be created,
		  even if they will be empty.
		"""
		self._seqfile_paths[filter_type] = (file1, file2)
		if force_create and isinstance(file1, str) and file1 != "-":
			self._force_create[file1] = False
			if file2 is not None:
				self._force_create[file2] = False
	
	def has_writer(self, filter_type):
		return filter_type in self._seqfile_paths
	
	def get_writer(self, filter_type):
		paths = self._seqfile_paths[filter_type]
		if paths not in self._writers:
			self._writers[paths] = self._create_writer(*paths)
		return self._writers[paths]
	
	def get_multiplex_writer(self, name):
		"""
		If the program is being run in multiplexed mode, get the
		output file for the specific name.
		"""
		assert self.multiplexed
		path = self.name_pattern.format(name=name)
		if path not in self._writers:
			self._writers[path] = self._create_writer(path)
		return self._writers[path]
	
	def _create_writer(self, file1, file2=None):
		seqfile = seqio.open(file1, file2, mode='w', **self.seqio_open_args)
		if file1 in self._force_create:
			self._force_create[file1] = True
		if file2 is not None and file2 in self._force_create:
			self._force_create[file2] = True
		if isinstance(seqfile, seqio.SingleRecordWriter):
			return SingleEndWriter(seqfile)
		elif isinstance(seqfile, seqio.PairRecordWriter):
			return PairedEndWriter(seqfile)
		else:
			raise Exception("Unrecognized type of writer {}".format(writer.__class__))
	
	@property
	def writers(self):
		"""
		Returns all the open output files.
		"""
		return self._writers.values()
	
	def summary(self):
		"""
		Returns a tuple (total number of records written, 
		(total number of read1 bp written, total number of read2 bp written)).
		"""
		writers = self.writers
		written = sum(w.written for w in writers)
		written_bp = (
			sum(w.read1_bp for w in writers),
			sum(w.read2_bp for w in writers),
		)
		return (written, written_bp)
	
	def write(self, filter_type, read1, read2=None):
		"""
		Write read(s) to the correct output file(s) (if any)
		for the given filter type.
		filter_type -- class of filter that was triggered by read(s).
		"""
		writer = None
		
		if filter_type is None and self.multiplexed and read1.match:
			name = read1.match.adapter.name
			writer = self.get_multiplex_writer(name)
		
		elif self.has_writer(filter_type):
			writer = self.get_writer(filter_type)
		
		if writer is not None:
			writer.write(read1, read2)
		else:
			self.discarded += 1
	
	def close(self):
		# touch any files in force_create that haven't been created
		for path, created in self._force_create.items():
			if not created:
				with open(path, 'w'): os.utime(path, None)
		# close any open writers
		for writer in self._writers.values():
			if writer is not None:
				writer.close()
