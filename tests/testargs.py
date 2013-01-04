from __future__ import print_function, division, absolute_import

from cutadapt.scripts import cutadapt

from nose.tools import raises


@raises(SystemExit)
def test_no_args():
	cutadapt.main()
