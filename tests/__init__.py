"""We'd like to unit-test the code after merely running 'setup.py build' which
   requires finding the Cython-compiled extension modules in their temporary
   location.
   This enables that to happen, and is an acceptable kludge (I think) for the
   purposes of unit testing.
"""
from __future__ import print_function, division, absolute_import

import sys, os
from glob import glob
from warnings import warn

# Add build directories (I'd only expect one!) to the front of the path.
build_dirs = glob('build/lib.*/cutadapt')
if build_dirs:
    warn("Testing code found in " + str(build_dirs))
    sys.path[:0] = map(os.path.dirname, build_dirs)
    if 'cutadapt' in sys.modules:
        reload(sys.modules['cutadapt'])
