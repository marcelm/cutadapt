#!/bin/bash
#
# Build manylinux wheels. Based on the example at
# <https://github.com/pypa/python-manylinux-demo>
#
# It is best to run this in a fresh clone of the repository!
#
# Run this within the repository root:
#   ./buildwheels.sh
#
# The wheels will be put into the wheelhouse/ subdirectory.
#
# For interactive tests:
#   docker run -it -v $(pwd):/io quay.io/pypa/manylinux2010_x86_64 /bin/bash

set -xeuo pipefail

manylinux=quay.io/pypa/manylinux2010_x86_64

# For convenience, if this script is called from outside of a docker container,
# it starts a container and runs itself inside of it.
if ! grep -q docker /proc/1/cgroup; then
  # We are not inside a container
  docker pull ${manylinux}
  exec docker run --rm -v $(pwd):/io ${manylinux} /io/$0
fi

# Strip binaries (copied from multibuild)
STRIP_FLAGS=${STRIP_FLAGS:-"-Wl,-strip-all"}
export CFLAGS="${CFLAGS:-$STRIP_FLAGS}"
export CXXFLAGS="${CXXFLAGS:-$STRIP_FLAGS}"

# We require Python 3.5+
rm /opt/python/cp27* /opt/python/cp34*

PYBINS="/opt/python/*/bin"
HAS_CYTHON=0
for PYBIN in ${PYBINS}; do
#    ${PYBIN}/pip install -r /io/requirements.txt
    ${PYBIN}/pip wheel --no-deps /io/ -w wheelhouse/
done
ls wheelhouse/

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" --plat manylinux1_x86_64 -w repaired/
done

# Created files are owned by root, so fix permissions.
chown -R --reference=/io/setup.py repaired/
mv repaired/*.whl /io/dist/

# TODO Install packages and test them
#for PYBIN in ${PYBINS}; do
#    ${PYBIN}/pip install cutadapt --no-index -f /io/wheelhouse
#    (cd $HOME; ${PYBIN}/nosetests ...)
#done
