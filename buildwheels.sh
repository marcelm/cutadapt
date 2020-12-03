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
# The wheels will be put into the dist/ subdirectory.

set -xeuo pipefail

manylinux=quay.io/pypa/manylinux2010_x86_64

# For convenience, if this script is called from outside of a docker container,
# it starts a container and runs itself inside of it.
if ! grep -q docker /proc/1/cgroup && ! test -d /opt/python; then
  # We are not inside a container
  docker pull ${manylinux}
  exec docker run --rm -v $(pwd):/io ${manylinux} /io/$0
fi

if ! test -d /io/dist; then
  mkdir /io/dist
  chown --reference=/io/setup.py /io/dist
fi

# Strip binaries (copied from multibuild)
STRIP_FLAGS=${STRIP_FLAGS:-"-Wl,-strip-all"}
export CFLAGS="${CFLAGS:-$STRIP_FLAGS}"
export CXXFLAGS="${CXXFLAGS:-$STRIP_FLAGS}"

for PYBIN in /opt/python/cp3[6789]-*/bin; do
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
