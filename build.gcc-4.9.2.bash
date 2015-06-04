#!/usr/bin/env bash
set -e -o pipefail

# script to build agg on CHUK cluster using gcc 4.9.2

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ] ; do SOURCE="$(readlink "$SOURCE")"; done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

export PATH=/illumina/thirdparty/gcc/gcc-4.9.2/bin/:/bin:/usr/bin:/usr/X11R6/bin:/usr/local/bin:.
export CXX=/illumina/thirdparty/gcc/gcc-4.9.2/bin/g++
export CC=/illumina/thirdparty/gcc/gcc-4.9.2/bin/gcc
export LD_LIBRARY_PATH=/illumina/thirdparty/gcc/${EL}/gcc-4.9.2/lib

# we want to use our own version of boost, so unset all user-defined boost environment variables
unset BOOST_ROOT
unset BOOST_DIR
unset BOOST_ROOT
unset BOOST_INCLUDEDIR
unset BOOST_LIBRARYDIR

make CXX=${CXX} CC=${CC} 

exit $?

