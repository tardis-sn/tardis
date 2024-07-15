#!/bin/bash

set -e

export NUMBA_DEVELOPER_MODE=1
export NUMBA_DISABLE_ERROR_MESSAGE_HIGHLIGHTING=1
export PYTHONFAULTHANDLER=1

unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
  SEGVCATCH=catchsegv
  export CC="${CC} -pthread"
elif [[ "$unamestr" == 'Darwin' ]]; then
  SEGVCATCH=""
else
  echo Error
fi

# limit CPUs in use on PPC64LE and AARCH64, fork() issues
# occur on high core count systems
archstr=`uname -m`
if [[ "$archstr" == 'ppc64le' ]]; then
    TEST_NPROCS=1
elif [[ "$archstr" == 'aarch64' ]]; then
    TEST_NPROCS=4
else
    TEST_NPROCS=${CPU_COUNT}
fi

# Disable NumPy dispatching to AVX512_SKX feature extensions if the chip is
# reported to support the feature and NumPy >= 1.22 as this results in the use
# of low accuracy SVML libm replacements in ufunc loops.
_NPY_CMD='from numba.misc import numba_sysinfo;\
          sysinfo=numba_sysinfo.get_sysinfo();\
          print(sysinfo["NumPy AVX512_SKX detected"] and
                sysinfo["NumPy Version"]>="1.22")'
NUMPY_DETECTS_AVX512_SKX_NP_GT_122=$(python -c "$_NPY_CMD")
echo "NumPy >= 1.22 with AVX512_SKX detected: $NUMPY_DETECTS_AVX512_SKX_NP_GT_122"

if [[ "$NUMPY_DETECTS_AVX512_SKX_NP_GT_122" == "True" ]]; then
    export NPY_DISABLE_CPU_FEATURES="AVX512_SKX"
fi

# Validate Numba dependencies
python -m pip check

# Check Numba executables are there
numba -h

# run system info tool
numba -s

# Check test discovery works
python -m numba.tests.test_runtests

if [[ "$archstr" == 'aarch64' ]] || [[ "$archstr" == "ppc64le" ]]; then
	echo "Skipping numba test suite on $archstr"
	#echo 'Running only a random selection of tests'
	#$SEGVCATCH python -m numba.runtests -b --random='0.15' --exclude-tags='long_running' -m $TEST_NPROCS -- numba.tests
# Else run the whole test suite
else
	echo 'Running all the tests except long_running'
	echo "Running: $SEGVCATCH python -m numba.runtests -b -m $TEST_NPROCS -- $TESTS_TO_RUN"
$SEGVCATCH python -m numba.runtests -b --exclude-tags='long_running' -m $TEST_NPROCS -- $TESTS_TO_RUN
fi
