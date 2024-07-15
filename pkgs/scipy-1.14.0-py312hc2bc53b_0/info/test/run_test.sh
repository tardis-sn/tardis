

set -ex



test ! -d $SP_DIR/scipy/_lib/tests
test ! -d $SP_DIR/scipy/cluster/tests
test ! -d $SP_DIR/scipy/constants/tests
test ! -d $SP_DIR/scipy/datasets/tests
test ! -d $SP_DIR/scipy/fftpack/tests
test ! -d $SP_DIR/scipy/fft/tests
test ! -d $SP_DIR/scipy/integrate/tests
test ! -d $SP_DIR/scipy/interpolate/tests
test ! -d $SP_DIR/scipy/io/tests
test ! -d $SP_DIR/scipy/linalg/tests
test ! -d $SP_DIR/scipy/misc/tests
test ! -d $SP_DIR/scipy/ndimage/tests
test ! -d $SP_DIR/scipy/odr/tests
test ! -d $SP_DIR/scipy/optimize/tests
test ! -d $SP_DIR/scipy/signal/tests
test ! -d $SP_DIR/scipy/sparse/tests
test ! -d $SP_DIR/scipy/spatial/tests
test ! -d $SP_DIR/scipy/special/tests
test ! -d $SP_DIR/scipy/stats/tests
test ! -f $SP_DIR/scipy/_lib/_test_ccallback.cpython-312-x86_64-linux-gnu.so
test ! -f $SP_DIR/scipy/integrate/_test_multivariate.cpython-312-x86_64-linux-gnu.so
test ! -f $SP_DIR/scipy/io/_test_fortran.cpython-312-x86_64-linux-gnu.so
test ! -f $SP_DIR/scipy/ndimage/_ctest.cpython-312-x86_64-linux-gnu.so
test ! -f $SP_DIR/scipy/ndimage/_cytest.cpython-312-x86_64-linux-gnu.so
test ! -f $SP_DIR/scipy/special/_test_internal.cpython-312-x86_64-linux-gnu.so
python -c "import scipy; throw_away_the_return_value = scipy.test()" > testlog
python -c "import sys; lines=open('testlog').readlines(); sys.exit(0 if any('conda-forge builds of' in x for x in lines) else 1)"
(pytest --pyargs scipy || echo "failure was expected") > testlog
python -c "import sys; lines=open('testlog').readlines(); sys.exit(0 if any('conda-forge builds of' in x for x in lines) else 1)"
python -c "import sys; lines=open('testlog').readlines(); sys.exit(0 if any('======== 1 failed' in x for x in lines) else 1)"
exit 0
