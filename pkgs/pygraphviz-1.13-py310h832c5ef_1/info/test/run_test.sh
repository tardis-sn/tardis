

set -ex



pytest --doctest-modules --durations=10 --pyargs pygraphviz
exit 0
