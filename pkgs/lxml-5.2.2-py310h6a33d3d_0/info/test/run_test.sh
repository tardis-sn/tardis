

set -ex



pip check
test_file=$(jq '.files[] | select( . | endswith("/test.py"))' $CONDA_PREFIX/conda-meta/lxml-5.2.2-${PKG_BUILD_STRING}.json)
if [[ ${test_file} ]]; then echo "found test.py file being packaged ${test_file}"; exit 1; fi
exit 0
