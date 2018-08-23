git lfs install
if test -e $REF_DATA_HOME; then
    echo "Ref data available"
else
    git clone https://github.com/tardis-sn/tardis-refdata.git $REF_DATA_HOME
cd $REF_DATA_HOME
# Use the following to get the ref-data from the master;
git fetch origin
git checkout origin/master
# Use the following to get the ref-data from a specific pull request
#git fetch origin pull/11/head:thomson-ref
git lfs pull --include="atom_data/kurucz_cd23_chianti_H_He.h5" origin/master
git lfs pull --include="atom_data/chianti_He.h5" origin/master
git lfs pull --include="plasma_reference/" origin/master
git lfs pull --include="unit_test_data.h5" origin/master
cd $TRAVIS_BUILD_DIR
