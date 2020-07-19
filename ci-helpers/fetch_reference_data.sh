#!/usr/bin/env bash
git lfs install --skip-smudge

if test -e $REF_DATA_DIR; then
    echo "Reference data available"

else
    git clone $REF_DATA_URL $REF_DATA_DIR -o upstream
    cd $REF_DATA_DIR

    # Use the following to get the ref-data from the master;
    git fetch upstream
    git checkout upstream/master

    # Use the following to get the ref-data from a specific pull request
    # git fetch upstream pull/20/head:update-ref
    # git checkout update-ref

    git lfs pull --include="atom_data/kurucz_cd23_chianti_H_He.h5" upstream
    git lfs pull --include="atom_data/chianti_He.h5" upstream
    git lfs pull --include="plasma_reference/" upstream
    git lfs pull --include="unit_test_data.h5" upstream
    git lfs pull --include="packet_unittest.h5" upstream

    echo MD5: `md5sum unit_test_data.h5`
    echo MD5: `md5sum atom_data/kurucz_cd23_chianti_H_He.h5`

fi

exit 0