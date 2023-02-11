#!/usr/bin/env bash

set -e

REF_PATH="$GITHUB_WORKSPACE/tardis-refdata"
REPO_URL="https://dev.azure.com/tardis-sn/TARDIS/_apis/git/repositories/tardis-refdata"

FILES=('atom_data/kurucz_cd23_chianti_H_He.h5'
       'atom_data/chianti_He.h5'
       'montecarlo_1e5_compare_data.h5'
       'packet_unittest.h5'
       'sdec_ref.h5'
       'unit_test_data.h5'
       'arepo_data/arepo_snapshot.hdf5'
       'arepo_data/arepo_snapshot.json'
       'nlte_atom_data/TestNLTE_He_Ti.h5')

mkdir -p $REF_PATH/atom_data
mkdir -p $REF_PATH/arepo_data
mkdir -p $REF_PATH/nlte_atom_data

for FILE in "${FILES[@]}"
do
    wget -q "$REPO_URL/items?path=$FILE&resolveLfs=true" -O $REF_PATH/$FILE
done

exit 0
