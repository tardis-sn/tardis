#!/usr/bin/env bash

set -euo pipefail

MIRROR_URL="https://dev.azure.com/tardis-sn/TARDIS/_apis/git/repositories/tardis-refdata"
HDF5_FILES=("arepo_data/arepo_snapshot.hdf5"
            "arepo_data/arepo_snapshot.json"
            "atom_data/chianti_He.h5"
            "atom_data/kurucz_cd23_chianti_H_He.h5"
            "montecarlo_1e5_compare_data.h5"
            "nlte_atom_data/TestNLTE_He_Ti.h5"
            "packet_unittest.h5"
            "sdec_ref.h5"
            "unit_test_data.h5")

REFDATA_PATH="$GITHUB_WORKSPACE/tardis-refdata"
mkdir -p "$REFDATA_PATH/arepo_data"
mkdir -p "$REFDATA_PATH/atom_data"

if [[ -z $REFDATA_SHA ]]; then
    VERSION="master"
    TYPE="branch"

else
    VERSION="$REFDATA_SHA"
    TYPE="commit"

fi

for FILE in "${HDF5_FILES[@]}"
do
    curl -o "$REFDATA_PATH/$FILE" "$MIRROR_URL/items?path=$FILE&versionType=$TYPE&version=$VERSION&resolveLfs=true"
done

exit 0
