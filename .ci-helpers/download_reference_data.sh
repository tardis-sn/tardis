#!/usr/bin/env bash

set -e

REF_PATH="$GITHUB_WORKSPACE/tardis-refdata"
REPO_URL="https://dev.azure.com/tardis-sn/TARDIS/_apis/git/repositories/tardis-refdata"

mkdir -p $REF_PATH/atom_data
wget -q "$REPO_URL/items?path=atom_data/kurucz_cd23_chianti_H_He.h5&resolveLfs=true" -O $REF_PATH/atom_data/kurucz_cd23_chianti_H_He.h5
wget -q "$REPO_URL/items?path=atom_data/chianti_He.h5&resolveLfs=true" -O $REF_PATH/atom_data/chianti_He.h5
wget -q "$REPO_URL/items?path=unit_test_data.h5&resolveLfs=true" -O $REF_PATH/unit_test_data.h5&resolveLfs=true
wget -q "$REPO_URL/items?path=packet_unittest.h5&resolveLfs=true" -O $REF_PATH/packet_unittest.h5&resolveLfs=true
wget -q "$REPO_URL/items?path=montecarlo_1e5_compare_data.h5&resolveLfs=true" -O $REF_PATH/montecarlo_1e5_compare_data.h5

exit 0
