conda env create -f tardis_env27.yml
source activate tardis

if [[ $SETUP_CMD == *coverage* ]]; then
    conda install -c conda-forge git-lfs=2.2.1 -y
    git lfs install --skip-smudge
    if [ ! -d "$HOME/tardis-refdata" ]; then
        git clone $TARDIS_REF_DATA_URL $HOME/tardis-refdata
    fi
    cd $HOME/tardis-refdata

    # Checkout PR version of tardis-refdata, not master
    git fetch origin pull/3/head:carsus-ref
    git checkout carsus-ref

    # Pull the files
    git lfs pull --include="atom_data/kurucz_cd23_chianti_H_He.h5" origin
    git lfs pull --include="atom_data/chianti_He.h5" origin
    git lfs pull --include="plasma_reference/" origin
    git lfs pull --include="unit_test_data.h5" origin
    cd $TRAVIS_BUILD_DIR
fi
