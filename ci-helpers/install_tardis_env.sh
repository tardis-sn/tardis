cd $TRAVIS_BUILD_DIR
if test -e $HOME/miniconda/envs/tardis; then
    echo "TARDIS env already installed.";
else
   conda env create -f tardis_env27.yml
fi

source activate tardis