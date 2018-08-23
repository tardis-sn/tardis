export PYTHONIOENCODING=UTF8
if test -e $HOME/miniconda/bin; then
     echo "miniconda already installed.";
else
    wget $MINICONDA_URL -O miniconda.sh
    chmod +x miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH=$HOME/miniconda/bin:$PATH
    hash -r
    conda update --yes conda
fi