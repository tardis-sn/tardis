#!/usr/bin/env bash
set -xv

# NOTE: If the `pkgs` and `.asv/env` folders exists. It doesn't need to be deleted except,
# if you want to run `asv setup` to recreate the ASV environment.
time rm --force --recursive pkgs .asv/env
time asv setup

time asv machine --yes

rm --force --recursive .asv/results
time asv run

rm --force --recursive .asv/html
time asv publish

#mamba deactivate
#mamba activate tardis
#python setup.py develop
#pytest tardis || true

#mamba deactivate
#mamba activate benchmark
