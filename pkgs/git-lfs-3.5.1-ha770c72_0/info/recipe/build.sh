#!/usr/bin/env bash

set -o xtrace -o errexit

#bash -x install.sh

# from https://github.com/git-lfs/git-lfs/blob/main/script/install.sh
mkdir -p "${PREFIX}/bin"
for g in git*; do
  install $g "${PREFIX}/bin/${g}"
done
