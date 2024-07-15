#!/bin/bash

export JPEG_ROOT=$PREFIX
export JPEG2K_ROOT=$PREFIX
export ZLIB_ROOT=$PREFIX
# currently disabled, see meta.yaml
# export IMAGEQUANT_ROOT=$PREFIX
export TIFF_ROOT=$PREFIX
export FREETYPE_ROOT=$PREFIX
# export FRIBIDI_ROOT=$PREFIX
export LCMS_ROOT=$PREFIX
export WEBP_ROOT=$PREFIX
export XCB_ROOT=$PREFIX

# add --vendor-raqm to installation (cannot be passed through pip install)
echo "[build_ext]" >> setup.cfg
echo "vendor-raqm=1" >> setup.cfg
# sanity check
cat setup.cfg

$PYTHON -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv
