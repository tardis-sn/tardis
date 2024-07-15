#!/bin/bash
set -ex

echo $MESON_ARGS

if [[ -f "$BUILD_PREFIX/meson_cross_file.txt" ]];
then
    # See https://conda-forge.org/blog/posts/2020-10-29-macos-arm64/
    # about MESON_ARGS and cross file.
    cat $BUILD_PREFIX/meson_cross_file.txt

    # Based on scipy-feedstock:
    # HACK: extend $CONDA_PREFIX/meson_cross_file that's created in
    # https://github.com/conda-forge/ctng-compiler-activation-feedstock/blob/main/recipe/activate-gcc.sh
    # https://github.com/conda-forge/clang-compiler-activation-feedstock/blob/main/recipe/activate-clang.sh
    # to use host python; requires that [binaries] section is last in meson_cross_file
    echo "python = '${PREFIX}/bin/python'" >> $BUILD_PREFIX/meson_cross_file.txt
    cat $BUILD_PREFIX/meson_cross_file.txt

    $PYTHON -m pip install . -vv --no-build-isolation \
        --config-settings=setup-args=--cross-file=$BUILD_PREFIX/meson_cross_file.txt
else
    $PYTHON -m pip install . -vv --no-build-isolation
fi
