#!/bin/bash
set -ex

# Set a few environment variables that are not set due to
# https://github.com/conda/conda-build/issues/3993
export PIP_NO_BUILD_ISOLATION=True
export PIP_NO_DEPENDENCIES=True
export PIP_IGNORE_INSTALLED=True
export PIP_NO_INDEX=True
export PYTHONDONTWRITEBYTECODE=True

# need to use force to reinstall the tests the second time
# (otherwise pip thinks the package is installed already)
pip install dist/scipy*.whl --force-reinstall

# delete tests from baseline output "scipy"
if [[ "$PKG_NAME" == "scipy" ]]; then
    # verify $RECIPE_DIR/test_folders_to_delete.txt is up to date;
    # it's enough to do this on one platform
    if [[ "${target_platform}" == "linux-64" ]]; then
        # validating this is important because windows does not have
        # a good dynamic command like 'find ... -name tests -type d',
        # and we're using this file to do the deletions on windows.
        find ${SP_DIR}/scipy -name tests -type d -printf '%p\n' \
            | sort -k1 | sed "s:${SP_DIR}/scipy/::g" > testfolders
        echo "Test folders to be deleted:"
        cat testfolders
        # diff returns error code if there are differences
        diff $RECIPE_DIR/test_folders_to_delete.txt testfolders

        # same procedure for extra test DLLs/SOs; as above, but additionally, replace
        # ABI tag with a marker (here it's helpful this branch is only for linux-64)
        find ${SP_DIR}/scipy -regex ".*_c?y?test.*\.so" -printf '%p\n' \
            | sort -k1 | sed "s:${SP_DIR}/scipy/::g" \
            | sed "s:.cpython-${PY_VER//./}-x86_64-linux-gnu.so:SUFFIX_MARKER:g" \
            > testlibs
        echo "Test libraries to be deleted:"
        cat testlibs
        if [[ $python_impl == "cpython" ]] ; then
            # don't try on pypy; has different tag so the above sed doesn't apply
            diff $RECIPE_DIR/test_libraries_to_delete.txt testlibs
        fi
    fi

    # do the actual deletion
    find ${SP_DIR}/scipy -name tests -type d | xargs rm -r
    # different syntax for regex on osx ('\x2d\x45'=='-E'; else echo insists it's a flag)
    OSX_EXTRA="$([[ ${target_platform} == osx-* ]] && echo -e '\x2d\x45' || echo '')"
    find ${OSX_EXTRA} ${SP_DIR}/scipy -regex ".*_c?y?test.*\.so" | xargs rm

    # copy "test" with informative error message into installation
    cp $RECIPE_DIR/test_conda_forge_packaging.py $SP_DIR/scipy/_lib
fi
