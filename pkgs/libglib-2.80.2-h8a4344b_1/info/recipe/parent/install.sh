#!/bin/bash

set -exuo pipefail

unset _CONDA_PYTHON_SYSCONFIGDATA_NAME

ninja -C builddir install || (cat meson-logs/meson-log.txt; false)
# remove libtool files
find $PREFIX -name '*.la' -delete

# gdb folder has a nested folder structure similar to our host prefix
# (255 chars) which causes installation issues so remove it.
rm -rf $PREFIX/share/gdb

if [[ "$PKG_NAME" != glib ]]; then
    if [[ "$PKG_NAME" == glib-tools ]]; then
        mkdir .keep
        # We ship these binaries as part of the glib-tools package because
        # they can be depended on separately by other packages, e.g. gtk
        # (equivalent to Debian's libglib2.0-bin)
        mv $PREFIX/bin/gdbus .keep
        mv $PREFIX/bin/glib-compile-schemas .keep
    else
        rm -f $PREFIX/bin/gapplication
        rm $PREFIX/bin/gio*
        rm $PREFIX/bin/gresource
        rm $PREFIX/bin/gsettings
        rm $PREFIX/share/bash-completion/completions/{gapplication,gdbus,gio,gresource,gsettings}

        # Copy the [de]activate scripts to $PREFIX/etc/conda/[de]activate.d.
        # This will allow them to be run on environment activation.
        for CHANGE in "activate" "deactivate"
        do
            mkdir -p "${PREFIX}/etc/conda/${CHANGE}.d"
            cp "${RECIPE_DIR}/scripts/${CHANGE}.sh" "${PREFIX}/etc/conda/${CHANGE}.d/${PKG_NAME}_${CHANGE}.sh"
        done
    fi
    rm $PREFIX/bin/{gdbus*,glib-*,gobject*,gtester*,gi-*}
    if [[ "$PKG_NAME" == glib-tools ]]; then
        mv .keep/* $PREFIX/bin
    fi
    rm -r $PREFIX/include/gio-* $PREFIX/include/glib-*
    rm -r $PREFIX/lib/glib-*
    rm -r $PREFIX/lib/lib{gmodule,glib,gobject,gthread,gio}-2.0${SHLIB_EXT}
    rm -r $PREFIX/share/aclocal/{glib-*,gsettings*}
    rm -r $PREFIX/share/gettext/its
    rm -r $PREFIX/share/glib-*
    # Manually install introspection data during cross-compilation
    # These files are the only difference when running with a different setting of -Dintrospection
    if [[ "${CONDA_BUILD_CROSS_COMPILATION:-0}" == 1 ]]; then
        cp -ap introspection/lib/girepository-1.0 $PREFIX/lib
        cp -ap introspection/share/gir-1.0 $PREFIX/share
    fi
fi
