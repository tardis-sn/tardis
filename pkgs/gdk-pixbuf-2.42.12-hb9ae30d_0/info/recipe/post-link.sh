#!/bin/bash

# Make sure that gdk-pixbuf's loaders.cache is fully up-to-date.
#
# Packages that install gdk-pixbuf loaders (such as librsvg) should have
# post-link and post-unlink scripts that just execute this one, which will be
# available as `$PREFIX/bin/.gdk-pixbuf-post-link.sh`.

set -e

# When cross-compiling, or installing for a different platform, the gdk-pixbuf-query-loaders binary can't be executed
# https://github.com/conda-forge/gdk-pixbuf-feedstock/issues/23
"$PREFIX/bin/gdk-pixbuf-query-loaders" --update-cache 2>>"${PREFIX}/.messages.txt" || \
(
    echo "ERROR: Failed to update gdk-pixbuf's cache, some plugins may not be found."
    echo "To fix this, activate the environment and run:"
    echo "    gdk-pixbuf-query-loaders --update-cache"
) >> "${PREFIX}/.messages.txt"
