

set -ex



conda inspect linkages -p $PREFIX $PKG_NAME
fribidi -h
exit 0
