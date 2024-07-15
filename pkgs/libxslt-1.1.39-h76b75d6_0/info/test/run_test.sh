

set -ex



xsltproc --version
pkg-config --cflags libxslt libexslt
pkg-config --libs libxslt libexslt
exit 0
