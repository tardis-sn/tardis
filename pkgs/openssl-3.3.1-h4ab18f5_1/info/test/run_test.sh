

set -ex



touch checksum.txt
$PREFIX/bin/openssl sha256 checksum.txt
pkg-config --print-errors --exact-version "3.3.1" openssl
if [[ "$(pkg-config --variable=prefix openssl)" == "" ]]; then exit 1; else exit 0; fi
exit 0
