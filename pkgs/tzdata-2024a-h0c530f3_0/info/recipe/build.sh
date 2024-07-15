make -e \
  DESTDIR=./build \
  USRDIR='' \
  POSIXRULES='' \
  install

mkdir -p "${PREFIX}/share"
mv ./build/share/zoneinfo "${PREFIX}/share/"
