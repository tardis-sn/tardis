#!/bin/bash

set -ex

source ${RECIPE_DIR}/setup_compiler.sh

# ensure patch is applied
grep 'conda-forge:: allow' gcc/gcc.c*

GCC_CONFIGURE_OPTIONS=()

if [[ "$channel_targets" == *conda-forge* ]]; then
  GCC_CONFIGURE_OPTIONS+=(--with-pkgversion="conda-forge gcc ${gcc_version}-${PKG_BUILDNUM}")
  GCC_CONFIGURE_OPTIONS+=(--with-bugurl="https://github.com/conda-forge/ctng-compilers-feedstock/issues/new/choose")
fi

source $RECIPE_DIR/get_cpu_arch.sh

for tool in addr2line ar as c++filt cc c++ fc gcc g++ gfortran ld nm objcopy objdump ranlib readelf size strings strip; do
  tool_upper=$(echo $tool | tr a-z A-Z | sed "s/+/X/g")
  if [[ "$tool" == "cc" ]]; then
     tool=gcc
  elif [[ "$tool" == "fc" ]]; then
     tool=gfortran
  elif [[ "$tool" == "c++" ]]; then
     tool=g++
  fi
  eval "export ${tool_upper}_FOR_BUILD=\$BUILD_PREFIX/bin/\$BUILD-\$tool"
  eval "export ${tool_upper}=\$BUILD_PREFIX/bin/\$HOST-\$tool"
  eval "export ${tool_upper}_FOR_TARGET=\$BUILD_PREFIX/bin/\$TARGET-\$tool"
done

# workaround a bug in gcc build files when using external binutils
# and build != host == target
export gcc_cv_objdump=$OBJDUMP_FOR_TARGET

ls $BUILD_PREFIX/bin/

./contrib/download_prerequisites

# We want CONDA_PREFIX/usr/lib not CONDA_PREFIX/usr/lib64 and this
# is the only way. It is incompatible with multilib (obviously).
TINFO_FILES=$(find . -path "*/config/*/t-*")
for TINFO_FILE in ${TINFO_FILES}; do
  echo TINFO_FILE ${TINFO_FILE}
  sed -i.bak 's#^\(MULTILIB_OSDIRNAMES.*\)\(lib64\)#\1lib#g' ${TINFO_FILE}
  rm -f ${TINFO_FILE}.bak
  sed -i.bak 's#^\(MULTILIB_OSDIRNAMES.*\)\(libx32\)#\1lib#g' ${TINFO_FILE}
  rm -f ${TINFO_FILE}.bak
done

# workaround for https://gcc.gnu.org/bugzilla//show_bug.cgi?id=80196
if [[ "$gcc_version" == "11."* && "$build_platform" != "$target_platform" ]]; then
  sed -i.bak 's@-I$glibcxx_srcdir/libsupc++@-I$glibcxx_srcdir/libsupc++ -nostdinc++@g' libstdc++-v3/configure
fi

mkdir -p build
cd build

# We need to explicitly set the gxx include dir because previously
# with ct-ng, native build was not considered native because
# BUILD=HOST=x86_64-build_unknown-linux-gnu and TARGET=x86_64-conda-linux-gnu
# Depending on native or not, the include dir changes. Setting it explictly
# goes back to the original way.
# See https://github.com/gcc-mirror/gcc/blob/16e2427f50c208dfe07d07f18009969502c25dc8/gcc/configure.ac#L218

../configure \
  --prefix="$PREFIX" \
  --with-slibdir="$PREFIX/lib" \
  --libdir="$PREFIX/lib" \
  --mandir="$PREFIX/man" \
  --build=$BUILD \
  --host=$HOST \
  --target=$TARGET \
  --enable-default-pie \
  --enable-languages=c,c++,fortran,objc,obj-c++ \
  --enable-__cxa_atexit \
  --disable-libmudflap \
  --enable-libgomp \
  --disable-libssp \
  --enable-libquadmath \
  --enable-libquadmath-support \
  --enable-libsanitizer \
  --enable-lto \
  --enable-threads=posix \
  --enable-target-optspace \
  --enable-plugin \
  --enable-gold \
  --disable-nls \
  --disable-bootstrap \
  --disable-multilib \
  --enable-long-long \
  --with-sysroot=${PREFIX}/${TARGET}/sysroot \
  --with-build-sysroot=${PREFIX}/${TARGET}/sysroot \
  --with-gxx-include-dir="${PREFIX}/${TARGET}/include/c++/${gcc_version}" \
  "${GCC_CONFIGURE_OPTIONS[@]}"

make -j${CPU_COUNT} || (cat ${TARGET}/libgcc/config.log; false)
