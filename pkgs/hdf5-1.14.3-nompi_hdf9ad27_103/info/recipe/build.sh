#!/bin/bash
set -ex

export LIBRARY_PATH="${PREFIX}/lib"

# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/libtool/build-aux/config.* ./bin

# hdf5 autogen.sh doesn't find libtool on mac correctly
export HDF5_LIBTOOL=$BUILD_PREFIX/bin/libtool

HDF5_OPTIONS=

if [[ "$target_platform" == linux-* ]]; then
    # Direct Virtual File System (O_DIRECT)
    # is only valid for linux
    HDF5_OPTIONS="${HDF5_OPTIONS} --enable-direct-vfd"
fi

if [[ ! -z "$mpi" && "$mpi" != "nompi" ]]; then
  export HDF5_OPTIONS="${HDF5_OPTIONS} --enable-parallel"

  export CC=$PREFIX/bin/mpicc
  export CXX=$PREFIX/bin/mpic++
  export FC=$PREFIX/bin/mpifort

  if [[ "$CONDA_BUILD_CROSS_COMPILATION" == "1" ]]; then
    if [[ "$mpi" == "openmpi" ]]; then
      rm $PREFIX/bin/opal_wrapper
      echo '#!/bin/bash' > $PREFIX/bin/opal_wrapper
      echo "export OPAL_PREFIX=$PREFIX" >> $PREFIX/bin/opal_wrapper
      echo "$BUILD_PREFIX/bin/\$(basename \"\$0\") \"\$@\"" >> $PREFIX/bin/opal_wrapper
      chmod +x $PREFIX/bin/opal_wrapper
    fi
    export CC_FOR_BUILD=$BUILD_PREFIX/bin/mpicc
    export CXX_FOR_BUILD=$BUILD_PREFIX/bin/mpic++
    export FC_FOR_BUILD=$BUILD_PREFIX/bin/mpifort
  else
    export CC_FOR_BUILD=$PREFIX/bin/mpicc
    export CXX_FOR_BUILD=$PREFIX/bin/mpic++
    export FC_FOR_BUILD=$PREFIX/bin/mpifort
  fi

  if [[ "$target_platform" == linux-* ]]; then
    # --as-needed appears to cause problems with fortran compiler detection
    # due to missing libquadmath
    # unclear why required libs are stripped but still linked
    export FFLAGS="${FFLAGS:-} -Wl,--no-as-needed -Wl,--disable-new-dtags"
    export LDFLAGS="${LDFLAGS} -Wl,--no-as-needed -Wl,--disable-new-dtags"
  fi
else
  export CC=$(basename ${CC})
  export CXX=$(basename ${CXX})
  export F95=$(basename ${F95})
  export FC=$(basename ${FC})
fi

if [[ "$CONDA_BUILD_CROSS_COMPILATION" == 1 && $target_platform == "osx-arm64" ]]; then
  export ac_cv_sizeof_long_double=8
  export hdf5_cv_ldouble_to_long_special=no
  export hdf5_cv_long_to_ldouble_special=no
  export hdf5_cv_ldouble_to_llong_accurate=yes
  export hdf5_cv_llong_to_ldouble_correct=yes
  export hdf5_cv_disable_some_ldouble_conv=no
  export hdf5_cv_system_scope_threads=yes
  export hdf5_cv_printf_ll="l"

  export hdf5_cv_system_scope_threads=yes
  export hdf5_cv_printf_ll="l"
  export PAC_FC_MAX_REAL_PRECISION=15
  export PAC_C_MAX_REAL_PRECISION=17
  export PAC_FC_ALL_INTEGER_KINDS="{1,2,4,8,16}"
  export PAC_FC_ALL_REAL_KINDS="{4,8}"
  export H5CONFIG_F_NUM_RKIND="INTEGER, PARAMETER :: num_rkinds = 2"
  export H5CONFIG_F_NUM_IKIND="INTEGER, PARAMETER :: num_ikinds = 5"
  export H5CONFIG_F_RKIND="INTEGER, DIMENSION(1:num_rkinds) :: rkind = (/4,8/)"
  export H5CONFIG_F_IKIND="INTEGER, DIMENSION(1:num_ikinds) :: ikind = (/1,2,4,8,16/)"
  export PAC_FORTRAN_NATIVE_INTEGER_SIZEOF="                    4"
  export PAC_FORTRAN_NATIVE_INTEGER_KIND="           4"
  export PAC_FORTRAN_NATIVE_REAL_SIZEOF="                    4"
  export PAC_FORTRAN_NATIVE_REAL_KIND="           4"
  export PAC_FORTRAN_NATIVE_DOUBLE_SIZEOF="                    8"
  export PAC_FORTRAN_NATIVE_DOUBLE_KIND="           8"
  export PAC_FORTRAN_NUM_INTEGER_KINDS="5"
  export PAC_FC_ALL_REAL_KINDS_SIZEOF="{4,8}"
  export PAC_FC_ALL_INTEGER_KINDS_SIZEOF="{1,2,4,8,16}"
  # Do we need this below?
  export hdf5_cv_szlib_can_encode=yes

  HDF5_OPTIONS="${HDF5_OPTIONS} --enable-tests=no"
fi

# regen config after patches to configure.ac
./autogen.sh

./configure --prefix="${PREFIX}" \
            --with-pic \
            --host="${HOST}" \
            --build="${BUILD}" \
            --with-zlib="${PREFIX}" \
            --with-szlib="${PREFIX}" \
            --with-pthread=yes  \
            ${HDF5_OPTIONS} \
            --enable-cxx \
            --enable-fortran \
            --with-default-plugindir="${PREFIX}/lib/hdf5/plugin" \
            --enable-threadsafe \
            --enable-build-mode=production \
            --enable-unsupported \
            --enable-hlgiftools=yes \
            --enable-using-memchecker \
            --enable-static=no \
            --enable-ros3-vfd \
            || (cat config.log; false)

# allow oversubscribing with openmpi in make check
export OMPI_MCA_rmaps_base_oversubscribe=yes

if [[ "$CONDA_BUILD_CROSS_COMPILATION" == 1 ]]; then
  (
    # Make a native build of the hdetect and H5make_libsettings executables
    mkdir -p native-build/bin
    pushd native-build/bin

    # MACOSX_DEPLOYMENT_TARGET is for the target_platform and not for build_platform
    unset MACOSX_DEPLOYMENT_TARGET

    $CC_FOR_BUILD ${SRC_DIR}/fortran/src/H5match_types.c -I ${SRC_DIR}/src/ -o H5match_types
    $FC_FOR_BUILD ${SRC_DIR}/fortran/src/H5_buildiface.F90 -I ${SRC_DIR}/fortran/src/ -L $BUILD_PREFIX/lib -o H5_buildiface
    $FC_FOR_BUILD ${SRC_DIR}/hl/fortran/src/H5HL_buildiface.F90 -I ${SRC_DIR}/hl/fortran/src -I ${SRC_DIR}/fortran/src -L $BUILD_PREFIX/lib -o H5HL_buildiface
    popd
  )
  export PATH=`pwd`/native-build/bin:$PATH
fi

make -j "${CPU_COUNT}" ${VERBOSE_AT}

make install V=1

if [[ ${mpi} == "mpich" || (${mpi} == "openmpi" && "$(uname)" == "Darwin") ]]; then
  # ph5diff hangs on darwin with openmpi, skip the test
  # ph5diff also crashes on mpich 4.1
  echo <<EOF > tools/test/h5diff/testph5diff.sh
#!/bin/sh
exit 0
EOF
fi
if [[ ("$target_platform" != "linux-ppc64le") && \
      ("$target_platform" != "linux-aarch64") && \
      ("$target_platform" != "osx-arm64") ]]; then
  # https://github.com/h5py/h5py/issues/817
  # https://forum.hdfgroup.org/t/hdf5-1-10-long-double-conversions-tests-failed-in-ppc64le/4077
  make check RUNPARALLEL="${RECIPE_DIR}/mpiexec.sh -n 2"
fi
