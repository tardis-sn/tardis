if [[ "$blas_impl" == "openblas" ]]; then
  ln -s $PREFIX/lib/pkgconfig/openblas.pc $PREFIX/lib/pkgconfig/blas.pc
  ln -s $PREFIX/lib/pkgconfig/openblas.pc $PREFIX/lib/pkgconfig/cblas.pc
  ln -s $PREFIX/lib/pkgconfig/openblas.pc $PREFIX/lib/pkgconfig/lapack.pc
fi
