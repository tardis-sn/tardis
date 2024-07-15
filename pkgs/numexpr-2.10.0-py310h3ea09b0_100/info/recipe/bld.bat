if "%mkl_support%" == "mkl" (
  (
  echo [mkl]
  echo library_dirs = %LIBRARY_LIB%
  echo include_dirs = %LIBRARY_INC%
  echo libraries = mkl_rt
  ) > site.cfg
)

%PYTHON% -m pip install . --no-deps -vv
