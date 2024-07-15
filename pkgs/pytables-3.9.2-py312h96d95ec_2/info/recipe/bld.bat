set BLOSC_DIR=%LIBRARY_PREFIX%
set BLOSC2_DIR=%LIBRARY_PREFIX%
set BZIP2_DIR=%LIBRARY_PREFIX%
set HDF5_DIR=%LIBRARY_PREFIX%
set LZO_DIR=%LIBRARY_PREFIX%
set COPY_DLLS=FALSE
set PYTABLES_NO_BLOSC2_WHEEL=TRUE

del /F /Q .\tables\*.c

%PYTHON% setup.py install --single-version-externally-managed --record record.txt
if errorlevel 1 exit 1
