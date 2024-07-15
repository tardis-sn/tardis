set "PKG_CONFIG_PATH=%LIBRARY_LIB%\pkgconfig;%LIBRARY_PREFIX%\share\pkgconfig"
%PYTHON% setup.py install --with-cython
if errorlevel 1 exit 1
