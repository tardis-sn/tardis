set HDF5_VERSION=%hdf5%
set HDF5_DIR="%LIBRARY_PREFIX%"

REM tell setup.py to not 'pip install' exact package requirements
set H5PY_SETUP_REQUIRES=0

"%PYTHON%" -m pip install . --no-deps --ignore-installed --no-cache-dir -vv
if errorlevel 1 exit 1
