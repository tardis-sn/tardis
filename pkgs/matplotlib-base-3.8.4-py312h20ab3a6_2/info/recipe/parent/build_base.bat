set LIB=%LIBRARY_LIB%;%LIB%
set LIBPATH=%LIBRARY_LIB%;%LIBPATH%
set INCLUDE=%LIBRARY_INC%;%INCLUDE%

ECHO [directories] > mplsetup.cfg
ECHO basedirlist = %LIBRARY_PREFIX% >> mplsetup.cfg
ECHO [packages] >> mplsetup.cfg
ECHO tests = False >> mplsetup.cfg
ECHO sample_data = False >> mplsetup.cfg
ECHO toolkits_tests = False >> mplsetup.cfg
ECHO [libs] >> mplsetup.cfg
ECHO system_freetype = True >> mplsetup.cfg

set MPLSETUPCFG=mplsetup.cfg

%PYTHON% -m pip install --no-deps --no-build-isolation -vv .
if errorlevel 1 exit 1
