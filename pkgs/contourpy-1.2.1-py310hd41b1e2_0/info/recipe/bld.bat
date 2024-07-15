@echo on

mkdir builddir
%PYTHON% -m pip install . -vv --no-build-isolation -Cbuilddir=builddir
if %ERRORLEVEL% neq 0 (type builddir\meson-logs\meson-log.txt && exit 1)
