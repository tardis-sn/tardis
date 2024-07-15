set LIB=%LIBRARY_LIB%;%LIB%
set INCLUDE=%LIBRARY_INC%;%INCLUDE%
set LIBRARY_DIRS=%LIBRARY_BIN%;%LIBRARY_LIB%

set JPEG_ROOT=%LIBRARY_PREFIX%
set JPEG2K_ROOT=%LIBRARY_PREFIX%
set ZLIB_ROOT=%LIBRARY_PREFIX%
:: currently disabled, see meta.yaml
:: set IMAGEQUANT_ROOT=%LIBRARY_PREFIX%
set TIFF_ROOT=%LIBRARY_PREFIX%
set FREETYPE_ROOT=%LIBRARY_PREFIX%
:: as above
:: set FRIBIDI_ROOT=%LIBRARY_PREFIX%
set LCMS_ROOT=%LIBRARY_PREFIX%
set WEBP_ROOT=%LIBRARY_PREFIX%
set XCB_ROOT=%LIBRARY_PREFIX%

:: add --vendor-raqm to installation (cannot be passed through pip install)
echo [build_ext] >> setup.cfg
echo vendor-raqm=1 >> setup.cfg
:: sanity check
type setup.cfg

:: debug
echo "Search for imagequant.dll, part 1"
find "D:\bld" -iname imagequant.dll
echo "Search for imagequant.dll, part 2"
dir imagequant.dll /s
echo "Search for imagequant.dll, part 3"
dir %LIBRARY_PREFIX%
echo "Search for imagequant.dll, part 4"
dir %LIBRARY_LIB%
echo "Search for imagequant.dll, part 5"
dir %LIBRARY_BIN%
echo "Search complete"

%PYTHON% -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv
if errorlevel 1 exit 1
