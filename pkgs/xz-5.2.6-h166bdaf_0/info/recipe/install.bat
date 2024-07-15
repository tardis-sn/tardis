ninja install
if errorlevel 1 exit /b 1

if %PKG_NAME%==xz-static (
	echo "Keeping all files, conda will dedupe"
) else (
	DEL src\common\inttypes.h
	DEL src\common\stdint.h
	DEL %LIBRARY_LIB%\liblzma_static.lib

	cd %SRC_DIR%
	MOVE src\liblzma\api\lzma %LIBRARY_INC%\
	COPY src\liblzma\api\lzma.h %LIBRARY_INC%\
	exit /b 0
)
