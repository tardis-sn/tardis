@echo on

for /f "delims=" %%a in ('pkg-config --libs --msvc-syntax glib-2.0') do (
    %CC% /I "%LIBRARY_PREFIX%\include\glib-2.0" /I %LIBRARY_PREFIX%\lib\glib-2.0\include "/utf-8" "-D_GNU_SOURCE" "-DUNICODE" "-D_UNICODE" "-DG_DISABLE_CAST_CHECKS" "/wo4057" "/wd4068" "/wo4090" "/wd4100" "/wd4116" "/wo4125" "/wd4127" "/wd4146" "/wd4152" "/wd4201" "/wd4232" "/wo4245" "/wo4267" "/wd4334" "/wo4389" "/wo4702" "/wd4706" /Fe:output.exe test.c /link /MACHINE:x64 "/release" "/nologo" "/OPT:REF" %%a "/SUBSYSTEM:CONSOLE" "kernel32.lib" "user32.lib" "gdi32.lib" "winspool.lib" "shell32.lib" "ole32.lib" "oleaut32.lib" "uuid.lib" "comdlg32.lib" "advapi32.lib"
    if errorlevel 1 exit 1
)

output.exe
if errorlevel 1 exit 1
