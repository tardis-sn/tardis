@echo off

"%PREFIX%\Scripts\jupyter-nbclassic-extension.exe" enable qgrid --py --sys-prefix > NUL 2>&1 && if errorlevel 1 exit 1
