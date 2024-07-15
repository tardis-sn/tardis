@echo on
setlocal enabledelayedexpansion

REM Set a few environment variables that are not set due to
REM https://github.com/conda/conda-build/issues/3993
set PIP_NO_BUILD_ISOLATION=True
set PIP_NO_DEPENDENCIES=True
set PIP_IGNORE_INSTALLED=True
set PIP_NO_INDEX=True
set PYTHONDONTWRITEBYTECODE=True

:: `pip install dist\numpy*.whl` does not work on windows,
:: so use a loop; there's only one wheel in dist/ anyway
for /f %%f in ('dir /b /S .\dist') do (
    REM need to use force to reinstall the tests the second time
    REM (otherwise pip thinks the package is installed already)
    pip install %%f --force-reinstall
    if %ERRORLEVEL% neq 0 exit 1
)

FOR /F "tokens=* USEBACKQ" %%F IN (`python -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"`) DO (
    SET EXT_SUFFIX=%%F
)

REM delete tests from baseline output "scipy"
if "%PKG_NAME%"=="scipy" (
    REM folders in test_folders_to_delete.txt are relative to %SP_DIR%\scipy
    REM the contents of this are validated in build-output.sh
    for /F %%f in (%RECIPE_DIR%\test_folders_to_delete.txt) do (
        set "g=%%f"
        rmdir /s /q %SP_DIR%\scipy\!g:/=\!
    )
    REM additionally delete folder not found on linux
    rmdir /s /q %SP_DIR%\scipy\_build_utils\tests

    REM same for test_libraries_to_delete.txt
    for /F %%f in (%RECIPE_DIR%\test_libraries_to_delete.txt) do (
        set "g=%%f"
        REM replace suffix marker with python ABI
        set "h=!g:SUFFIX_MARKER=%EXT_SUFFIX%!"
        del /f %SP_DIR%\scipy\!h:/=\!
    )

    REM copy "test" with informative error message into installation
    copy %RECIPE_DIR%\test_conda_forge_packaging.py %SP_DIR%\scipy\_lib
)
