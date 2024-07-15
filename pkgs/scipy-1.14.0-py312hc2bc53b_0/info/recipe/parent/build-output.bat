@echo on
setlocal enabledelayedexpansion

REM Set a few environment variables that are not set due to
REM https://github.com/conda/conda-build/issues/3993
set PIP_NO_BUILD_ISOLATION=True
set PIP_NO_DEPENDENCIES=True
set PIP_IGNORE_INSTALLED=True
set PIP_NO_INDEX=True
set PYTHONDONTWRITEBYTECODE=True

REM delete tests from baseline output "scipy"
if "%PKG_NAME%"=="scipy" (
    REM `pip install dist\numpy*.whl` does not work on windows,
    REM so use a loop; there's only one wheel in dist/ anyway
    for /f %%f in ('dir /b /S .\dist') do (
        pip install %%f
        if %ERRORLEVEL% neq 0 exit 1
    )

    REM copy "test" with informative error message into installation
    copy %RECIPE_DIR%\test_conda_forge_packaging.py %SP_DIR%\scipy\_lib

    REM clean up dist folder for building tests
    rmdir /s /q dist
) else (
    REM copy of build.sh, except different build tags and fortran setup;
    REM instead of using the same script (lightly templated on tags) per output,
    REM we keep the global build to reuse the cache when building the tests.

    REM need to rename project as well; for more details see
    REM https://scipy.github.io/devdocs/building/redistributable_binaries.html
    sed -i "s:name = \"scipy\":name = \"scipy-tests\":g" pyproject.toml

    REM set compilers to clang-cl
    set "CC=clang-cl"
    set "CXX=clang-cl"

    :: flang 17 still uses "temporary" name
    set "FC=flang-new"

    :: need to read clang version for path to compiler-rt
    FOR /F "tokens=* USEBACKQ" %%F IN (`clang.exe -dumpversion`) DO (
        SET "CLANG_VER=%%F"
    )

    :: attempt to match flags for flang as we set them for clang-on-win, see
    :: https://github.com/conda-forge/clang-win-activation-feedstock/blob/main/recipe/activate-clang_win-64.sh
    :: however, -Xflang --dependent-lib=msvcrt currently fails as an unrecognized option, see also
    :: https://github.com/llvm/llvm-project/issues/63741
    set "FFLAGS=-D_CRT_SECURE_NO_WARNINGS -D_MT -D_DLL --target=x86_64-pc-windows-msvc -nostdlib"
    set "LDFLAGS=--target=x86_64-pc-windows-msvc -nostdlib -Xclang --dependent-lib=msvcrt -fuse-ld=lld"
    set "LDFLAGS=%LDFLAGS% -Wl,-defaultlib:%BUILD_PREFIX%/Library/lib/clang/!CLANG_VER:~0,2!/lib/windows/clang_rt.builtins-x86_64.lib"

    :: see explanation here:
    :: https://github.com/conda-forge/scipy-feedstock/pull/253#issuecomment-1732578945
    set "MESON_RSP_THRESHOLD=320000"

    :: -wnx flags mean: --wheel --no-isolation --skip-dependency-check
    %PYTHON% -m build -w -n -x ^
        -Cbuilddir=builddir ^
        -Cinstall-args=--tags=tests ^
        -Csetup-args=-Dblas=blas ^
        -Csetup-args=-Dcpp_std=c++17 ^
        -Csetup-args=-Dlapack=lapack ^
        -Csetup-args=-Dfortran_std=none ^
        -Csetup-args=-Duse-g77-abi=true
    if %ERRORLEVEL% neq 0 (type builddir\meson-logs\meson-log.txt && exit 1)

    for /f %%f in ('dir /b /S .\dist') do (
        pip install %%f
        if %ERRORLEVEL% neq 0 exit 1
    )
)
