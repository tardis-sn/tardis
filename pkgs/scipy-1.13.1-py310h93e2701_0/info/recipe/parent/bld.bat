@echo on
setlocal enabledelayedexpansion

mkdir builddir

:: check if clang-cl is on path as required
clang-cl.exe --version
if %ERRORLEVEL% neq 0 exit 1

:: set compilers to clang-cl
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
    -Csetup-args=-Dblas=blas ^
    -Csetup-args=-Dcpp_std=c++17 ^
    -Csetup-args=-Dlapack=lapack ^
    -Csetup-args=-Dfortran_std=none ^
    -Csetup-args=-Duse-g77-abi=true
if %ERRORLEVEL% neq 0 (type builddir\meson-logs\meson-log.txt && exit 1)
