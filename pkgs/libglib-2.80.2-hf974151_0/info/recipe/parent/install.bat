cd forgebuild
ninja install
if errorlevel 1 exit 1

if NOT [%PKG_NAME%] == [glib] (
  if [%PKG_NAME%] == [libglib] (
      del %LIBRARY_PREFIX%\bin\gdbus.exe
      if errorlevel 1 exit 1
      del %LIBRARY_PREFIX%\bin\gio-querymodules.exe
      if errorlevel 1 exit 1
      del %LIBRARY_PREFIX%\bin\gio.exe
      if errorlevel 1 exit 1
      del %LIBRARY_PREFIX%\bin\glib-compile-schemas.exe
      if errorlevel 1 exit 1
      del %LIBRARY_PREFIX%\bin\gresource.exe
      if errorlevel 1 exit 1
      del %LIBRARY_PREFIX%\bin\gsettings.exe
      if errorlevel 1 exit 1
  )
  del %LIBRARY_PREFIX%\bin\gdbus-codegen
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\bin\glib-compile-resources.exe
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\bin\glib-genmarshal
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\bin\glib-gettextize
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\bin\glib-mkenums
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\bin\gobject-query.exe
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\bin\gi-compile-repository.exe
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\bin\gi-inspect-typelib.exe
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\bin\gi-decompile-typelib.exe
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\bin\gtester*
  if errorlevel 1 exit 1

  rmdir /s /q %LIBRARY_PREFIX%\include\gio-win32-2.0
  if errorlevel 1 exit 1
  rmdir /s /q %LIBRARY_PREFIX%\include\glib-2.0
  if errorlevel 1 exit 1

  rmdir /s /q %LIBRARY_PREFIX%\lib\glib-2.0\include
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\lib\pkgconfig\gio-*
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\lib\pkgconfig\glib-*
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\lib\pkgconfig\gmodule-*
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\lib\pkgconfig\gobject-*
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\lib\pkgconfig\gthread-*
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\lib\pkgconfig\girepository-*
  if errorlevel 1 exit 1

  del %LIBRARY_PREFIX%\share\aclocal\glib-*
  if errorlevel 1 exit 1
  del %LIBRARY_PREFIX%\share\aclocal\gsettings.m4
  if errorlevel 1 exit 1
  rmdir /s /q %LIBRARY_PREFIX%\share\gettext\its
  if errorlevel 1 exit 1
  rmdir /s /q %LIBRARY_PREFIX%\share\glib-2.0
  if errorlevel 1 exit 1
) else (
  @rem intl.lib is currently statically linked and thus should not be mentioned as a hard dependency.
  @rem Remove this once we are linking it dynamically.
  for %%f in (%LIBRARY_PREFIX%\lib\pkgconfig\*.pc) do (
    sed -i "s/-lintl//g" %%f
  )
)

rem We don't have bash as a dependency so these shouldn't exist, but
rem sometimes a system bash will be picked up and they will get installed.
rem Just delete them, but don't check for errors in case they do not exist.
rem If we do want them in the future, they should go in glib-tools.
del %LIBRARY_PREFIX%\share\bash-completion\completions\gapplication
del %LIBRARY_PREFIX%\share\bash-completion\completions\gdbus
del %LIBRARY_PREFIX%\share\bash-completion\completions\gio
del %LIBRARY_PREFIX%\share\bash-completion\completions\gresource
del %LIBRARY_PREFIX%\share\bash-completion\completions\gsettings

del %LIBRARY_PREFIX%\bin\*.pdb
