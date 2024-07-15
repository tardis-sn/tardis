#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Snappy::snappy" for configuration "Release"
set_property(TARGET Snappy::snappy APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Snappy::snappy PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsnappy.so.1.2.0"
  IMPORTED_SONAME_RELEASE "libsnappy.so.1"
  )

list(APPEND _cmake_import_check_targets Snappy::snappy )
list(APPEND _cmake_import_check_files_for_Snappy::snappy "${_IMPORT_PREFIX}/lib/libsnappy.so.1.2.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
