#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ioc_gaialib" for configuration ""
set_property(TARGET ioc_gaialib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(ioc_gaialib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libioc_gaialib.a"
  )

list(APPEND _cmake_import_check_targets ioc_gaialib )
list(APPEND _cmake_import_check_files_for_ioc_gaialib "${_IMPORT_PREFIX}/lib/libioc_gaialib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
