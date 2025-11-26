#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "AstDyn::astdyn" for configuration "Release"
set_property(TARGET AstDyn::astdyn APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(AstDyn::astdyn PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libastdyn.a"
  )

list(APPEND _cmake_import_check_targets AstDyn::astdyn )
list(APPEND _cmake_import_check_files_for_AstDyn::astdyn "${_IMPORT_PREFIX}/lib/libastdyn.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
