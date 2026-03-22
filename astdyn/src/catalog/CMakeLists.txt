add_library(astdyn_catalog STATIC
    GaiaDR3Catalog.cpp
    CatalogIntegration.cpp
)

target_include_directories(astdyn_catalog
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(astdyn_catalog
    PUBLIC
        astdyn
    PRIVATE
        CURL::libcurl
        SQLite::SQLite3
        LibXml2::LibXml2
        ZLIB::ZLIB
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(astdyn_catalog PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS astdyn_catalog
    EXPORT AstDynTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/astdyn/catalog
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/astdyn
    FILES_MATCHING PATTERN "*.hpp"
)
