if(ASTDYN_BUILD_PYTHON)
    pybind11_add_module(pyastdyn bindings.cpp)
    target_link_libraries(pyastdyn PRIVATE astdyn Eigen3::Eigen)

    # Set output location to root for easy testing
    set_target_properties(pyastdyn PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}"
    )
endif()
