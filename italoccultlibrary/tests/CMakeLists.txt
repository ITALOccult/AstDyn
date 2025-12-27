# Tests for ITALOccultLibrary

find_package(GTest REQUIRED)
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/../astdyn/include
    ${EIGEN3_INCLUDE_DIR}
    ${GTEST_INCLUDE_DIRS}
)

# Test for ChebyshevApproximation
add_executable(test_chebyshev test_chebyshev.cpp)
target_link_libraries(test_chebyshev 
    italoccultlib 
    GTest::GTest 
    GTest::Main
    Eigen3::Eigen
)

# Test for AstDynWrapper
add_executable(test_astdyn_wrapper test_astdyn_wrapper.cpp)
target_link_libraries(test_astdyn_wrapper 
    italoccultlib 
    GTest::GTest 
    GTest::Main
    Eigen3::Eigen
)

# Test AstDynWrapper precision (34713 Ilse)
add_executable(test_34713_precision test_34713_precision.cpp)
# Link against the local astdyn library explicitly
target_link_libraries(test_34713_precision PRIVATE 
    italoccultlib 
    ${PROJECT_SOURCE_DIR}/../astdyn/build/src/libastdyn.a
    ${Boost_LIBRARIES}
)
add_test(NAME test_34713_precision COMMAND test_34713_precision)

add_test(NAME ChebyshevTest COMMAND test_chebyshev)
add_test(NAME AstDynWrapperTest COMMAND test_astdyn_wrapper)
