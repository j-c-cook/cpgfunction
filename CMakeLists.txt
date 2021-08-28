cmake_minimum_required(VERSION 3.6)
project(CpGFunction)

set(CMAKE_CXX_STANDARD 11)
# Compiler requests -std=c++11 be set for compiler option
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Everything will need to access this
include_directories(include)
include_directories(third_party)
include_directories(third_party/LinearAlgebra)
include_directories(third_party/LinearAlgebra/include)
add_subdirectory(third_party/LinearAlgebra)

# Create a DLL with an API
add_library(cpgfunction SHARED src/boreholes.cpp
        src/gfunction.cpp
        src/heat_transfer.cpp
        src/interpolation.cpp
        src/utilities.cpp
        third_party/nlohmann/json.hpp
        src/coordinates.cpp
        src/statistics.cpp
        third_party/LinearAlgebra/src/axpy.cpp
        third_party/LinearAlgebra/src/spmv.cpp)

# It depends on boost (stripped)
find_package(Boost REQUIRED)
target_link_libraries(cpgfunction Boost::boost)

# It depends on threads
find_package(Threads REQUIRED)
target_link_libraries(cpgfunction ${CMAKE_THREAD_LIBS_INIT})

#find_package(LAPACK)
#target_link_libraries(cpgfunction ${LAPACK_LIBRARIES})

# OpenBLAS is specified because the default (Netlib BLAS) is much slower for LU decomposition (dgesv_)
set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    target_link_libraries(cpgfunction ${BLAS_LIBRARIES})
else()
    # ERROR
endif()

OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
    FIND_PACKAGE(OpenMP)
    IF(OPENMP_FOUND)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    ENDIF()
ENDIF()

# Create a test binary that depends on the DLL
add_executable(test1 test/gFunction_minimal.cpp)
add_executable(test2 test/interpolation.cpp)
# add_executable(test3 test/configurations.cpp)
add_executable(test4 test/borefield_definition.cpp)
add_executable(test5 test/time_definition.cpp)
add_executable(test6 test/compute_UBHWT_gFunction.cpp)
add_executable(test7 test/packed.cpp)

target_link_libraries(test1 cpgfunction)
target_link_libraries(test2 cpgfunction)
# target_link_libraries(test3 cpgfunction)
target_link_libraries(test4 cpgfunction)
target_link_libraries(test5 cpgfunction)
target_link_libraries(test6 cpgfunction)
target_link_libraries(test7 cpgfunction)

# target_compile_definitions(cpgfunction PUBLIC TEST_RESOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
# Copy validation files to build directory so tests can open
file(GLOB JsonValidation test/validation/*.json)
file(GLOB CustomLayout test/custom_configurations/*.json)

file(COPY ${JsonValidation} DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CustomLayout} DESTINATION ${CMAKE_BINARY_DIR})
message(STATUS "The binary directory: ${CMAKE_BINARY_DIR}")

# Turn on testing
enable_testing()

# Test the test binary itself
add_test(NAME RunTest1 COMMAND "${CMAKE_BINARY_DIR}/test1")
add_test(NAME RunTest2 COMMAND "${CMAKE_BINARY_DIR}/test2")
# add_test(NAME RunTest3 COMMAND "${CMAKE_BINARY_DIR}/test3")
add_test(NAME RunTest4 COMMAND "${CMAKE_BINARY_DIR}/test4")
add_test(NAME RunTest5 COMMAND "${CMAKE_BINARY_DIR}/test5")
# Pass variable path into test 6 for json files
add_test(NAME RunTest6 COMMAND ${CMAKE_BINARY_DIR}/test6)
add_test(NAME RunTest7 COMMAND ${CMAKE_BINARY_DIR}/test7)
