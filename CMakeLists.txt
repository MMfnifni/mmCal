cmake_minimum_required(VERSION 3.20)
project(mmCal VERSION 1.5.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(GLOB_RECURSE MMCAL_CORE_SOURCES CONFIGURE_DEPENDS "src/*.cpp")
list(REMOVE_ITEM MMCAL_CORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/cli/main.cpp")

add_library(mmCal.Core STATIC ${MMCAL_CORE_SOURCES})
target_include_directories(mmCal.Core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_executable(mmCal src/cli/main.cpp)
target_link_libraries(mmCal PRIVATE mmCal.Core)
target_include_directories(mmCal PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")


add_executable(mmCal.Benchmarks EXCLUDE_FROM_ALL benchmarks/benchmark_main.cpp)
target_link_libraries(mmCal.Benchmarks PRIVATE mmCal.Core)
target_include_directories(mmCal.Benchmarks PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks")

file(GLOB_RECURSE MMCAL_TEST_SOURCES CONFIGURE_DEPENDS "tests/*.cpp")
add_executable(mmCal.Tests ${MMCAL_TEST_SOURCES})
target_link_libraries(mmCal.Tests PRIVATE mmCal.Core)
target_include_directories(mmCal.Tests PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests")

if(MSVC)
    target_compile_options(mmCal.Core PRIVATE /W4 /utf-8)
    target_compile_options(mmCal PRIVATE /W4 /utf-8)
    target_compile_options(mmCal.Tests PRIVATE /W4 /utf-8)
    target_compile_options(mmCal.Benchmarks PRIVATE /W4 /utf-8)
else()
    target_compile_options(mmCal.Core PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(mmCal PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(mmCal.Tests PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(mmCal.Benchmarks PRIVATE -Wall -Wextra -Wpedantic)
endif()

enable_testing()
add_test(NAME mmCal.Tests COMMAND mmCal.Tests)
