cmake_minimum_required(VERSION 3.20)
project(mmCal VERSION 1.5.3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(MMCAL_ENABLE_LTO "Enable interprocedural optimization/LTO in Release builds" ON)
option(MMCAL_NATIVE_OPTIMIZATION "Optimize for the CPU performing the build (-march=native)" OFF)
option(MMCAL_LINUX_OUT_SUFFIX "Name the Linux executable mmCal.out" ON)

file(GLOB_RECURSE MMCAL_CORE_SOURCES CONFIGURE_DEPENDS "src/*.cpp")
list(REMOVE_ITEM MMCAL_CORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/cli/main.cpp")

add_library(mmCal.Core STATIC ${MMCAL_CORE_SOURCES})
target_include_directories(mmCal.Core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_executable(mmCal src/cli/main.cpp)
target_link_libraries(mmCal PRIVATE mmCal.Core)
target_include_directories(mmCal PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_executable(mmCal.Benchmarks EXCLUDE_FROM_ALL
    benchmarks/benchmark_main.cpp
    benchmarks/random_expression_fuzzer.cpp)
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

set(MMCAL_TARGETS mmCal.Core mmCal mmCal.Tests mmCal.Benchmarks)

if(MSVC)
    foreach(target IN LISTS MMCAL_TARGETS)
        target_compile_options(${target} PRIVATE /W4 /utf-8)
    endforeach()
else()
    foreach(target IN LISTS MMCAL_TARGETS)
        target_compile_options(${target} PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            $<$<CONFIG:Release>:-O3>)
    endforeach()

    # Keep floating-point semantics conservative. No -ffast-math is used.
    if(MMCAL_NATIVE_OPTIMIZATION)
        foreach(target IN LISTS MMCAL_TARGETS)
            target_compile_options(${target} PRIVATE $<$<CONFIG:Release>:-march=native>)
        endforeach()
    endif()
endif()

if(MMCAL_ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT MMCAL_IPO_SUPPORTED OUTPUT MMCAL_IPO_ERROR)
    if(MMCAL_IPO_SUPPORTED)
        foreach(target IN LISTS MMCAL_TARGETS)
            set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
        endforeach()
    else()
        message(WARNING "LTO/IPO is not supported by this toolchain: ${MMCAL_IPO_ERROR}")
    endif()
endif()

if(UNIX AND NOT APPLE AND MMCAL_LINUX_OUT_SUFFIX)
    set_target_properties(mmCal PROPERTIES OUTPUT_NAME "mmCal.out")
endif()

enable_testing()
add_test(NAME mmCal.Tests COMMAND mmCal.Tests)
