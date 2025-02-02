cmake_minimum_required(VERSION 3.22)

project(PicassoEQ)

set(LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libs)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "Loading CPM.cmake...")
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cpm.cmake)
message(STATUS "CPM.cmake loaded.")

CPMAddPackage(
    NAME JUCE
    GITHUB_REPOSITORY juce-framework/JUCE
    GIT_TAG 8.0.4
    VERSION 8.0.4
    SOURCE_DIR ${LIB_DIR}/juce
)

CPMAddPackage(
    NAME eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen
    GIT_TAG 3.4
    VERSION 3.4
    SOURCE_DIR ${LIB_DIR}/eigen
)

CPMAddPackage(
    NAME lbfgspp
    GITHUB_REPOSITORY yixuan/LBFGSpp
    VERSION 0.3.0
    SOURCE_DIR ${LIB_DIR}/lbfgspp
)

CPMAddPackage(
    NAME ollamahpp
    GITHUB_REPOSITORY jmont-dev/ollama-hpp
    VERSION 0.9.4
    SOURCE_DIR ${LIB_DIR}/ollamahpp
)

# CPMAddPackage(
#     NAME gtest
#     GITHUB_REPOSITORY google/googletest
#     VERSION 1.15.2
#     SOURCE_DIR ${LIB_DIR}/gtest
#     OPTIONS
#         "INSTALL_GTEST OFF"
#         "gtest_force_shared_crt ON"
# )

enable_testing()

if (MSVC)
    # add_compile_options(/Wall /WX)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_subdirectory(source)
# add_subdirectory(test)