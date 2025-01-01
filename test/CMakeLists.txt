cmake_minimum_required(VERSION 3.22)

project(PicassoEQTest)

enable_testing()

add_executable(${PROJECT_NAME}
  source/PicassoEQ.t.cpp
)

target_include_directories(${PROJECT_NAME}
  PRIVATE
    ${LIB_DIR}/gtest/googletest/include
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  PicassoEQ
  gtest
  gtest_main
)

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME})