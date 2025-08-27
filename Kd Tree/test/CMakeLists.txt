cmake_minimum_required(VERSION 3.11)
project(cs350-test)

# Testing
add_executable(${PROJECT_NAME} Common.hpp Common.cpp main.cpp test-kdtree.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC cs350-lib)

# GTest
enable_testing()
find_package(GTest CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE GTest::gtest)
include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME})
