cmake_minimum_required(VERSION 3.16)
project(eof_toolbox LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)
include(FetchContent)

find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
  )
  FetchContent_GetProperties(eigen)
  if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
  endif()
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  target_include_directories(Eigen3::Eigen INTERFACE "${eigen_SOURCE_DIR}")
endif()

add_library(eof src/eof.cpp)
add_library(eof::eof ALIAS eof)
target_include_directories(eof
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(eof PUBLIC Eigen3::Eigen)
target_compile_features(eof PUBLIC cxx_std_17)

enable_testing()
add_executable(test_eof tests/test_eof.cpp)
target_link_libraries(test_eof PRIVATE eof::eof)
add_test(NAME eof COMMAND test_eof)

install(TARGETS eof EXPORT eofTargets)
install(DIRECTORY include/eof DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
