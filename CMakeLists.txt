cmake_minimum_required(VERSION 3.0.0)
project(WDP CXX)

add_executable(im im.cpp)
add_executable(moran moran.cpp)

target_compile_features(im PRIVATE cxx_std_14)
target_compile_features(moran PRIVATE cxx_std_14)
