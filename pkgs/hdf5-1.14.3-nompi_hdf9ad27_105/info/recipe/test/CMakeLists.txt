cmake_minimum_required(VERSION 3.19)
project(TEST_HDF5_CMAKE LANGUAGES C)


set(HDF5_PREFER_PARALLEL TRUE)
set(HDF5_FIND_DEBUG TRUE)
find_package(HDF5 REQUIRED COMPONENTS C)

add_executable(h5_cmprss)
target_sources(
  h5_cmprss
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/h5_cmprss.c"
)
target_link_libraries(h5_cmprss PUBLIC hdf5::hdf5)
