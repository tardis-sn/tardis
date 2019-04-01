# Config file for MOAB; use the CMake find_package() function to pull this into
# your own CMakeLists.txt file.
#
# This file defines the following variables:
# MOAB_FOUND        - boolean indicating that MOAB is found
# MOAB_INCLUDE_DIRS - include directories from which to pick up MOAB includes
# MOAB_LIBRARIES    - libraries need to link to MOAB; use this in target_link_libraries for MOAB-dependent targets
# MOAB_CXX, MOAB_CC, MOAB_F77, MOAB_FC - compilers used to compile MOAB
# MOAB_CXXFLAGS, MOAB_CCFLAGS, MOAB_FFLAGS, MOAB_FCFLAGS - compiler flags used to compile MOAB; possibly need to use these in add_definitions or CMAKE_<LANG>_FLAGS_<MODE> 

SET(MOAB_FOUND 1)

# Compilers used by MOAB

SET(MOAB_CXX "g++")
SET(MOAB_CC "gcc")
SET(MOAB_F77 "gfortran")
SET(MOAB_FC "gfortran")


# Compiler flags used by MOAB

SET(MOAB_CXXFLAGS " -DBOOST_MATH_DISABLE_FLOAT128 -m64  -ftree-vectorize -O2 -DNDEBUG @AM_CXXFLAGS@")
SET(MOAB_CFLAGS " -m64  -ftree-vectorize -O2 -DNDEBUG @AM_CFLAGS@")
SET(MOAB_FORTRAN_FLAGS "  -ftree-vectorize -O2 -fcray-pointer @AM_FFLAGS@")

# Library and include defs
set(MOAB_INCLUDE_DIRS, "/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/include")
set(MOAB_LIBRARY_DIRS, "/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib")
SET(MOAB_INCLUDE_DIRS "/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/include")
SET(MOAB_LIBRARIES "         -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib  -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib -lMOAB    -lhdf5  -lz  -lz  -lz -ldl -lm -lm -lm     -lz -ldl -lm -lm -lm")
SET(DAGMC_LIBRARIES "         -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib  -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib -ldagmc -lMOAB    -lhdf5  -lz  -lz  -lz -ldl -lm -lm -lm     -lz -ldl -lm -lm -lm")
