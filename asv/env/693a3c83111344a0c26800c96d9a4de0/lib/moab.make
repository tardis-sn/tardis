# The values below are for an un-installed copy of MOAB used directly
# from its build build directory.  These values will be overridden below
# for installed copies of MOAB.
MOAB_LIBDIR = /feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/src/.libs
MOAB_INCLUDES = -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/src \
                -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/src \
                -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/src/oldinc \
                -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/src/parallel \
                -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/src/parallel \
                -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/src/LocalDiscretization \
                -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/src/RefineMesh

MOAB_INCLUDES += 

MOAB_CPPFLAGS =       -I/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/include -isystem /home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/include  
MOAB_CXXFLAGS =  -DBOOST_MATH_DISABLE_FLOAT128 -m64  -ftree-vectorize -O2 -DNDEBUG 
MOAB_CFLAGS =  -m64  -ftree-vectorize -O2 -DNDEBUG 
MOAB_FFLAGS =   -ftree-vectorize -O2 -fcray-pointer
MOAB_FCFLAGS =   -ftree-vectorize -O2 -fcray-pointer
MOAB_LDFLAGS =          -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib   -Wl,-rpath,/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib   -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib   -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib

MOAB_LIBS_LINK = ${MOAB_LDFLAGS} -L${MOAB_LIBDIR} -lMOAB    -lhdf5  -lz  -lz  -lz -ldl -lm -lm -lm       -lz -ldl -lm -lm -lm
DAGMC_LIBS_LINK = ${MOAB_LDFLAGS} -L${MOAB_LIBDIR} -ldagmc -lMOAB    -lhdf5  -lz  -lz  -lz -ldl -lm -lm -lm       -lz -ldl -lm -lm -lm

MOAB_CXX = g++
MOAB_CC  = gcc
MOAB_FC  = gfortran
MOAB_F77  = gfortran

# Feature list
MOAB_MPI_ENABLED = no
MOAB_FORTRAN_ENABLED = yes
MOAB_HDF5_ENABLED = yes
MOAB_NETCDF_ENABLED = no
MOAB_PNETCDF_ENABLED = no
MOAB_IGEOM_ENABLED = no
MOAB_IMESH_ENABLED = yes
MOAB_IREL_ENABLED = no

# Override MOAB_LIBDIR and MOAB_INCLUDES from above with the correct
# values for the installed MOAB.

MOAB_LIBDIR=/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib
MOAB_INCLUDES=-I/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/include
