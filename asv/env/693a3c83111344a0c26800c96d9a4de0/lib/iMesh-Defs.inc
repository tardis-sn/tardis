IMESH_CXXFLAGS =  -DBOOST_MATH_DISABLE_FLOAT128 -m64  -ftree-vectorize -O2 -DNDEBUG 
IMESH_CFLAGS =  -m64  -ftree-vectorize -O2 -DNDEBUG
IMESH_LDFLAGS =          -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib  -Wl,-rpath,/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib   -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib   -L/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib
IMESH_FC = gfortran
IMESH_FCFLAGS =   -ftree-vectorize -O2 -fcray-pointer
FC_DEFINE = -D
IMESH_FCDEFS = -DTEMPLATE_SPECIALIZATION -DTEMPLATE_FUNC_SPECIALIZATION

IMESH_CPPFLAGS =       -I/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/include -isystem /home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/include  -I$(IMESH_INCLUDEDIR) $(IMESH_INCLUDEDIR2)

# This definition is deprecated; use IMESH_CPPFLAGS instead
IMESH_INCLUDES = $(IMESH_CPPFLAGS)

IMESH_LIBS = $(IMESH_LDFLAGS) -L$(MOAB_LIBDIR) \
             -L$(IMESH_LIBDIR) -liMesh -lMOAB \
               -lhdf5  -lz  -lz  -lz -ldl -lm -lm -lm  \
                   \
              -lstdc++  -lz -ldl -lm -lm -lm

IMESH_INCLUDEDIR = /feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/itaps/imesh
IMESH_INCLUDEDIR2 = -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/itaps/imesh -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/itaps/imesh/.. -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/itaps/imesh/.. -I/feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/src
IMESH_LIBDIR = /feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/itaps/imesh/.libs
MOAB_LIBDIR = /feedstock_root/build_artefacts/moab_1510674180998/work/moab-4.9.1/src/.libs
IMESH_INCLUDEDIR=/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/include
IMESH_INCLUDEDIR2=
IMESH_LIBDIR=/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib
MOAB_LIBDIR=/home/youssef/Tardis/Tardis/tardis/asv/env/693a3c83111344a0c26800c96d9a4de0/lib
