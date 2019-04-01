/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */

#ifndef __partitioner_base_hpp__
#define __partitioner_base_hpp__

#include <stdlib.h>
#include <vector>

#include "moab/MOABConfig.h"
#include "moab/Range.hpp"
#include "moab/Types.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#endif
namespace moab {

  class Interface;
}

using namespace moab;

  class PartitionerBase 
  {

  public:
    PartitionerBase( Interface *impl = NULL,
                      const bool use_coords = false);
    
    virtual ~PartitionerBase();

    virtual ErrorCode partition_mesh_and_geometry(const double part_geom_mesh_size,
                                                  const int nparts,
                                                  const char *zmethod,
                                                  const char *other_method,
                                                  double imbal_tol,
                                                  const int part_dim = 3,
                                                  const bool write_as_sets = true,
                                                  const bool write_as_tags = false,
                                                  const int obj_weight = 0,
                                                  const int edge_weight = 0,
                                                  const bool part_surf = false,
                                                  const bool ghost = false,
                                                  const bool spherical_coords = false,
                                                  const bool print_time = false) = 0;

    virtual ErrorCode partition_mesh( const int nparts,
                                      const char *method,
                                      const int part_dim = 3, 
                                      const bool write_as_sets = true,
                                      const bool write_as_tags = false,
                                      const bool partition_tagged_sets = false,
                                      const bool partition_tagged_ents = false,
                                      const char *aggregating_tag = NULL,
                                      const bool print_time = false) = 0;

    virtual ErrorCode write_partition(const int nparts, Range &elems, 
                                const int *assignment,
                                const bool write_as_sets,
                                const bool write_as_tags) = 0;

      // put closure of entities in the part sets too
    virtual ErrorCode include_closure() = 0;

    Range &part_sets() {return partSets;};
    
    const Range &part_sets() const {return partSets;};

  protected:

    Interface *mbImpl;
#ifdef MOAB_HAVE_MPI
    ParallelComm *mbpc;
#endif
    bool write_output;
    bool useCoords;
    bool newComm;

    Range partSets;

  };

inline
PartitionerBase::PartitionerBase(Interface *impl,
                                  const bool use_coords)
    : mbImpl(impl), useCoords(use_coords), newComm(false)
{
#ifdef MOAB_HAVE_MPI
  mbpc = ParallelComm::get_pcomm(mbImpl, 0);
  if (!mbpc) {
    mbpc = new ParallelComm(impl, MPI_COMM_WORLD, 0);
    newComm = true;
  }
#endif
}

inline
PartitionerBase::~PartitionerBase()
{
#ifdef MOAB_HAVE_MPI
  if (newComm)
    delete mbpc;
#endif
  mbImpl = NULL;
}

#endif

