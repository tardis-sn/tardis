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


#ifndef MOAB_WRITER_IFACE_HPP
#define MOAB_WRITER_IFACE_HPP

#include <vector>
#include <string>
#include "moab/Types.hpp"

namespace moab {

class FileOptions;

/**
 *\brief Interface for mesh writer implementations.
 *\version 1.00
 *\date 2004-4-23
 *\author Jason Kraftcheck
 */
class WriterIface
{
  public:
  
    virtual ~WriterIface() {}
    
    /**
     *\brief Export mesh to a file.
     *
     * Method all writers must provide to export a mesh.
     *
     *\param file_name      The name of the file to create.
     *\param overwrite      If false, reader should fail if the file already
     *                      exists.  
     *\param meshset_list   A list of meshsets to export, or NULL if the
     *                      whole mesh is to be exported.
     *\param num_sets       The length of <code>meshset_list</code> or zero
     *                      if the whole mesh is to be exported.
     *\param qa_records     File history metadata
     *\param tag_list       Array of handles for tags to write.  If null,
     *                      write all tags.  If non-NULL but num_tags is
     *                      zero, write no tags.
     *\param requseted_output_dimension  The geometric dimension of the
     *                      output mesh (coord values per vertex.)  If
     *                      zero, the dimension of the mesh as returned
     *                      from Interface should be used.
     *\author Jason Kraftcheck
     */
    virtual ErrorCode write_file( const char* file_name,
                                    const bool overwrite,
                                    const FileOptions& opts,
                                    const EntityHandle* meshset_list,
                                    const int num_sets,
                                    const std::vector<std::string>& qa_records,
                                    const Tag* tag_list = NULL,
                                    int num_tags = 0,
                                    int requested_output_dimension = 3 ) = 0;
};

} // namespace moab 

#endif

    
