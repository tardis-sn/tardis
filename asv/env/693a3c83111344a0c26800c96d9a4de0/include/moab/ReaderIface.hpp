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


#ifndef MOAB_READER_IFACE_HPP
#define MOAB_READER_IFACE_HPP

#include "moab/Types.hpp"

#include <vector>

namespace moab {

class FileOptions;

/**
 *\brief Interface for mesh reader implementations.
 *\version 1.00
 *\date 2004-4-23
 *\author Jason Kraftcheck
 */
class ReaderIface
{
  public:
  
    virtual ~ReaderIface() {}
    
      /** Struct used to specify subset of file to read */
    struct IDTag {
      const char* tag_name;  //!< Name of tag containing integer IDs
      const int* tag_values; //!< Array of integer ID values
      int num_tag_values;    //!< Length of tag_values array
    };

    struct SubsetList {
       /** An array of tag name and value sets specifying
        *  the subset of the file to read.  If multiple
        *  tags are specified, the sets that match all
        *  tags (intersection) should be read.
        */
      IDTag* tag_list;
      int tag_list_length;   //!< Length of tag_list array
      int num_parts;         //!< If non-zero, load 1/num_parts of the matching sets
      int part_number;       //!< If num_parts is non-zero, load part_number-th fraction of the sets
    };

    
    /**
     *\brief Load mesh from a file.
     *
     * Method all readers must provide to import a mesh.
     *
     *\param file_name      The file to read.
     *\param file_set       Optional pointer to entity set representing
     *                      file.  If this is not NULL, reader may optionally
     *                      tag the pointed-to set with format-specific
     *                      meta-data.
     *\param subset_list    An optional struct pointer specifying the tags identifying
     *                      entity sets to be read.
     *\param file_id_tag    If specified, reader should store for each entity
     *                      it reads, a unique integer ID for this tag.
     *\author Jason Kraftcheck
     */
    virtual ErrorCode load_file( const char* file_name,
                                 const EntityHandle* file_set,
                                 const FileOptions& opts,
                                 const SubsetList* subset_list = 0,
                                 const Tag* file_id_tag = 0 ) = 0;


    /**
     *\brief Read tag values from a file.
     *
     * Read the list if all integer tag values from the file for
     * a tag that is a single integer value per entity.
     *
     *\param file_name      The file to read.
     *\param tag_name       The tag for which to read values
     *\param tag_values_out Output: The list of tag values.
     *\param subset_list    An array of tag name and value sets specifying
     *                      the subset of the file to read.  If multiple
     *                      tags are specified, the sets that match all
     *                      tags (intersection) should be read.
     *\param subset_list_length The length of the 'subset_list' array.
     */
    virtual ErrorCode read_tag_values( const char* file_name,
                                         const char* tag_name,
                                         const FileOptions& opts,
                                         std::vector<int>& tag_values_out,
                                         const SubsetList* subset_list = 0 ) = 0;
};

} // namespace moab 

#endif

    
