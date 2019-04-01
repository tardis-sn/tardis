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


#ifndef MOAB_READER_WRITER_SET_HPP
#define MOAB_READER_WRITER_SET_HPP

#include <list>
#include <string>
#include "moab/Types.hpp"

namespace moab {

class ReaderIface;
class WriterIface;
class Core;

/**
 *\brief Maintain list of readers and writers.
 *\version 1.00
 *\date 2004-4-23
 *\author Jason Kraftcheck
 */
class ReaderWriterSet
{

  public:
    
    typedef ReaderIface* (*reader_factory_t)( Interface* );
    typedef WriterIface* (*writer_factory_t)( Interface* );
  
    ReaderWriterSet( Core* mdb );
  
    ~ReaderWriterSet();
    
    /**
     * Regiseter a reader and/or writer
     * Either factory function may be NULL, but not both.
     *
     *\param reader_fact  A factory method to create an instance of the reader
     *\param writer_fact  A factory method to create an instance of the reader
     *\param description  A short description of the file format.
     *\param extensions   A null-terminated list of file extensions
     *\param name         File format identifier string.
     */
    ErrorCode register_factory( reader_factory_t reader_fact,
                                  writer_factory_t writer_fact,
                                  const char* description,
                                  const char* const* extensions,
                                  const char* name );
    ErrorCode register_factory( reader_factory_t reader_fact,
                                  writer_factory_t writer_fact,
                                  const char* description,
                                  const char* extension,
                                  const char* name );
  
    /** 
     * Create a reader object for the passed file name 
     * according to the dot-extension of the file name.
     * Caller must delete the object when finished.
     * Returns null if no matching file extension.
     */
    ReaderIface* get_file_extension_reader( const std::string& filename ) const;

    /** 
     * Create a writer object for the passed file name 
     * according to the dot-extension of the file name.
     * Caller must delete the object when finished.
     * Returns null if no matching file extension.
     */
    WriterIface* get_file_extension_writer( const std::string& filename ) const;
    
    /**
     * Create a reader object for the passed file format type.
     * Caller is responsible for deletion of returned object.
     * Returns NULL if no match.
     */
    ReaderIface* get_file_reader( const char* format_name ) const; 
     
    /**
     * Create a writer object for the passed file format type.
     * Caller is responsible for deletion of returned object.
     * Returns NULL if no match.
     */
    WriterIface* get_file_writer( const char* format_name ) const; 
    
    /** 
     * Get the file extension from a file name
     */
    static std::string extension_from_filename( const std::string& filename );
  
    class Handler {
      
      friend class ReaderWriterSet;
      
      public:
      
      Handler( reader_factory_t read_f,
               writer_factory_t write_f,
               const char* name,
               const char* desc, 
               const char* const* ext, 
               int num_ext );
      
      inline const std::string& name() const { return mName; }
      inline const std::string& description() const { return mDescription; }
      inline void get_extensions( std::vector<std::string>& list_out ) const
        { list_out = mExtensions; }
      
      inline bool have_reader() const { return NULL != mReader; }
      inline bool have_writer() const { return NULL != mWriter; }
      
      inline ReaderIface* make_reader( Interface* iface ) const
        { return have_reader() ? mReader(iface) : NULL; }
      
      inline WriterIface* make_writer( Interface* iface ) const
        { return have_writer() ? mWriter(iface) : NULL; }

      bool reads_extension(const char *ext) const;
      bool writes_extension(const char *ext) const;
      
      bool operator==( const char* name ) const;
      
      private:
      
      reader_factory_t mReader;
      writer_factory_t mWriter;
      
      std::string mName, mDescription;
      std::vector<std::string> mExtensions;
    };
    
    typedef std::list<Handler>::const_iterator iterator;
    
    inline iterator begin() const { return handlerList.begin(); }
    
    inline iterator end()   const { return handlerList.end();   }
    
    
    iterator handler_from_extension( const std::string& extension,
                                      bool with_reader = false, 
                                      bool with_writer = false) const;
    
    iterator handler_by_name( const char* name ) const;
    
  private:
  
    Core* mbCore;
  
    std::list<Handler> handlerList;
};

} // namespace moab 

#endif
