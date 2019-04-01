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


/*  
 *
 *  File:      Error.hpp
 *
 *  Purpose:   To keep track of detail information about errors that occur
 *             in MB.
 *
 *  Date:      9-26-2002
 *
 *  Author:    Clinton Stimpson
 *
 */



#ifndef MOAB_ERROR_HPP
#define MOAB_ERROR_HPP

#ifndef IS_BUILDING_MB
#error "Error.hpp isn't supposed to be included into an application"
#endif

#include <string>
#include <stdarg.h>
#include <stdio.h>

#include "moab/Types.hpp"
#include "moab/Compiler.hpp"

#ifdef WIN32
#define VSNPRINTF _vsnprintf
#else
#define VSNPRINTF vsnprintf
#endif

namespace moab {

class Error
{
  //! string to hold the last error that occurred in MB
  std::string mLastError;

public:

  Error() {}
  ~Error(){}

  ErrorCode set_last_error(const std::string& error) 
  { 
    mLastError = error; 
    return MB_SUCCESS; 
  }

  inline ErrorCode set_last_error(const char* fmt, ...) MB_PRINTF(1);
  
  ErrorCode set_last_error( const char* fmt, va_list args )
  {
    char text[1024];
    VSNPRINTF( text, sizeof(text), fmt, args );
    mLastError = text;
    return MB_SUCCESS;
  }

  ErrorCode get_last_error(std::string& error) const
  { 
    error = mLastError; 
    return MB_SUCCESS;
  }

};

inline ErrorCode Error::set_last_error(const char* fmt, ...)
{
  ErrorCode result = MB_FAILURE;
  if (fmt)
  {
    va_list args;
    va_start( args, fmt );
    result = set_last_error( fmt, args );
    va_end( args );
  }
  return result;
}

} // namespace moab 

#endif


