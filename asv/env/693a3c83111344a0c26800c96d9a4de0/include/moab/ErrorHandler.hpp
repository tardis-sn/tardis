#ifndef MOAB_ERROR_HANDLER_HPP
#define MOAB_ERROR_HANDLER_HPP

#ifdef _WIN32
#define __func__ __FUNCTION__
#endif

#include "moab/Types.hpp"

#include <sstream>
#include <string.h>

namespace moab {

//! ErrorType - passed to the error handling routines indicating whether this is a new error (globally
//! fatal or per-processor relevant) to be created, or an existing one to be handled
enum ErrorType {MB_ERROR_TYPE_NEW_GLOBAL = 0, MB_ERROR_TYPE_NEW_LOCAL = 1, MB_ERROR_TYPE_EXISTING = 2};

//! Initialize MOAB error handler (e.g. create a utility object for printing error output)
void MBErrorHandler_Init();

//! Finalize MOAB error handler (e.g. delete the utility object for printing error output)
void MBErrorHandler_Finalize();

//! Indicates whether MBErrorHandler_Init has been called
bool MBErrorHandler_Initialized();

//! Get information about the last error
void MBErrorHandler_GetLastError(std::string& error);

//! Routine that is called to create a new error or handle an existing one
ErrorCode MBError(int line, const char* func, const char* file, const char* dir,
                  ErrorCode err_code, const char* err_msg, ErrorType err_type);


} // namespace moab

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define MBSTRINGIFY_(X) #X
#define MBSTRINGIFY(X) MBSTRINGIFY_(X)

#ifdef LOCDIR
#define __MBSDIR__ MBSTRINGIFY(LOCDIR)
#else
#define __MBSDIR__ ""
#endif

//! Set a new error with the given error message (a string or a stream) and return the given error code
//! Used in functions which return ErrorCode
#define MB_SET_ERR(err_code, err_msg) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg; \
    return moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, err_code, err_ostr.str().c_str(), moab::MB_ERROR_TYPE_NEW_LOCAL); \
  } while (false)

//! Set a new error with the given error message (a string or a stream) and return
//! Used in functions which return void types (or have no return types at all, e.g. constructors)
#define MB_SET_ERR_RET(err_msg) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg; \
    moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, moab::MB_FAILURE, err_ostr.str().c_str(), moab::MB_ERROR_TYPE_NEW_LOCAL); \
    return; \
  } while (false)

//! Set a new error with the given error message (a string or a stream) and return the given value
//! Used in functions which return any data type
#define MB_SET_ERR_RET_VAL(err_msg, ret_val) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg; \
    moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, moab::MB_FAILURE, err_ostr.str().c_str(), moab::MB_ERROR_TYPE_NEW_LOCAL); \
    return ret_val; \
  } while (false)

//! Set a new error with the given error message (a string or a stream) and continue
//! Used in functions which return any data type
#define MB_SET_ERR_CONT(err_msg) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg; \
    moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, moab::MB_FAILURE, err_ostr.str().c_str(), moab::MB_ERROR_TYPE_NEW_LOCAL); \
  } while (false)

//! Similar to MB_SET_ERR except that the error is considered globally fatal
#define MB_SET_GLB_ERR(err_code, err_msg) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg; \
    return moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, err_code, err_ostr.str().c_str(), moab::MB_ERROR_TYPE_NEW_GLOBAL); \
  } while (false)

//! Similar to MB_SET_ERR_RET except that the error is considered globally fatal
#define MB_SET_GLB_ERR_RET(err_msg) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg; \
    moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, moab::MB_FAILURE, err_ostr.str().c_str(), moab::MB_ERROR_TYPE_NEW_GLOBAL); \
    return; \
  } while (false)

//! Similar to MB_SET_ERR_RET_VAL except that the error is considered globally fatal
#define MB_SET_GLB_ERR_RET_VAL(err_msg, ret_val) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg; \
    moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, moab::MB_FAILURE, err_ostr.str().c_str(), moab::MB_ERROR_TYPE_NEW_GLOBAL); \
    return ret_val; \
  } while (false)

//! Similar to MB_SET_ERR_CONT except that the error is considered globally fatal
#define MB_SET_GLB_ERR_CONT(err_msg) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg; \
    moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, moab::MB_FAILURE, err_ostr.str().c_str(), moab::MB_ERROR_TYPE_NEW_GLOBAL); \
  } while (false)

//! Check error code, if not MB_SUCCESS, call the error handler and return the given error code
//! Used in functions which return ErrorCode
#define MB_CHK_ERR(err_code) \
  do { \
    if (moab::MB_SUCCESS != err_code) \
      return moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, err_code, "", moab::MB_ERROR_TYPE_EXISTING); \
  } while (false)

//! Check error code, if not MB_SUCCESS, call the error handler and return
//! Used in functions which return void types (or have no return types at all, e.g. constructors)
#define MB_CHK_ERR_RET(err_code) \
  do { \
    if (moab::MB_SUCCESS != err_code) { \
      moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, err_code, "", moab::MB_ERROR_TYPE_EXISTING); \
      return; \
    } \
  } while (false)

//! Check error code, if not MB_SUCCESS, call the error handler and return the given value
//! Used in functions which return any data type
#define MB_CHK_ERR_RET_VAL(err_code, ret_val) \
  do { \
    if (moab::MB_SUCCESS != err_code) { \
      moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, err_code, "", moab::MB_ERROR_TYPE_EXISTING); \
      return ret_val; \
    } \
  } while (false)

//! Check error code, if not MB_SUCCESS, call the error handler and continue
//! Used in functions which return any data type
#define MB_CHK_ERR_CONT(err_code) \
  do { \
    if (moab::MB_SUCCESS != err_code) { \
      moab::MBError(__LINE__, __func__, __FILENAME__, __MBSDIR__, err_code, "", moab::MB_ERROR_TYPE_EXISTING); \
    } \
  } while (false)

//! Check error code, if not MB_SUCCESS, set a new error with the given error message and return the given error code
//! Used in functions which return ErrorCode
#define MB_CHK_SET_ERR(err_code, err_msg) \
  do { \
    if (moab::MB_SUCCESS != err_code) \
      MB_SET_ERR(err_code, err_msg); \
  } while (false)

//! Check error code, if not MB_SUCCESS, set a new error with the given error message and return
//! Used in functions which return void types (or have no return types at all, e.g. constructors)
#define MB_CHK_SET_ERR_RET(err_code, err_msg) \
  do { \
    if (moab::MB_SUCCESS != err_code) \
      MB_SET_ERR_RET(err_msg); \
  } while (false)

//! Check error code, if not MB_SUCCESS, set a new error with the given error message and return the given value
//! Used in functions which return any data type
#define MB_CHK_SET_ERR_RET_VAL(err_code, err_msg, ret_val) \
  do { \
    if (moab::MB_SUCCESS != err_code) \
      MB_SET_ERR_RET_VAL(err_msg, ret_val); \
  } while (false)

//! Check error code, if not MB_SUCCESS, set a new error with the given error message and continue
//! Used in functions which return any data type
#define MB_CHK_SET_ERR_CONT(err_code, err_msg) \
  do { \
    if (moab::MB_SUCCESS != err_code) \
      MB_SET_ERR_CONT(err_msg); \
  } while (false)


#endif
