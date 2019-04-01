/** \file   Compiler.hpp
 *  \author Jason Kraftcheck 
 *  \date   2010-12-16
 *
 * Provide pre-processor macros for compiler-specific features.  All
 * defined macros should expand to nothing if not supported by the
 * compiler.
 */

#ifndef moab_COMPILER_HPP
#define moab_COMPILER_HPP

#ifdef IS_BUILDING_MB

/**\defgroup PRIVCOMP Private Compiler-Specifc Pre-Processor Macros */
/*@{*/

/**\def __restrict__
 *\brief Provide functionality similar to C99 \c restrict keyword
 *
 * Tell the compiler that a pointer is not aliased.  This means that
 * programmer guarantees that no other pointer will be used to reference
 * memory that is referenced through the designated pointer unless it
 * is obivous to the compiler in the relevant function.  A typical use
 * for this is to specify that two pointer arguments to a function will
 * never be used to reference overlapping memory.  For example:
 *\code
 *  void* memcpy(void* __restrict__ dest, const void* __restrict__ src, size_t len);
 *\endcode
 * Says that the memory locations indicated by the \c dest and \c src pointers
 * will never be used to reference overlapping memory, including offsets up to
 * \c len.
 *
 * Notifying the compiler about lack of pointer aliasing allows it to make
 * better optimizations.  However, the behavior is undefined (and probably
 * broken in platform-specific ways) if designated pointers are aliased.
 */
#ifdef __cplusplus
# if !defined __GNUC__ || __GNUC__ < 4 || __GNUC_MINOR__ < 5
#   define __restrict__
# endif
#endif

/*@}*/

#endif

/**\defgroup PUBCOMP Public Compiler-Specifc Pre-Processor Macros */
/*@{*/

/**\def PRINT_FORMAT(start)
 *\brief Give a hint to the compiler the function is like \c printf
 *
 * Tell the compiler that the function involves a printf-style format
 * string and varargs list.  This gives the compiler the opportunity
 * to warn if the argument types do not match the format string.
 * This macro should be inluded after the complete function declaration,
 * but before the closing semi-colon.
 *
 *\param START The position of the format string in the argument list, where
 *             the first argument is 1.
 *\NOTE This macro is designed to be used with member functions of C++ classes,
 *      and therefore explicitly accounts for the implicit \c this pointer
 *      in the argument list.  It will not work correctly with static or 
 *      non-member functions.
 *\NOTE This macro assumes that the arguments referenced in the format string
 *      begin immediately after the format string itself.
 */
#ifdef __GNUC__
  #define MB_PRINTF(START) __attribute__((format(printf,(START)+1,(START)+2)))
#else
  #define MB_PRINTF(START)
#endif
 
/**\def MB_DLL_EXPORT
 *\brief Declare a function or class to be visible in shared library.
 */
/**\def MB_DLL_HIDDEN
 *\brief Declare a function or class to be internal to a shared library.
 */
#if defined _MSC_VER || defined __CYGWIN__ || defined __MINGW32__ \
 || defined __MINGW64__ || defined _WIN32
  #if !defined IS_BUILDING_MB || !defined MB_EXPORTS
    #define MB_DLL_EXPORT __dllspec(dllexport)
  #elif !defined MB_WIN_DLL
    #define MB_DLL_EXPORT __dllspec(dllimport)
  #else
    #define MB_DLL_EXPORT
  #endif
  #define MB_DLL_HIDDEN
#elif defined __GNUC__ && __GNUC__ > 3
  #define MB_DLL_EXPORT __attribute__((visibility("default")))
  #define MB_DLL_HIDDEN __attribute__((visibility("hidden")))
#else
  #define MB_DLL_EXPORT
  #define MB_DLL_HIDDEN
#endif

/**\def MB_DEPRECATED
 *\brief Mark function or API as deprecated
 */
#if defined(__GNUC__) && (1000 * __GNUC__ + __GNUC_MINOR__ ) > 3000
#  define MB_DEPRECATED __attribute__((__deprecated__))
#else
#  define MB_DEPRECATED
#endif

/*@}*/

#endif // moab_COMPILER_HPP
