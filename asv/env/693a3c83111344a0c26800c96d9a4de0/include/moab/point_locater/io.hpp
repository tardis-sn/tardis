#ifndef GET_FILE_OPTIONS_H
#define GET_FILE_OPTIONS_H

#include "moab/Core.hpp"

#include <iostream>

namespace moab {
namespace point_locator{
namespace io{
	
template< typename String = std::string, 
	  typename String_vector = std::vector< std::string>,
	  typename Char_vector = std::vector< const char *> >
class File_options{
	public:
	File_options(): meshFiles(), interpTag(), 
			gNormTag(), ssNormTag(),
			readOpts(), outFile(),
			writeOpts(), dbgFile(),
			ssTagNames(), ssTagValues(),
			help(true)
			{}
	String_vector meshFiles;
	String interpTag;
	String gNormTag;
	String ssNormTag;
	String readOpts;
	String outFile;
	String writeOpts;
	String dbgFile;                
	Char_vector ssTagNames;
	Char_vector ssTagValues;
	bool help;
};

// Check first character for a '-'.
// Return true if one is found.  False otherwise.
bool check_for_flag(const char *str) {
  if (str[0] == '-')
    return true;
  else
    return false;
}

// New get_file_options() function with added possibilities for mbcoupler_test.
ErrorCode get_file_options(int argc, char **argv, 
                           std::vector<std::string> &meshFiles,
                           std::string &interpTag,
                           std::string &gNormTag,
                           std::string &ssNormTag,
                           std::vector<const char *> &ssTagNames,
                           std::vector<const char *> &ssTagValues,
                           std::string &readOpts,
                           std::string &outFile,
                           std::string &writeOpts,
                           std::string &dbgFile,
                           bool &help)
{
  #ifdef MESHDIR
  std::string TestDir( STRINGIFY(MESHDIR) );
  #else
  std::string TestDir(".");
  #endif

  // Initialize some of the outputs to null values indicating not present
  // in the argument list.
  gNormTag = "";
  ssNormTag = "";
  readOpts = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=3.0.1;CPUTIME";
  outFile = "";
  writeOpts = "PARALLEL=WRITE_PART;CPUTIME";
  dbgFile = "";
  std::string defaultDbgFile = argv[0];  // The executable name will be the default debug output file.

  // These will indicate if we've gotten our required parameters at the end of parsing.
  bool haveMeshes = false;
  bool haveInterpTag = false;

  // Loop over the values in argv pulling out an parsing each one
  int npos = 1;

  if (argc > 1 && argv[1] == std::string("-h")) {
    help = true;
    return MB_SUCCESS;
  }

  while (npos < argc) {
    if (argv[npos] == std::string("-meshes")) {
      // Parse out the mesh filenames
      npos++;
      int numFiles = 2;
      meshFiles.resize(numFiles);
      for (int i = 0; i < numFiles; i++) {
        if ((npos < argc) && (!check_for_flag(argv[npos])))
          meshFiles[i] = argv[npos++];
        else {
          std::cerr << "    ERROR - missing correct number of mesh filenames" << std::endl;
          return MB_FAILURE;
        }
      }

      haveMeshes = true;
    }
    else if (argv[npos] == std::string("-itag")) {
      // Parse out the interpolation tag
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        interpTag = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <interp_tag>" << std::endl;
        return MB_FAILURE;
      }

      haveInterpTag = true;
    }
    else if (argv[npos] == std::string("-gnorm")) {
      // Parse out the global normalization tag
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        gNormTag = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <gnorm_tag>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-ssnorm")) {
      // Parse out the subset normalization tag and selection criteria
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        ssNormTag = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <ssnorm_tag>" << std::endl;
        return MB_FAILURE;
      }

      if ((npos < argc) && (!check_for_flag(argv[npos]))) {
        char* opts = argv[npos++];
        char sep1[1] = {';'};
        char sep2[1] = {'='};
        bool end_vals_seen = false;
        std::vector<char *> tmpTagOpts;

        // first get the options
        for (char* i = strtok(opts, sep1); i; i = strtok(0, sep1)) {
          tmpTagOpts.push_back(i);
        }

        // parse out the name and val or just name.
        for (unsigned int j = 0; j < tmpTagOpts.size(); j++) {
          char* e = strtok(tmpTagOpts[j], sep2);
          ssTagNames.push_back(e);
          e = strtok(0, sep2);
          if (e != NULL) {
            // We have a value
            if (end_vals_seen) {
              // ERROR we should not have a value after none are seen
              std::cerr << "    ERROR - new value seen after end of values in <ssnorm_selection>" << std::endl;
              return MB_FAILURE;
            }
            // Otherwise get the value string from e and convert it to an int
            int *valp = new int;
            *valp = atoi(e);
            ssTagValues.push_back((const char *) valp);
          }
          else {
            // Otherwise there is no '=' so push a null on the list
            end_vals_seen = true;
            ssTagValues.push_back((const char *) 0);
          }
        }
      }
      else {
        std::cerr << "    ERROR - missing <ssnorm_selection>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-ropts")) {
      // Parse out the mesh file read options
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        readOpts = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <roptions>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-outfile")) {
      // Parse out the output file name
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        outFile = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <out_file>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-wopts")) {
      // Parse out the output file write options
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        writeOpts = argv[npos++];
      else {
        std::cerr << "    ERROR - missing <woptions>" << std::endl;
        return MB_FAILURE;
      }
    }
    else if (argv[npos] == std::string("-dbgout")) {
      // Parse out the debug output file name.
      // If no name then use the default.
      npos++;
      if ((npos < argc) && (!check_for_flag(argv[npos])))
        dbgFile = argv[npos++];
      else
        dbgFile = defaultDbgFile;
    }
    else {
      // Unrecognized parameter.  Skip it and move along.
      std::cerr << "    ERROR - Unrecognized parameter:" << argv[npos] << std::endl;
      std::cerr << "            Skipping..." << std::endl;
      npos++;
    }
  }

  if (!haveMeshes) {
    meshFiles.resize(2);
    meshFiles[0] = std::string(TestDir + "/64bricks_1khex.h5m");
    meshFiles[1] = std::string(TestDir + "/64bricks_12ktet.h5m");
    std::cout << "Mesh files not entered; using default files " 
              << meshFiles[0] << " and " << meshFiles[1] << std::endl;
  }
  
  if (!haveInterpTag) {
    interpTag = "vertex_field";
    std::cout << "Interpolation field name not given, using default of " << interpTag << std::endl;
  }

#ifdef MOAB_HAVE_HDF5
  if (1 == argc) {
    std::cout << "No arguments given; using output file dum.h5m." << std::endl;
    outFile = "dum.h5m";
  }
#endif
    
  return MB_SUCCESS;
}

// "generic" get_file_options
template<typename Options>
ErrorCode get_file_options(int argc, char **argv, Options & o){
	return	get_file_options( argc, argv, 
			  o.meshFiles, o.interpTag,
			  o.gNormTag, o.ssNormTag,
			  o.ssTagNames, o.ssTagValues,
			  o.readOpts, o.outFile,
			  o.writeOpts, o.dbgFile,
			  o.help);
}



} //namespace io
} //namespace point_locator
} //namespace moab

#endif // GET_FILE_OPTIONS_H
