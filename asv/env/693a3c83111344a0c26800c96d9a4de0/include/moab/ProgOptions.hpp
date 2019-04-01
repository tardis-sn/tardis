#ifndef MOAB_PROGRAM_OPTIONS_H
#define MOAB_PROGRAM_OPTIONS_H

#include <vector>
#include <map>
#include <string>
#include <iostream>


class ProgOpt;


/** A simple command-line option parser and help utility
 *
 * Utility class to specify a program's command-line options arguments, produce a help message
 * explaining how they work, and parse user's command line input (producing useful errors messages
 * if any problems arise).  Loosely (okay, very loosely) inspired by boost program_options.
 * 
 * Options are specified by a comma-separated namestring.  An option named "foo,f" can be specified
 * three ways on the command line: "-f val", "--foo val", or "--foo=val".  The types of options
 * and arguments are specified by function templates.  Valid template values for positional argument
 * and options are int, double, and std::string.  void may also be used in options, and it indicates
 * a command line option that does not take an argument.  
 *
 * Example usage:
 * ProgOptions po( "Example usage of ProgOptions" );
 * po.addOpt<void>( "verbose,v", "Turn on verbose messages" );
 * po.addOpt<std::string> ("foo", "Specify the foo string" );
 * int x = 0;
 * po.addOpt<int>( ",x", "Specify the x number", &x ); // x will be automatically set when options parsed
 * po.parseCommandLine( argc, argv );
 * bool verbose = po.numOptSet("verbose") > 0;
 * std::string foo;
 * if( !po.getOpt( "foo", &foo ) )
 *    foo = "default"; 
 * ...
 * 
 * See the file dagmc_preproc.cpp in the dagmc directory for a real-world example.
 */
class ProgOptions{

public:

  /** 
   * Flags for addOpt and addRequiredArg functions; may be combined with bitwise arithmetic
   * (though not all combinations make sense!)
   **/


  /// Set for a flag that, when detected, prints help text and halts program.
  /// Constructor creates such a flag by default, so the user shouldn't need to use this directly.
  static const int help_flag = 1<<0;

  /// Flag indicating that an option should be given a "cancel" flag.
  /// This creates, for option --foo, an additional option --no-foo that
  /// clears all previously read instances of the foo option
  static const int add_cancel_opt = 1<<1;

  /// When applied to a flag argument (one with template type void), indicate that the 
  /// value 'false' should be stored into the pointer that was given at option creation time.
  /// This overrides the default behavior, which is to store the value 'true'.
  static const int store_false = 1<<2 ;

  /// Specify a numerical flag where any positive integer is an acceptable
  /// value.  E.g. --dimension=3 is equivalent to -3.  Only values in the
  /// range [0,9] are accepted and the flag type must be integer.
  static const int int_flag = 1<<3 ;

  /** Substitue any occurance of the '%' symbol in a string with
   *  the the MPI rank of this process in MPI_COMM_WORLD.  This
   *  option has no effect if not compiled with MPI.  This flag
   *  has no effect for non-string options.
   */
  static const int rank_subst = 1<<4;

  /// Set for a flag that, when detected, will call printVersion() and halt the program.
  static const int version_flag = 1<<5;
  
  ///unimplemented flag for required arguments that may be given multiple times
  //const static int accept_multiple;


  /**
   * @param helptext A brief summary of the program's function, to be printed
   *        when the help flag is detected
   */
  ProgOptions( const std::string& helptext = "",
               const std::string& briefdesc = "" );
  ~ProgOptions();

  /** Specify the program version 
   *
   * Set the program version to a given string.  This will be printed when printVersion()
   * is called.  
   * @param version_string The version string
   * @param addflag If true, a default '--version' option will be added.  If false, 
   *        the version will be set, but no option will be added to the parser.
   */ 
  void setVersion( const std::string& version_string, bool addFlag = true );

  /** Specify a new command-line option
   *
   * Instruct the parser to accept a new command-line argument, as well as specifying
   * how the argument should be handled.  The template parameter indicates the type of 
   * command-line option being specified: acceptable types are void (indicating a flag
   * without an argument), int, double, and std::string.  
   * 
   * @param namestring The command-line options name(s).  Format is longname,shortname.  
   *        If the comma is omitted, or appears only at the end, this option will have
   *        no shortname; if the comma is the first letter of the namestring, the option 
   *        has no longname.  
   * @param helpstring The help information displayed for the option when the program is
   *        invoked with --help
   * @param value A pointer to memory in which to store the parsed value for this option.
   *        If NULL, then the value of the option must be queried using the getOpt function.
   *        If the template parameter is void and value is non-NULL, treat value as a bool* 
   *        and store 'true' into it when the flag is encountered.  (See also store_false, above)
   * @param flags Option behavior flags, which should come from static vars in the ProgOptions
   *        class
   */
  template <typename T>
  void addOpt( const std::string& namestring, const std::string& helpstring,
	       T* value, int flags = 0 );

  /** Specify a new command-line option
   * 
   * This funtion is identical to the 4-arg version, but omits the value parameter, which 
   * is assumed to be NULL
   */
  template <typename T>
  void addOpt( const std::string& namestring, const std::string& helpstring, int flags = 0 ){
    addOpt<T>( namestring, helpstring, NULL, flags );
  }

  /** Add a new line of help text to the option help printout
   *
   * Add a line of text to the option-related help.  Called between calls to addOpt(), 
   * this function can be used to divide the option list into groups of related options
   * to make the help text more convenient.
   */
  void addOptionHelpHeading( const std::string& );


  /** Add required positional argument
   *
   * Add a new required positional argument.  The order in which arguments are specified
   * is the order in which they will be expected on the command line.
   * The template parameter may be int, double, or std::string (but not void)
   * @param helpname The name to give the argument in the help text
   * @param helpstring The help text for the argument
   * @param value Pointer to where parsed value from command line should be stored.  
   *        If NULL, the value must be queried using getReqArg()
   */
  template <typename T>
  void addRequiredArg( const std::string& helpname, const std::string& helpstring, T* value = NULL, int flags = 0 );

  /** Add optional positional arguments
   *
   * Specify location in ordered argument list at which optional arguments
   * may occur.  Optional arguments are allowed at only one location
   * it argument list (this function may not be called more than once.). 
   * The template parameter may be int, double, or std::string (but not void)
   * @param count The maximum number of optional arguments.  Specify zero for unlimited.
   * @param helpname The name to give the argument in the help text
   * @param helpstring The help text for the arguments
   */
  template <typename T>
  void addOptionalArgs( unsigned max_count, const std::string& helpname, const std::string& helpstring, int flags = 0 );

  /** 
   * Print the full help to the given stream
   */
  void printHelp( std::ostream& str = std::cout );

  /**
   * Print only the usage message to the given stream
   */
  void printUsage( std::ostream& str = std::cout );
  
  /**
   * Print the version string to the given stream
   */
  void printVersion( std::ostream& str = std::cout );

  /**
   * Parse command-line inputs as given to main()
   */
  void parseCommandLine( int argc, char* argv[] );

  /**
   *
   * Get the value of the named option.
   * @param namestring The name string given when the option was created.  This need not be 
   *        idential to the created name; only the longname, or the shortname (with comma prefix),
   *        will also work.
   * @param value Pointer to location to store option argument, if any is found
   * @return True if the option was set and its argument was stored into value; false otherwise.
   */
  template <typename T>
  bool getOpt( const std::string& namestring, T* value );

  /**
   * Get a list of values for the named option-- one value for each time it was 
   * given on the command line.
   *
   * This function cannot be called with void as the template parameter; 
   * compilers will reject vector<void> as a type.  This means it cannot be
   * called for flag-type options.  To count the number of times a given flag
   * was specified, use numOptSet()
   * @param namestring See similar argument to getOpt()
   * @param values Reference to list to store values into.  Will have as many entries
   *        as there were instances of this option on the command line
   */
  template <typename T>
  void getOptAllArgs( const std::string& namestring, std::vector<T>& values );

  
  /**
   * @param namestring See similar argument to getOpt()
   * @return The number of times the named option appeared on the command line. 
   */
  int numOptSet( const std::string& namestring );

  /**
   * Retrieve the value of a required command-line argument by name
   * @param namestring The helpname that was given to addRequiredArg when the
   *        desired argument was created
   */
  template <typename T>
  T getReqArg( const std::string& namestring );
  
  /**
   * Append the values of any required or optional arguments
   * @param namestring The helpname that was given to addRequiredArg or
   *                   addOptionalArgs.
   */
  template <typename T>
  void getArgs( const std::string& namestring, std::vector<T>& values );

  /** 
   * Prints an error message to std::cerr, along with a brief usage message, 
   * then halts the program.  Used throughout ProgramOptions implementation. 
   * Users may call this directly if they detect an incorrect usage of program
   * options that the ProgramOptions wasn't able to detect itself.
   * @param message The error message to print before program halt.
   */
  void error( const std::string& message );

  /**
   * Write help data formatted for use as a unix man page.
   */
  void write_man_page( std::ostream& to_this_stream );
  
protected:

  std::string get_option_usage_prefix( const  ProgOpt& option );

  void get_namestrings( const std::string& input, std::string* l, std::string* s );

  ProgOpt* lookup( const std::map<std::string, ProgOpt* >&, const std::string& );
  ProgOpt* lookup_option( const std::string& );

  bool evaluate( const ProgOpt& opt, void* target, const std::string& option, unsigned* arg_idx = NULL);
  bool process_option( ProgOpt* opt, std::string arg, const char* value = 0 );
  
  std::map< std::string, ProgOpt* > long_names;
  std::map< std::string, ProgOpt* > short_names;
  std::map< std::string, ProgOpt* > required_args;

  typedef std::pair<ProgOpt*, std::string> help_line; 
  std::vector< help_line > option_help_strings;
  std::vector< help_line > arg_help_strings;
  std::vector< std::string > main_help;
  std::string brief_help;
  
  bool expect_optional_args;
  unsigned optional_args_position, max_optional_args;
  
  std::string progname;
  std::string progversion;
  
    // if an option was specified with the int_flag, this
    // will contain the long name of the option
  std::string number_option_name;
  
};

#endif /* MOAB_PROGRAM_OPTIONS_H */
