/*  Filename   :     UnkonwnInterface.h
 *  Creator    :     Clinton Stimpson
 *
 *  Date       :     10 Jan 2002
 *
 *  Owner      :     Clinton Stimpson
 *
 *  Description:     Contains declarations for MBuuid which keeps
 *                   track of different interfaces.
 *                   Also contains the declaration for the base class
 *                   UknownInterface from which all interfaces are
 *                   derived from
 */

#ifndef MOAB_UNKNOWN_INTERFACE_HPP
#define MOAB_UNKNOWN_INTERFACE_HPP

#include <memory.h>

namespace moab {

//!  struct that handles universally unique id's for the Mesh Database

// note: this MBuuid is compliant with the windows GUID.  
// It is possible to do a memcpy() to copy the data from a MBuuid to a GUID
// if we want to support dll registration
struct MBuuid
{
   //! default constructor that initializes to zero
   MBuuid()
   {
      memset( this, 0, sizeof(MBuuid) );
   }
   //! constructor that takes initialization arguments
   MBuuid( unsigned l, unsigned short w1, unsigned short w2, 
         unsigned char b1, unsigned char b2, unsigned char b3, 
         unsigned char b4, unsigned char b5, unsigned char b6, 
         unsigned char b7, unsigned char b8 )
   {
      data1 = l;
      data2 = w1;
      data3 = w2;
      data4[0] = b1;
      data4[1] = b2;
      data4[2] = b3;
      data4[3] = b4;
      data4[4] = b5;
      data4[5] = b6;
      data4[6] = b7;
      data4[7] = b8;
   }
   //! copy constructor
   MBuuid( const MBuuid& mdbuuid )
   {
      memcpy( this, &mdbuuid, sizeof(MBuuid));
   }
   //! sets this uuid equal to another one
   MBuuid &operator=(const MBuuid& orig)
   {
      memcpy( this, &orig, sizeof(MBuuid));
      return *this;
   }
   //! returns whether two uuid's are equal
   bool operator==(const MBuuid& orig) const
   {
      return !memcmp(this, &orig, sizeof(MBuuid));
   }
   //! returns whether two uuid's are not equal
   bool operator!=(const MBuuid& orig) const
   {
      return!(*this == orig);
   }

   //! uuid data storage
   unsigned   data1;
   unsigned short data2;
   unsigned short data3;
   unsigned char  data4[8];
};
  
//! uuid for an unknown interface
//! this can be used to either return a default interface
//! or a NULL interface
static const MBuuid IDD_MBUnknown = MBuuid( 0xf4f6605e, 0x2a7e, 0x4760, 
   0xbb, 0x06, 0xb9, 0xed, 0x27, 0xe9, 0x4a, 0xec );


//! base class for all interface classes
class UnknownInterface
{
public:
   virtual int QueryInterface
      ( const MBuuid&, UnknownInterface** ) = 0;
   virtual ~UnknownInterface() {}
};

} // namespace moab 

#endif  // MOAB_UNKNOWN_INTERFACE_HPP

