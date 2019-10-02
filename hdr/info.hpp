#ifndef INFO_H_
#define INFO_H_
//-----------------------------------------------------------------------------
//
// This header is part of the VAMPIRE open source package under the
// Free BSD licence (see licence file for details).
//
// (c) R F L Evans 2018. All rights reserved.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Functions to provide code version and revision information for vampire
//-----------------------------------------------------------------------------
namespace vinfo{

   //-------------------------------------------------------------------------
   // Shared variables
   //-------------------------------------------------------------------------

   //-------------------------------------------------------------------------
   // Function declarations
   //-------------------------------------------------------------------------
   std::string githash(); // return code revision hash
   std::string version(); // return code revision hash

}

#endif /* INFO_H_ */
