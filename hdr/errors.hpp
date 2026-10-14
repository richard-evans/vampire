//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#include <string>

#ifndef ERRORS_H_
#define ERRORS_H_
//==========================================================
// Namespace errors
//==========================================================
namespace err
{
   extern bool check;
   extern void vexit();
   extern void zexit(std::string);
   extern void v_parallel_all_exit(std::string message = "");

   // print message (red) to screen and log, then exit via vexit
   extern void fatal(std::string const message);
   // print message (yellow) to screen on rank 0 and to the log on every rank
   extern void warning(std::string const message);
}

#endif /*ERRORS_H_*/
