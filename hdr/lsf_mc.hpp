//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) David R Papp 2024. All rights reserved.
//
//------------------------------------------------------------------------------
//

// Program headers
#ifndef LSF_MC_H_
#define LSF_MC_H_

namespace montecarlo{

   // Declare spin vectors
   extern std::vector<double> mod_S;
   extern bool mc_set;
   extern void mcinit();

}
#endif /*LSF_MC_H_*/
