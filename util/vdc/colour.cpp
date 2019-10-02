//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// program header
#include "vdc.hpp"

namespace vdc{

   void rgb( const double& ireal, double &red, double &green, double &blue){

      if(ireal>0.8){
         red = 0.0;
         green = 0.0;
         blue = 1.0;
      }
      else if(ireal>=0.0){
         red = 1.0-ireal*1.2;
         green = 1.0-ireal*1.2;
         blue = 1.0;
      }
      else if(ireal>=-0.8){
         red = 1.0;
         green = 1.0+ireal*1.2;
         blue = 1.0+ireal*1.2;
      }
      else if(ireal<-0.8){
         red = 1.0;
         green = 0.0;
         blue = 0.0;
      }
      else{
         red = 1.0;
         green = 1.0;
         blue = 1.0;
      }

      if(blue<0.0) blue=0.0;
      if(red<0.0) red=0.0;
      if(green<0.0) green=0.0;

      return;

   }

} // end of namespace vdc
