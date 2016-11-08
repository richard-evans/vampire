/// Output to povray

#include "cfgConverter.hpp"

void rgb( double ireal, double &red, double &green, double &blue){

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
   
}