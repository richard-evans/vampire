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

void coulourpallet(){
   
         const vec3 pallete[11] = vec3[11](
         vec3(0,				0,				0.564705882),
         vec3(0,				0.058823529,	1),
         vec3(0,				0.564705882,	1),
         vec3(0.058823529,	1,				0.933333333),
         vec3(0.564705882,	1,				0.439215686),
         vec3(1,				0.933333333,	0),
         vec3(1,				0.439215686,	0),
         vec3(0.933333333,	0,				0),
         vec3(0.498039216,	0,				0),
         vec3(0.498039216,	0,				0),
         vec3(0.498039216,	0,				0));

      const vec3 pallete2[7] = vec3[7](
         vec3(1,				0,				0),
         vec3(1,				1,				0),
         vec3(0,				1,				0),
         vec3(0,				1,				1),
         vec3(0,				0,				1),
         vec3(1,				0,				1),
         vec3(1,				0,				0));


      float red;
      float green;
      float blue;

      int index = int(floor((spin.z*0.5+0.5)*10));

      red   = pallete[index].r + ((spin.z*0.5+0.5)*10-index)*(pallete[index+1].r-pallete[index].r);
      green = pallete[index].g + ((spin.z*0.5+0.5)*10-index)*(pallete[index+1].g-pallete[index].g);
      blue  = pallete[index].b + ((spin.z*0.5+0.5)*10-index)*(pallete[index+1].b-pallete[index].b);


      float theta = (atan(spin.y,spin.x)+3.1415926)/(3.1415926*2);

      //index = int(floor((theta)*7));

      //red   = pallete2[index].r + (theta*7-index)*(pallete2[index+1].r-pallete2[index].r);
      //green = pallete2[index].g + (theta*7-index)*(pallete2[index+1].g-pallete2[index].g);
      //blue  = pallete2[index].b + (theta*7-index)*(pallete2[index+1].b-pallete2[index].b);
}