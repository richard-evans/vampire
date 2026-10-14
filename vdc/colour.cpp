//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans, Daniel Meilak 2017-2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations
double scale( double scale_min, double scale_max, double start_min, double start_max, double x );
void rgb2hsl( double& red, double& green, double& blue, double& hue, double& light, double& saturation );
void rgb2hsi( double& red, double& green, double& blue, double& h2, double& intensity, double& saturation );
void hsl2rgb( double& red, double& green, double& blue, double& hue, double& light, double& saturation );
void hsi2rgb( double& red, double& green, double& blue, double& hue, double& intensity, double& saturation );
void interpolate(double xy_angle, double& red, double& green, double& blue );


void rgb( const double& sx, const double& sy, const double& sz, double& red, double& green, double& blue){

   // parameters
   const double pi  = 3.14159265358979323846;

   // variables
   double yz_angle, x_angle, hue, light, saturation;
   double sx2, sy2, sz2; // spin coordinates after change of axes

   // perform change of coordinates
   sx2 = vdc::vector_x[0]*sx + vdc::vector_x[1]*sy + vdc::vector_x[2]*sz;
   sy2 = vdc::vector_y[0]*sx + vdc::vector_y[1]*sy + vdc::vector_y[2]*sz;
   sz2 = vdc::vector_z[0]*sx + vdc::vector_z[1]*sy + vdc::vector_z[2]*sz;

   // in y,z plane, angle = the spin direction in the range [0-2pi]
   yz_angle = std::acos(sz2);

   // to apply colourmap, need value between 0-255
   yz_angle = scale(0.0, 255.0, 0.0, 2.0*pi, yz_angle);

   // find rgb values
   // this is a linear interpolation. In future, may want to implement cubic interp
   interpolate( yz_angle, red, green, blue );

   // by default the x-component is ignored and the system is coloured as 2D
   // if --3D is used, colours are brighter/darker according to the x-component
   // of the magnetisation
   if ( vdc::x_axis_colour ){

      // temp values for hue, saturation and light
      hue = 0.0;
      saturation = 0.0;
      light = 0.0;

      // we have resolved xy plane colours. to adjust bightness (x-axis) we
      // convert to HSL
      rgb2hsl(red, green, blue, hue, light, saturation);

      // adjust lightness to reflect x-axis orientation
      // find angle w.r.t. x-axis [value -pi/2 to pi/2]
      x_angle = std::atan(sx2/std::sqrt(sz2*sz2 + sy2*sy2));
      if ( x_angle >= 0.0 ){
         light = light + (1.0 - light)*(2.0*x_angle/pi);
      }
      else {
         light = light + light*(2.0*x_angle/pi);
      }

      // convert back to rgb to get final colour
      hsl2rgb(red, green, blue, hue, light, saturation);

   }

   // Print out spin colours (for debugging)
   //std::cout << sx << "\t" << sy << "\t" << sz << "\t" << red << "\t" << green << "\t" << blue << std::endl;

   if ( red < 0.0 || red > 1.0 ){
      std::cout << "Error red = " << red << std::endl;
   }
   if ( green < 0.0 || green > 1.0 ){
      std::cout << "Error green = " << green << std::endl;
   }
   if ( blue < 0.0 || blue > 1.0 ){
      std::cout << "Error blue = " << blue << std::endl;
   }
   //std::exit(EXIT_FAILURE);
   return;
}

//------------------------------------------------------------------------------
// Scale a data set from range [start_min - start_max]
// to [scale_min - scale_max] with even bias
//------------------------------------------------------------------------------
double scale( double scale_min, double scale_max, double start_min, double start_max, double x ){

   return ((scale_max - scale_min)*(x - start_min)/(start_max - start_min)) + scale_min;
}

//------------------------------------------------------------------------------
// Convert RGB colour to HSL
//------------------------------------------------------------------------------
void rgb2hsl( double& red, double& green, double& blue, double& hue, double& light, double& saturation ){

   // function specific variables
   double maximum, minimum, chroma;

   // variables used to work out hue
   maximum = std::max(std::max(red,green),blue);
   minimum = std::min(std::min(red,green),blue);
   chroma  = maximum - minimum;

   // hue, depends on maximum value
   if ( (chroma <= 0.000001) && (chroma >= -0.000001) ){
      hue = 0.0;
   }
   else if ( (maximum <= (red + 0.000001)) && (maximum >=  (red - 0.000001)) ){
      hue = std::fmod(std::fmod( (green - blue)/chroma, 6.0 ) + 6.0, 6.0);
   }
   else if ( (maximum <= (green + 0.000001)) && (maximum >=  (green - 0.000001)) ){
      hue = (blue - red)/chroma + 2.0;
   }
   else if ( (maximum <= (blue + 0.000001)) && (maximum >=  (blue - 0.000001)) ){
      hue = (red - green)/chroma + 4.0;
   }
   hue = hue*60; // should be unnecessary as undone later

   // lightness
   light = 0.5*(maximum + minimum);

   // saturation
   if ( (light <= 1.000001) && (light >= (1.0 - 0.000001)) ){
      saturation = 0.0;
   }
   else {
      saturation = chroma/(1.0 - std::abs(2.0*light - 1.0));
   }

   return;
}

//------------------------------------------------------------------------------
// Convert RGB to HSI
//------------------------------------------------------------------------------
void rgb2hsi( double& red, double& green, double& blue, double& h2, double& intensity, double& saturation ){

   const double pi  = 3.14159265358979323846;
   double alpha, beta;//, c2;

   // calculate hue (h2) and chroma (c2)
   alpha = 0.5*( 2.0*red - green - blue );
   beta  = 0.5*std::sqrt(3.0)*( green - blue );

   h2 = std::atan2(beta,alpha);
   h2 = std::fmod(h2 + 2*pi, 2*pi);
   h2 = h2*180.0/pi;

   //c2 = std::sqrt( alpha*alpha + beta*beta );

   // calculate intensity
   intensity = ( red + green + blue )/3.0;

   // calculate saturation
   if ( (intensity <= 0.000001) && (intensity >= -0.000001) ){
      saturation = 0.0;
   }
   else {
      saturation = 1.0 - std::min(std::min(red,green),blue)/intensity;
   }

   return;
}

//------------------------------------------------------------------------------
// Convert HSL colour to RGB
//------------------------------------------------------------------------------
void hsl2rgb( double& red, double& green, double& blue, double& hue, double& light, double& saturation ){

   // function specific variables
   double chroma, temp_x, modifier;

   // first find new chroma
   chroma = (1.0 - std::abs(2.0*light - 1.0))*saturation;

   // hue
   hue    = hue/60;
   temp_x = chroma*(1.0 - std::abs( std::fmod(hue, 2.0) - 1.0));
   if ( (hue < 0.0) || (hue > 6.0) ){
      red   = 0.0;
      green = 0.0;
      blue  = 0.0;
   }
   else if ( (0.0 <= hue) && (hue <= 1.0) ){
      red   = chroma;
      green = temp_x;
      blue  = 0.0;
   }
   else if ( (1.0 <= hue) && (hue <= 2.0) ){
      red   = temp_x;
      green = chroma;
      blue  = 0.0;
   }
   else if ( (2.0 <= hue) && (hue <= 3.0) ){
      red   = 0.0;
      green = chroma;
      blue  = temp_x;
   }
   else if ( (3.0 <= hue) && (hue <= 4.0) ){
      red   = 0.0;
      green = temp_x;
      blue  = chroma;
   }
   else if ( (4.0 <= hue) && (hue <= 5.0) ){
      red   = temp_x;
      green = 0.0;
      blue  = chroma;
   }
   else if ( (5.0 <= hue) && (hue <= 6.0) ){
      red   = chroma;
      green = 0.0;
      blue  = temp_x;
   }

   // finally we adjust the values by a constant
   modifier = light - 0.5*chroma;
   red      = red   + modifier;
   green    = green + modifier;
   blue     = blue  + modifier;

   return;
}

//------------------------------------------------------------------------------
// Convert HSL colour to RGB
//------------------------------------------------------------------------------
void hsi2rgb( double& red, double& green, double& blue, double& hue, double& intensity, double& saturation ){

   double z, chroma, temp_x, modifier;

   hue    = hue/60;
   z      = 1 - std::abs( std::fmod(hue, 2.0) - 1.0 );
   chroma = (intensity*saturation)/(1.0 + z);
   temp_x = chroma*z;

   if ( (hue < 0.0) || (hue > 6.0) ){
      red   = 0.0;
      green = 0.0;
      blue  = 0.0;
   }
   else if ( (0.0 <= hue) && (hue <= 1.0) ){
      red   = chroma;
      green = temp_x;
      blue  = 0.0;
   }
   else if ( (1.0 <= hue) && (hue <= 2.0) ){
      red   = temp_x;
      green = chroma;
      blue  = 0.0;
   }
   else if ( (2.0 <= hue) && (hue <= 3.0) ){
      red   = 0.0;
      green = chroma;
      blue  = temp_x;
   }
   else if ( (3.0 <= hue) && (hue <= 4.0) ){
      red   = 0.0;
      green = temp_x;
      blue  = chroma;
   }
   else if ( (4.0 <= hue) && (hue <= 5.0) ){
      red   = temp_x;
      green = 0.0;
      blue  = chroma;
   }
   else if ( (5.0 <= hue) && (hue <= 6.0) ){
      red   = chroma;
      green = 0.0;
      blue  = temp_x;
   }

   modifier = intensity*(1.0 - saturation)/3.0;
   red      = 3*(red   + modifier);
   green    = 3*(green + modifier);
   blue     = 3*(blue  + modifier);

   return;
}

//------------------------------------------------------------------------------
// Find point between integer colorvals using y=mx+c (linear interpolation)
// x-values [0-255], y-values rgb [0-1]
//------------------------------------------------------------------------------
void interpolate( double angle, double& red, double& green, double& blue ){

   std::vector<double> m(3), c(3), ymin(3), ymax(3);
   int xmin, xmax;

   // values between which we interpolate
   xmin = int(std::floor(angle));
   xmax = int(std::ceil(angle));

   // find colorval associated with min and max
   ymin = vdc::colourmap[xmin];
   ymax = vdc::colourmap[xmax];

   // if different: work out gradients and intersects
   if ( ymin == ymax ) {
      red   = ymin[0];
      green = ymin[1];
      blue  = ymin[2];
   }
   else {
      for ( int i=0; i < 3; i++){
         m[i] = (ymax[i] - ymin[i])/(xmax - xmin);
         c[i] = ymin[i] - m[i]*xmin;
      }

      // input angle into y=mx+c
      red   = m[0]*angle + c[0];
      green = m[1]*angle + c[1];
      blue  = m[2]*angle + c[2];
   }

   return;
}

} // end of namespace vdc
