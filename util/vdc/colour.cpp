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
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations
void colourwheel( std::vector<std::vector<double>>& colourmap );
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
   double xy_angle, z_angle, hue, light, saturation;
   double sx2, sy2, sz2; //spin coordinates after change of axes

   // perform change of coordinates
   sx2 = vdc::vector_x[0]*sx + vdc::vector_x[1]*sy + vdc::vector_x[2]*sz;
   sy2 = vdc::vector_y[0]*sx + vdc::vector_y[1]*sy + vdc::vector_y[2]*sz;
   sz2 = vdc::vector_z[0]*sx + vdc::vector_z[1]*sy + vdc::vector_z[2]*sz;

   // in x,y plane, angle = the spin direction
   xy_angle = std::atan2(sy2,sx2);
   xy_angle = std::fmod(xy_angle + 2*pi, 2*pi);  // range [0-2pi]

   // to apply colourmap, need value between 0-255
   xy_angle = scale(0.0, 255.0, 0.0, 2.0*pi, xy_angle);

   // find rgb values
   interpolate( xy_angle, red, green, blue );

   // temp values for hue, saturation and light
   hue = 0.0;
   saturation = 0.0;
   light = 0.0;

   // we have resolved xy plane colours. to adjust bightness (x-axis) we
   // convert to HSL
   rgb2hsl(red, green, blue, hue, light, saturation);

   // adjust lightness to reflect z-axis orientation
   // find angle w.r.t. z-axis [value -pi/2 to pi/2]
   z_angle = std::atan(sz2/std::sqrt(sx2*sx2 + sy2*sy2));
   if ( z_angle >= 0.0 ){
      light = light + (1.0 - light)*(2.0*z_angle/pi);
   }
   else {
      light = light + light*(2.0*z_angle/pi);
   }

   // convert back to rgb to get final colour
   hsl2rgb(red, green, blue, hue, light, saturation);

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
   double alpha, beta, c2;

   // calculate hue (h2) and chroma (c2)
   alpha = 0.5*( 2.0*red - green - blue );
   beta  = 0.5*std::sqrt(3.0)*( green - blue );

   h2 = std::atan2(beta,alpha);
   h2 = std::fmod(h2 + 2*pi, 2*pi);
   h2 = h2*180.0/pi;

   c2 = std::sqrt( alpha*alpha + beta*beta );

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
// Find point between integer colorvals using y=mx+c
// x-values [0-255], y-values rgb [0-1]
//------------------------------------------------------------------------------
void interpolate(double xy_angle, double& red, double& green, double& blue ){

   std::vector<std::vector<double>> colourmap(256, std::vector<double>(3));
   std::vector<double> m(3), c(3), ymin(3), ymax(3);
   int xmin, xmax;

   // values between which we interpolate
   xmin = int(std::floor(xy_angle));
   xmax = int(std::ceil(xy_angle));
   // initialise colourwheel
   colourwheel( colourmap );

   // find colorval associated with min and max
   ymin = colourmap[xmin];
   ymax = colourmap[xmax];

   // work out gradients and intersects
   for ( int i=0; i < 3; i++){
      m[i] = (ymax[i] - ymin[i])/(xmax - xmin);
      c[i] = ymin[i] - m[i]*xmin;
   }
   // input angle into y=mx+c
   red   = m[0]*xy_angle + c[0];
   green = m[1]*xy_angle + c[1];
   blue  = m[2]*xy_angle + c[2];

   return;
}

//------------------------------------------------------------------------------
// initialise colourwheel vector
//------------------------------------------------------------------------------
void colourwheel( std::vector<std::vector<double>>& colourmap ){
   colourmap[0] = {0.937688,0.333524,0.948091};
   colourmap[1] = {0.943832,0.342828,0.942393};
   colourmap[2] = {0.949390,0.352750,0.936128};
   colourmap[3] = {0.954392,0.363232,0.929311};
   colourmap[4] = {0.958860,0.374216,0.921979};
   colourmap[5] = {0.962829,0.385575,0.914155};
   colourmap[6] = {0.966340,0.397272,0.905883};
   colourmap[7] = {0.969436,0.409209,0.897208};
   colourmap[8] = {0.972159,0.421299,0.888172};
   colourmap[9] = {0.974539,0.433497,0.878829};
   colourmap[10] = {0.976626,0.445734,0.869226};
   colourmap[11] = {0.978446,0.457972,0.859397};
   colourmap[12] = {0.980040,0.470168,0.849387};
   colourmap[13] = {0.981421,0.482281,0.839222};
   colourmap[14] = {0.982617,0.494300,0.828938};
   colourmap[15] = {0.983645,0.506196,0.818543};
   colourmap[16] = {0.984527,0.517989,0.808076};
   colourmap[17] = {0.985266,0.529642,0.797533};
   colourmap[18] = {0.985878,0.541157,0.786940};
   colourmap[19] = {0.986372,0.552554,0.776288};
   colourmap[20] = {0.986756,0.563819,0.765594};
   colourmap[21] = {0.987036,0.574940,0.754857};
   colourmap[22] = {0.987225,0.585954,0.744073};
   colourmap[23] = {0.987329,0.596830,0.733261};
   colourmap[24] = {0.987359,0.607585,0.722406};
   colourmap[25] = {0.987323,0.618208,0.711527};
   colourmap[26] = {0.987233,0.628715,0.700603};
   colourmap[27] = {0.987102,0.639094,0.689645};
   colourmap[28] = {0.986942,0.649357,0.678645};
   colourmap[29] = {0.986769,0.659488,0.667616};
   colourmap[30] = {0.986595,0.669503,0.656545};
   colourmap[31] = {0.986436,0.679384,0.645426};
   colourmap[32] = {0.986305,0.689155,0.634284};
   colourmap[33] = {0.986211,0.698787,0.623088};
   colourmap[34] = {0.986165,0.708312,0.611844};
   colourmap[35] = {0.986172,0.717716,0.600553};
   colourmap[36] = {0.986236,0.727001,0.589203};
   colourmap[37] = {0.986356,0.736180,0.577794};
   colourmap[38] = {0.986531,0.745264,0.566320};
   colourmap[39] = {0.986753,0.754257,0.554758};
   colourmap[40] = {0.987015,0.763154,0.543105};
   colourmap[41] = {0.987306,0.771982,0.531342};
   colourmap[42] = {0.987613,0.780738,0.519460};
   colourmap[43] = {0.987915,0.789439,0.507436};
   colourmap[44] = {0.988199,0.798073,0.495260};
   colourmap[45] = {0.988446,0.806659,0.482912};
   colourmap[46] = {0.988635,0.815182,0.470388};
   colourmap[47] = {0.988746,0.823650,0.457646};
   colourmap[48] = {0.988752,0.832060,0.444700};
   colourmap[49] = {0.988628,0.840393,0.431506};
   colourmap[50] = {0.988343,0.848634,0.418045};
   colourmap[51] = {0.987865,0.856767,0.404327};
   colourmap[52] = {0.987155,0.864764,0.390333};
   colourmap[53] = {0.986163,0.872591,0.376069};
   colourmap[54] = {0.984851,0.880197,0.361516};
   colourmap[55] = {0.983160,0.887543,0.346700};
   colourmap[56] = {0.981048,0.894578,0.331631};
   colourmap[57] = {0.978455,0.901226,0.316372};
   colourmap[58] = {0.975339,0.907443,0.300918};
   colourmap[59] = {0.971656,0.913155,0.285373};
   colourmap[60] = {0.967365,0.918306,0.269826};
   colourmap[61] = {0.962437,0.922845,0.254315};
   colourmap[62] = {0.956863,0.926719,0.238952};
   colourmap[63] = {0.950639,0.929902,0.223876};
   colourmap[64] = {0.943771,0.932362,0.209203};
   colourmap[65] = {0.936285,0.934099,0.195058};
   colourmap[66] = {0.928203,0.935119,0.181532};
   colourmap[67] = {0.919585,0.935443,0.168836};
   colourmap[68] = {0.910473,0.935102,0.156970};
   colourmap[69] = {0.900913,0.934145,0.146103};
   colourmap[70] = {0.890971,0.932627,0.136225};
   colourmap[71] = {0.880697,0.930612,0.127473};
   colourmap[72] = {0.870144,0.928148,0.119769};
   colourmap[73] = {0.859367,0.925314,0.113150};
   colourmap[74] = {0.848400,0.922162,0.107415};
   colourmap[75] = {0.837288,0.918740,0.102569};
   colourmap[76] = {0.826051,0.915112,0.098502};
   colourmap[77] = {0.814723,0.911316,0.095091};
   colourmap[78] = {0.803324,0.907388,0.092162};
   colourmap[79] = {0.791866,0.903355,0.089659};
   colourmap[80] = {0.780359,0.899246,0.087518};
   colourmap[81] = {0.768811,0.895079,0.085510};
   colourmap[82] = {0.757233,0.890864,0.083837};
   colourmap[83] = {0.745619,0.886624,0.082243};
   colourmap[84] = {0.733976,0.882352,0.080673};
   colourmap[85] = {0.722295,0.878067,0.079194};
   colourmap[86] = {0.710601,0.873765,0.077792};
   colourmap[87] = {0.698863,0.869457,0.076415};
   colourmap[88] = {0.687099,0.865139,0.075063};
   colourmap[89] = {0.675294,0.860818,0.073757};
   colourmap[90] = {0.663455,0.856488,0.072319};
   colourmap[91] = {0.651565,0.852155,0.071005};
   colourmap[92] = {0.639628,0.847804,0.069678};
   colourmap[93] = {0.627641,0.843452,0.068313};
   colourmap[94] = {0.615598,0.839095,0.066946};
   colourmap[95] = {0.603505,0.834729,0.065602};
   colourmap[96] = {0.591341,0.830355,0.064284};
   colourmap[97] = {0.579114,0.825970,0.063016};
   colourmap[98] = {0.566812,0.821583,0.061599};
   colourmap[99] = {0.554426,0.817184,0.060374};
   colourmap[100] = {0.541951,0.812777,0.059088};
   colourmap[101] = {0.529396,0.808369,0.057695};
   colourmap[102] = {0.516724,0.803948,0.056522};
   colourmap[103] = {0.503949,0.799524,0.055189};
   colourmap[104] = {0.491057,0.795083,0.053903};
   colourmap[105] = {0.478012,0.790640,0.052644};
   colourmap[106] = {0.464840,0.786193,0.051424};
   colourmap[107] = {0.451514,0.781731,0.050257};
   colourmap[108] = {0.438032,0.777261,0.049220};
   colourmap[109] = {0.424375,0.772770,0.048120};
   colourmap[110] = {0.410564,0.768269,0.047322};
   colourmap[111] = {0.396557,0.763754,0.046740};
   colourmap[112] = {0.382394,0.759216,0.046427};
   colourmap[113] = {0.368078,0.754661,0.046596};
   colourmap[114] = {0.353614,0.750060,0.047299};
   colourmap[115] = {0.339062,0.745430,0.048740};
   colourmap[116] = {0.324426,0.740751,0.050897};
   colourmap[117] = {0.309837,0.736016,0.054069};
   colourmap[118] = {0.295322,0.731217,0.058336};
   colourmap[119] = {0.281049,0.726343,0.063783};
   colourmap[120] = {0.267165,0.721371,0.070322};
   colourmap[121] = {0.253870,0.716308,0.077992};
   colourmap[122] = {0.241338,0.711121,0.086870};
   colourmap[123] = {0.229810,0.705805,0.096608};
   colourmap[124] = {0.219614,0.700353,0.107413};
   colourmap[125] = {0.210918,0.694751,0.118989};
   colourmap[126] = {0.203999,0.689009,0.131291};
   colourmap[127] = {0.198943,0.683096,0.144224};
   colourmap[128] = {0.195928,0.677034,0.157682};
   colourmap[129] = {0.194864,0.670809,0.171614};
   colourmap[130] = {0.195597,0.664426,0.185943};
   colourmap[131] = {0.197953,0.657910,0.200497};
   colourmap[132] = {0.201632,0.651248,0.215333};
   colourmap[133] = {0.206391,0.644457,0.230304};
   colourmap[134] = {0.211828,0.637562,0.245381};
   colourmap[135] = {0.217714,0.630567,0.260517};
   colourmap[136] = {0.223814,0.623480,0.275653};
   colourmap[137] = {0.229923,0.616316,0.290713};
   colourmap[138] = {0.235925,0.609097,0.305738};
   colourmap[139] = {0.241618,0.601816,0.320644};
   colourmap[140] = {0.246935,0.594504,0.335432};
   colourmap[141] = {0.251836,0.587158,0.350083};
   colourmap[142] = {0.256216,0.579800,0.364592};
   colourmap[143] = {0.260111,0.572415,0.378970};
   colourmap[144] = {0.263461,0.565010,0.393203};
   colourmap[145] = {0.266243,0.557617,0.407304};
   colourmap[146] = {0.268465,0.550209,0.421268};
   colourmap[147] = {0.270134,0.542804,0.435128};
   colourmap[148] = {0.271218,0.535387,0.448859};
   colourmap[149] = {0.271734,0.527982,0.462471};
   colourmap[150] = {0.271701,0.520573,0.476000};
   colourmap[151] = {0.271122,0.513171,0.489416};
   colourmap[152] = {0.269969,0.505756,0.502754};
   colourmap[153] = {0.268234,0.498366,0.516003};
   colourmap[154] = {0.265977,0.490956,0.529190};
   colourmap[155] = {0.263162,0.483520,0.542278};
   colourmap[156] = {0.259815,0.476095,0.555286};
   colourmap[157] = {0.255941,0.468642,0.568216};
   colourmap[158] = {0.251617,0.461165,0.581055};
   colourmap[159] = {0.246802,0.453652,0.593797};
   colourmap[160] = {0.241605,0.446094,0.606425};
   colourmap[161] = {0.236047,0.438486,0.618943};
   colourmap[162] = {0.230168,0.430809,0.631314};
   colourmap[163] = {0.224104,0.423032,0.643534};
   colourmap[164] = {0.217927,0.415169,0.655603};
   colourmap[165] = {0.211702,0.407184,0.667494};
   colourmap[166] = {0.205504,0.399058,0.679192};
   colourmap[167] = {0.199449,0.390795,0.690701};
   colourmap[168] = {0.193673,0.382343,0.702012};
   colourmap[169] = {0.188179,0.373723,0.713124};
   colourmap[170] = {0.182998,0.364881,0.724042};
   colourmap[171] = {0.178289,0.355845,0.734770};
   colourmap[172] = {0.173922,0.346566,0.745302};
   colourmap[173] = {0.169985,0.337070,0.755664};
   colourmap[174] = {0.166387,0.327320,0.765861};
   colourmap[175] = {0.163120,0.317346,0.775902};
   colourmap[176] = {0.160049,0.307119,0.785807};
   colourmap[177] = {0.157238,0.296666,0.795575};
   colourmap[178] = {0.154570,0.285952,0.805220};
   colourmap[179] = {0.152024,0.275051,0.814744};
   colourmap[180] = {0.149664,0.263945,0.824140};
   colourmap[181] = {0.147444,0.252636,0.833404};
   colourmap[182] = {0.145543,0.241182,0.842515};
   colourmap[183] = {0.144017,0.229604,0.851446};
   colourmap[184] = {0.143121,0.217998,0.860186};
   colourmap[185] = {0.143053,0.206386,0.868690};
   colourmap[186] = {0.144045,0.194808,0.876912};
   colourmap[187] = {0.146305,0.183361,0.884837};
   colourmap[188] = {0.150065,0.172187,0.892412};
   colourmap[189] = {0.155365,0.161400,0.899600};
   colourmap[190] = {0.162302,0.151032,0.906369};
   colourmap[191] = {0.170751,0.141359,0.912696};
   colourmap[192] = {0.180618,0.132437,0.918560};
   colourmap[193] = {0.191735,0.124462,0.923961};
   colourmap[194] = {0.203887,0.117652,0.928885};
   colourmap[195] = {0.216808,0.112135,0.933354};
   colourmap[196] = {0.230285,0.107942,0.937390};
   colourmap[197] = {0.244173,0.105248,0.941008};
   colourmap[198] = {0.258283,0.104034,0.944254};
   colourmap[199] = {0.272456,0.104173,0.947160};
   colourmap[200] = {0.286585,0.105597,0.949767};
   colourmap[201] = {0.300591,0.108073,0.952118};
   colourmap[202] = {0.314397,0.111532,0.954258};
   colourmap[203] = {0.327961,0.115651,0.956218};
   colourmap[204] = {0.341257,0.120310,0.958033};
   colourmap[205] = {0.354259,0.125435,0.959731};
   colourmap[206] = {0.366979,0.130836,0.961339};
   colourmap[207] = {0.379418,0.136365,0.962880};
   colourmap[208] = {0.391595,0.142096,0.964362};
   colourmap[209] = {0.403518,0.147862,0.965807};
   colourmap[210] = {0.415202,0.153628,0.967221};
   colourmap[211] = {0.426686,0.159414,0.968604};
   colourmap[212] = {0.437993,0.165152,0.969962};
   colourmap[213] = {0.449134,0.170885,0.971293};
   colourmap[214] = {0.460135,0.176498,0.972606};
   colourmap[215] = {0.471061,0.182036,0.973886};
   colourmap[216] = {0.481915,0.187511,0.975137};
   colourmap[217] = {0.492732,0.192804,0.976352};
   colourmap[218] = {0.503569,0.197982,0.977523};
   colourmap[219] = {0.514456,0.203022,0.978633};
   colourmap[220] = {0.525406,0.207899,0.979689};
   colourmap[221] = {0.536474,0.212587,0.980667};
   colourmap[222] = {0.547669,0.217088,0.981567};
   colourmap[223] = {0.559028,0.221368,0.982381};
   colourmap[224] = {0.570568,0.225422,0.983096};
   colourmap[225] = {0.582302,0.229249,0.983716};
   colourmap[226] = {0.594235,0.232843,0.984237};
   colourmap[227] = {0.606358,0.236227,0.984653};
   colourmap[228] = {0.618678,0.239354,0.984965};
   colourmap[229] = {0.631160,0.242265,0.985178};
   colourmap[230] = {0.643789,0.244961,0.985296};
   colourmap[231] = {0.656562,0.247499,0.985325};
   colourmap[232] = {0.669430,0.249819,0.985273};
   colourmap[233] = {0.682370,0.252030,0.985147};
   colourmap[234] = {0.695371,0.254083,0.984956};
   colourmap[235] = {0.708386,0.256007,0.984709};
   colourmap[236] = {0.721386,0.257863,0.984407};
   colourmap[237] = {0.734366,0.259680,0.984054};
   colourmap[238] = {0.747284,0.261432,0.983653};
   colourmap[239] = {0.760121,0.263191,0.983207};
   colourmap[240] = {0.772859,0.264941,0.982712};
   colourmap[241] = {0.785479,0.266761,0.982159};
   colourmap[242] = {0.797943,0.268669,0.981533};
   colourmap[243] = {0.810234,0.270687,0.980822};
   colourmap[244] = {0.822330,0.272847,0.980006};
   colourmap[245] = {0.834182,0.275246,0.979052};
   colourmap[246] = {0.845756,0.277908,0.977932};
   colourmap[247] = {0.857023,0.280889,0.976608};
   colourmap[248] = {0.867940,0.284304,0.975035};
   colourmap[249] = {0.878448,0.288192,0.973174};
   colourmap[250] = {0.888519,0.292629,0.970965};
   colourmap[251] = {0.898101,0.297682,0.968371};
   colourmap[252] = {0.907167,0.303395,0.965334};
   colourmap[253] = {0.915672,0.309845,0.961820};
   colourmap[254] = {0.923611,0.316999,0.957797};
   colourmap[255] = {0.930949,0.324889,0.953227};

   return;

}


} // end of namespace vdc
