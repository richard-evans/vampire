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
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

// program header
#include "vdc.hpp"

// openmp header
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

namespace vdc{

//------------------------------------------------------------------------------
// Function to output cells.inc file compatible with povray
//------------------------------------------------------------------------------
void output_cells_inc_file(unsigned int spin_file_id){

   // Open Povray Include File
	std::stringstream incpov_file_sstr;
	incpov_file_sstr << "cells-";
	incpov_file_sstr << std::setfill('0') << std::setw(8) << spin_file_id;
	incpov_file_sstr << ".inc";
	std::string incpov_file = incpov_file_sstr.str();

   // output informative message to user
   if(vdc::verbose) std::cout << "   Writing povray file " << incpov_file << "..." << std::flush;

   // open incfile
   std::ofstream incfile;
   incfile.open(incpov_file.c_str());

   const unsigned int tmid = vdc::materials.size();

   //---------------------------------------------------------------------------
   // parallelise stream formatting for better performance
   // step 1: parallel formatted output to stringstream in memory
   // step 2: binary write of formatted text to output file (awesomely fast!)
   //---------------------------------------------------------------------------
   #pragma omp parallel
   {

      std::stringstream otext;

      const double ux = vdc::cell_size[0]*0.5;
      const double uy = vdc::cell_size[1]*0.5;
      const double uz = vdc::cell_size[2]*0.5;

      // write to output text stream in parallel
      #pragma omp for
      for( unsigned int cell = 0; cell < total_cells; cell++){

         // get magnetization for colour contrast
         double mx = vdc::cell_magnetization[cell][tmid][0];
         double my = vdc::cell_magnetization[cell][tmid][1];
         double mz = vdc::cell_magnetization[cell][tmid][2];

         // temporary thread private variables defining spin colours
         double red=0.0, green=0.0, blue=1.0;

         // calculate rgb components based on spin orientation
         vdc::rgb(mx, my, mz, red, green, blue);

         double cx = vdc::cell_coords[3*cell + 0] - vdc::system_centre[0];
         double cy = vdc::cell_coords[3*cell + 1] - vdc::system_centre[1];
         double cz = vdc::cell_coords[3*cell + 2] - vdc::system_centre[2];

         // format text for povray file
         otext << "cell"<< "(" <<
                  cx << "," << cy << "," << cz << "," <<
                  cx-ux << "," << cy-uy << "," << cz-uz << "," <<
                  cx+ux << "," << cy+uy << "," << cz+uz << "," <<
                  mx << "," << my << "," << mz << "," <<
                  red << "," << green << "," << blue << ")\n";

      } // end of parallel for

      // force each thread to write to file in order
      #pragma omp critical
      incfile << otext.str();

   } // end of parallel region

   // flush data to include file and close
   incfile << std::flush;
   incfile.close();

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

//------------------------------------------------------------------------------
// Function to output spins.pov file compatible with povray
//------------------------------------------------------------------------------
void output_povray_cells_file(){

	std::ofstream pfile;
	pfile.open("cells.pov");

   // Calculate location of camera
   std::vector<double> dim = {vdc::system_size[0]+0.001, vdc::system_size[1]+0.001, vdc::system_size[2]+0.001};
   std::vector<double> vec(3);

   // direction camera looks (normalised)
   // technically only position if lookat is not (0,0,0)
	double size = sqrt(dim[0]*dim[0] + dim[1]*dim[1] + dim[2]*dim[2]);
   if (default_camera_pos){
      vec[0] = (1.0/dim[0]);
      vec[1] = (1.0/dim[1]);
      vec[2] = (1.0/dim[2]);
   }
   else { vec = vdc::camera_pos; }

   // normalise camera position vector
   double mag_vec = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
   vec[0]/=mag_vec;
   vec[1]/=mag_vec;
   vec[2]/=mag_vec;

   pfile << "//-------------------------------------------------------------------------" << std::endl;
   pfile << "// Povray file generated using vampire" << std::endl;
   pfile << "//-------------------------------------------------------------------------" << std::endl;
   pfile << "#version 3.5;"            << std::endl;
	pfile << "#include \"colors.inc\""  << std::endl;
	pfile << "#include \"metals.inc\""	<< std::endl;
	pfile << "#include \"screen.inc\""	<< std::endl;
   // look at position
	pfile << "#declare LX=" << vdc::camera_look_at[0]*dim[0]/2.0 << ";" << std::endl;
	pfile << "#declare LY=" << vdc::camera_look_at[1]*dim[1]/2.0 << ";" << std::endl;
	pfile << "#declare LZ=" << vdc::camera_look_at[2]*dim[2]/2.0 << ";" << std::endl;
   // camera position
	pfile << "#declare CX=" << size*vec[0]*6.0*vdc::camera_zoom << ";" << std::endl;
	pfile << "#declare CY=" << size*vec[1]*6.0*vdc::camera_zoom << ";" << std::endl;
	pfile << "#declare CZ=" << size*vec[2]*6.0*vdc::camera_zoom << ";" << std::endl;
	pfile << "#declare ref=0.05;" << std::endl;
	pfile << "global_settings { assumed_gamma 2.0 }" << std::endl;
   // background colour
	pfile << "background { color " << vdc::background_colour << " }" << std::endl;

	pfile << "Set_Camera(<CX,CY,CZ>, <LX,LY,LZ>, 15)" << std::endl;
	pfile << "Set_Camera_Aspect(4,3)"  << std::endl;
	pfile << "Set_Camera_Sky(<0,0,1>)" << std::endl;
	pfile << "light_source { <2*CX, 2*CY, 2*CZ> color White}" << std::endl;

   pfile << "#declare Initial_Frame = " << vdc::start_file_id << ";" << std::endl;
   pfile << "#declare Final_Frame = "   << vdc::final_file_id << ";" << std::endl;

   // sscale affects the spin arrow
   pfile << "#declare mscale" << "=" << 2.0 << ";" << std::endl;
   // rscale affects sphere(atom) size
   // cscale affects cube size
   pfile << "#declare spx = 0.1; // spacing between cells" << std::endl;
   pfile << "#declare spy = 0.1;" << std::endl;
   pfile << "#declare spz = 0.1;" << std::endl;
   pfile << "#declare cones"   << " = true;" << std::endl;
   pfile << "#declare arrows"  << " = false;" << std::endl;
   pfile << "#declare cubes"   << " = true;" << std::endl;
   pfile << "#declare mcolors" << " = true;" << std::endl;
   pfile << "#declare mcolor"  << " = pigment {color rgb < 0.1 0.1 0.1 >};" << std::endl;
   pfile << "#macro cell"<< "(cx,cy,cz,sx,sy,sz,ex,ey,ez,mx,my,mz,cr,cg,cb)" << std::endl;
   pfile << "union{"      << std::endl;
   pfile << "#if(cubes) box {<sx+spx*0.5,sy+spy*0.5,sz+spz*0.5>,<ex-spx*0.5,ey-spy*0.5,ez-spz*0.5>} #end" << std::endl;
   pfile << "#if(cones) cone {<cx+0.5*mx*mscale,cy+0.5*my*mscale,cz+0.5*mz*mscale>,0.0 <cx-0.5*mx*mscale,cy-0.5*my*mscale,cz-0.5*mz*mscale>,mscale*0.5} #end" << std::endl;
   pfile << "#if(arrows) cylinder {<cx+mx*0.5*mscale,cy+my*0.5*mscale,cz+mz*0.5*mscale>," <<
                                  "<cx-mx*0.5*mscale,cy-my*0.5*mscale,cz-mz*0.5*mscale>,mscale*0.12}";
   pfile << "            cone {<cx+mx*0.5*1.6*mscale, cy+my*0.5*1.6*mscale, cz+mz*0.5*1.6*mscale>,mscale*0.0" <<
                              "<cx+mx*0.5*mscale,     cy+my*0.5*mscale,     cz+mz*0.5*mscale    >,mscale*0.2} #end" << std::endl;
   pfile << "#if(mcolors) texture { pigment {color rgb <cr cg cb>}finish {reflection {ref} diffuse 1 ambient 0}}" << std::endl;
   pfile << "#else texture { mcolor finish {reflection {ref} diffuse 1 ambient 0}} #end" << std::endl;
   pfile << "}"    << std::endl;
   pfile << "#end" << std::endl;

   // frame specific povray output
	pfile << "#include concat(\"cells-\", str(frame_number, -8, 0) \".inc\")" << std::endl;

   // close output file
	pfile.close();

   //---------------------------------------------------------------------------
   // output povray ini file for rendering all files by default
   //---------------------------------------------------------------------------
   std::ofstream pifile;
	pifile.open("spins.ini");

   pifile << "Input_File_Name = \"spins.pov\"" << std::endl;
   pifile << "Width = 800"                     << std::endl;
   pifile << "Height = 600"                    << std::endl;
   pifile << "Antialias = On"                  << std::endl;
   pifile << "Antialias_Threshold = 0.3"       << std::endl;
   pifile << "Output_File_Type = N"            << std::endl;
   pifile << "Initial_Frame = " << vdc::start_file_id << std::endl;
   pifile << "Final_Frame = "   << vdc::final_file_id << std::endl;

   pifile.close();

   return;

}

}
