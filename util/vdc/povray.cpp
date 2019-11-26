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

// forward function declarations

//------------------------------------------------------------------------------
// Function to output spins.inc file compatible with povray
//------------------------------------------------------------------------------
void output_inc_file(unsigned int spin_file_id){

   // Open Povray Include File
	std::stringstream incpov_file_sstr;
	incpov_file_sstr << "spins-";
	incpov_file_sstr << std::setfill('0') << std::setw(8) << spin_file_id;
	incpov_file_sstr << ".inc";
	std::string incpov_file = incpov_file_sstr.str();

   // output informative message to user
   if(vdc::verbose) std::cout << "   Writing povray file " << incpov_file << "..." << std::flush;

   // temporary variables defining spin colours
   double red=0.0, green=0.0, blue=1.0;
   
   // open incfile
   std::ofstream incfile;
   incfile.open(incpov_file.c_str());

   //---------------------------------------------------------------------------
   // parallelise stream formatting for better performance
   // step 1: parallel formatted output to stringstream in memory
   // step 2: binary write of formatted text to output file (awesomely fast!)
   //---------------------------------------------------------------------------
   #pragma omp parallel
   {

      std::stringstream otext;

      // write to output text stream in parallel
      #pragma omp for
      for( auto &atom : vdc::sliced_atoms_list ){

         // flip antiferromagnetic spins if required
         if (std::find(afm_materials.begin(), afm_materials.end(), vdc::type[atom]+1) != afm_materials.end() ){
            spins[3*atom+0] = -spins[3*atom+0];
            spins[3*atom+1] = -spins[3*atom+1];
            spins[3*atom+2] = -spins[3*atom+2];
         }

         // get magnetization for colour contrast
         const double sx = spins[3*atom+0];
         const double sy = spins[3*atom+1];
         const double sz = spins[3*atom+2];

         // calculate rgb components based on spin orientation
         vdc::rgb(sx, sy, sz, red, green, blue);

         // format text for povray file
         otext << "spinm"<< type[atom] << "(" <<
                  coordinates[3*atom+0]-vdc::system_centre[0] << "," << coordinates[3*atom+1]-vdc::system_centre[1] << "," << coordinates[3*atom+2]-vdc::system_centre[2] << "," <<
                  spins[3*atom+0] << "," << spins[3*atom+1] << "," << spins[3*atom+2] << "," <<
                  red << "," << green << "," << blue << ")\n";

      } // end of parallel for

      // force each thread to write to file in order
      #pragma omp critical
      incfile << otext.str();

   } // end of parallel region

   //---------------------------------------------------------------------------
   // write non-magnetic atoms to inc file
   //---------------------------------------------------------------------------
   // parallelise stream formatting for better performance
   //---------------------------------------------------------------------------
   #pragma omp parallel
   {

      std::stringstream otext;

      // write to output text stream in parallel
      #pragma omp for
      for( auto &atom : vdc::sliced_nm_atoms_list ){

         // format text for povray file
         otext << "spinm"<< nm_type[atom] << "(" <<
                  nm_coordinates[3*atom+0]-vdc::system_centre[0] << "," <<
                  nm_coordinates[3*atom+1]-vdc::system_centre[1] << "," <<
                  nm_coordinates[3*atom+2]-vdc::system_centre[2] << "," <<
                  0.0 << "," << 0.0 << "," << 0.0 << "," << // no spin
                  0.3 << "," << 0.3 << "," << 0.3 << ")\n"; // grey colour by default

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
void output_povray_file(){

	std::ofstream pfile;
	pfile.open("spins.pov");

   // Calculate location of camera
   double dim[3] = {vdc::system_size[0]+0.001, vdc::system_size[1]+0.001, vdc::system_size[2]+0.001};
   double vec[3];

	double size = sqrt(dim[0]*dim[0] + dim[1]*dim[1] + dim[2]*dim[2]);
	vec[0] = (1.0/dim[0]);
	vec[1] = (1.0/dim[1]);
	vec[2] = (1.0/dim[2]);
	double mag_vec = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	vec[0]/=mag_vec;
	vec[1]/=mag_vec;
	vec[2]/=mag_vec;

   pfile << "//-------------------------------------------------------------------------" << std::endl;
   pfile << "// Povray file generated using vampire" << std::endl;
   pfile << "//-------------------------------------------------------------------------" << std::endl;
   pfile << "#version 3.5;" << std::endl;
	pfile << "#include \"colors.inc\"" << std::endl;
	pfile << "#include \"metals.inc\""	<< std::endl;
	pfile << "#include \"screen.inc\""	<< std::endl;
	pfile << "#declare LX=0.0;" << std::endl;
	pfile << "#declare LY=0.0;" << std::endl;
	pfile << "#declare LZ=0.0;" << std::endl;
	pfile << "#declare CX=" << size*vec[0]*6.0 << ";" << std::endl;
	pfile << "#declare CY=" << size*vec[1]*6.0 << ";" << std::endl;
	pfile << "#declare CZ=" << size*vec[2]*6.0 << ";" << std::endl;
	pfile << "#declare ref=0.05;" << std::endl;
	pfile << "global_settings { assumed_gamma 2.0 }" << std::endl;
	pfile << "background { color Gray30 }" << std::endl;

	pfile << "Set_Camera(<CX,CY,CZ>, <LX,LY,LZ>, 15)" << std::endl;
	pfile << "Set_Camera_Aspect(4,3)" << std::endl;
	pfile << "Set_Camera_Sky(<0,0,1>)" << std::endl;
	pfile << "light_source { <2*CX, 2*CY, 2*CZ> color White}" << std::endl;

   pfile << "#declare Initial_Frame = " << vdc::start_file_id << ";" << std::endl;
   pfile << "#declare Final_Frame = " << vdc::final_file_id << ";" << std::endl;

   // Determine non-magnetic materials looping over all non-magnetic atoms
   std::vector<bool> is_nm_mat(vdc::materials.size(),false);
   for( auto &atom : vdc::sliced_nm_atoms_list ){
      const int mat = vdc::nm_type[atom];
      is_nm_mat[mat] = true;
   }

   // Output material specific macros
	for(unsigned int imat=0; imat < vdc::materials.size(); imat++){
      if (std::find(remove_materials.begin(), remove_materials.end(), imat+1) == remove_materials.end() ){
         if(is_nm_mat[imat] == false){
      		pfile << "#declare sscale"<< imat << "=2.0;" << std::endl;
      		pfile << "#declare rscale"<< imat << "=1.2;" << std::endl;
      		pfile << "#declare cscale"<< imat << "=3.54;" << std::endl;
      		pfile << "#declare cones"<< imat << "=0;" << std::endl;
      		pfile << "#declare arrows"<< imat << "=1;" << std::endl;
      		pfile << "#declare spheres"<< imat << "=1;" << std::endl;
      		pfile << "#declare cubes" << imat << "=0;" << std::endl;
      		pfile << "#declare spincolors"<< imat << "=1;" << std::endl;
      		pfile << "#declare spincolor"<< imat << "=pigment {color rgb < 0.1 0.1 0.1 >};" << std::endl;
      		pfile << "#macro spinm"<< imat << "(cx,cy,cz,sx,sy,sz, cr,cg,cb)" << std::endl;
      		pfile << "union{" << std::endl;
      		pfile << "#if(spheres" << imat << ") sphere {<cx,cy,cz>,0.5*rscale"<< imat << "} #end" << std::endl;
      		pfile << "#if(cubes" << imat << ") box {<cx-cscale"<< imat << "*0.5,cy-cscale" << imat << "*0.5,cz-cscale"<< imat << "*0.5>,<cx+cscale"<< imat << "*0.5,cy+cscale" << imat << "*0.5,cz+cscale"<< imat << "*0.5>} #end" << std::endl;
      		pfile << "#if(cones"<< imat << ") cone {<cx+0.5*sx*sscale" << imat << ",cy+0.5*sy*sscale"<< imat << ",cz+0.5*sz*sscale"<< imat << ">,0.0 <cx-0.5*sx*sscale"<< imat << ",cy-0.5*sy*sscale"<< imat << ",cz-0.5*sz*sscale"<< imat << ">,sscale" << imat << "*0.5} #end" << std::endl;
      		pfile << "#if(arrows" << imat << ") cylinder {<cx+sx*0.5*sscale"<< imat <<",cy+sy*0.5*sscale"<< imat <<",cz+sz*0.5*sscale"<< imat <<
      					">,<cx-sx*0.5*sscale"<< imat <<",cy-sy*0.5*sscale"<< imat <<",cz-sz*0.5*sscale"<< imat <<">,sscale"<< imat <<"*0.12}";
      		pfile << "cone {<cx+sx*0.5*1.6*sscale"<< imat <<",cy+sy*0.5*1.6*sscale"<< imat <<",cz+sz*0.5*1.6*sscale"<< imat <<">,sscale"<< imat <<"*0.0 <cx+sx*0.5*sscale"<< imat <<
      					",cy+sy*0.5*sscale"<< imat <<",cz+sz*0.5*sscale"<< imat <<">,sscale"<< imat <<"*0.2} #end" << std::endl;
      		pfile << "#if(spincolors"<< imat << ") texture { pigment {color rgb <cr cg cb>}finish {reflection {ref} diffuse 1 ambient 0}}" << std::endl;
      		pfile << "#else texture { spincolor"<< imat << " finish {reflection {ref} diffuse 1 ambient 0}} #end" << std::endl;
      		pfile << "}" << std::endl;
      		pfile << "#end" << std::endl;
         }
         else{
            pfile << "#declare rscale"<< imat << "=1.2;" << std::endl;
            pfile << "#declare cscale"<< imat << "=0.1;" << std::endl;
            pfile << "#declare spheres"<< imat << "=1;" << std::endl;
            pfile << "#declare cubes" << imat << "=1;" << std::endl;
            pfile << "#declare spincolors"<< imat << "=1;" << std::endl;
            pfile << "#declare spincolor"<< imat << "=pigment {color rgb < 0.1 0.1 0.1 >};" << std::endl;
            pfile << "#macro spinm"<< imat << "(cx,cy,cz,sx,sy,sz,cr,cg,cb)" << std::endl;
            pfile << "union{" << std::endl;
            pfile << "#if(spheres" << imat << ") sphere {<cx,cy,cz>,0.5*rscale"<< imat << "} #end" << std::endl;
            pfile << "#if(cubes" << imat << ") box {<cx-cscale"<< imat << "*0.5,cy-cscale" << imat << "*0.5,cz-cscale"<< imat << "*0.5>,<cx+cscale"<< imat << "*0.5,cy+cscale" << imat << "*0.5,cz+cscale"<< imat << "*0.5>} #end" << std::endl;
            pfile << "#if(spincolors"<< imat << ") texture { pigment {color rgb <cr cg cb>}finish {reflection {ref} diffuse 1 ambient 0}}" << std::endl;
            pfile << "#else texture { spincolor"<< imat << " finish {reflection {ref} diffuse 1 ambient 0}} #end" << std::endl;
            pfile << "}" << std::endl;
            pfile << "#end" << std::endl;
         }
      }
	}
   // frame specific povray output
	pfile << "#include concat(\"spins-\", str(frame_number, -8, 0) \".inc\")" << std::endl;

   // close output file
	pfile.close();

   //---------------------------------------------------------------------------
   // output povray ini file for rendering all files by default
   //---------------------------------------------------------------------------
   std::ofstream pifile;
	pifile.open("spins.ini");

   pifile << "Input_File_Name = \"spins.pov\"" << std::endl;
   pifile << "Width = 800" << std::endl;
   pifile << "Height = 600" << std::endl;
   pifile << "Antialias = On" << std::endl;
   pifile << "Antialias_Threshold = 0.3" << std::endl;
   pifile << "Output_File_Type = N" << std::endl;
   pifile << "Initial_Frame = " << vdc::start_file_id << std::endl;
   pifile << "Final_Frame = " << vdc::final_file_id << std::endl;

   pifile.close();

   return;

}

}
