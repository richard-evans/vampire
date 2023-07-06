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
// Initialise Povray colourwheel
//------------------------------------------------------------------------------
void initialise_povray() {

   vdc::initialise_colourwheel();

}

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
      for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_atoms_list[i];

         // get magnetization for colour contrast
         double sx = spins[3*atom+0];
         double sy = spins[3*atom+1];
         double sz = spins[3*atom+2];

         // flip antiferromagnetic spins if required
         if (std::find(afm_materials.begin(), afm_materials.end(), vdc::type[atom]+1) != afm_materials.end() ){
            sx = -sx;
            sy = -sy;
            sz = -sz;
         }

         // temporary thread private variables defining spin colours
         double red=0.0, green=0.0, blue=1.0;

         // calculate rgb components based on spin orientation
         vdc::rgb(sx, sy, sz, red, green, blue);

         // format text for povray file
         otext << "spinm"<< type[atom]+1 << "(" <<
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
      for(size_t i=0; i < vdc::sliced_nm_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_nm_atoms_list[i];

         // format text for povray file
         otext << "spinm"<< nm_type[atom]+1 << "(" <<
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

   const double camr = size*6.0*vdc::camera_zoom;
   const double camp = 180.0*acos(vec[1]/sqrt(vec[0]*vec[0]+vec[1]*vec[1]))/M_PI;
   const double camt = 180.0*acos(vec[2]/sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]))/M_PI;

   pfile << "//-------------------------------------------------------------------------" << std::endl;
   pfile << "// Povray file generated using vampire" << std::endl;
   pfile << "//-------------------------------------------------------------------------" << std::endl;
   pfile << "#version 3.5;"            << std::endl;
	pfile << "#include \"colors.inc\""  << std::endl;
	pfile << "#include \"metals.inc\""	<< std::endl;
	pfile << "#include \"screen.inc\""	<< std::endl;
   // look at position
   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
   pfile << "// Camera parameters" << std::endl;
   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
	pfile << "#declare LX = " << vdc::camera_look_at[0]*dim[0]/2.0 << "; // camera looking position" << std::endl;
	pfile << "#declare LY = " << vdc::camera_look_at[1]*dim[1]/2.0 << ";" << std::endl;
	pfile << "#declare LZ = " << vdc::camera_look_at[2]*dim[2]/2.0 << ";" << std::endl;
   pfile << "// camera location" << std::endl;
   pfile << "#declare cam_theta  = " << camt << "; // angle from z in degrees" << std::endl;
   pfile << "#declare cam_phi    = " << camp << "; // angle from x in degrees" << std::endl;
   pfile << "#declare cam_radius = " << camr << "; // distance from origin" << std::endl;
   pfile << "#declare CX = cam_radius * cos(cam_phi*pi/180.0) * sin(cam_theta*pi/180.0);" << std::endl;
   pfile << "#declare CY = cam_radius * sin(cam_phi*pi/180.0) * sin(cam_theta*pi/180.0);" << std::endl;
   pfile << "#declare CZ = cam_radius * cos(cam_theta*pi/180.0);" << std::endl;
   // camera position
	//pfile << "#declare CX = " << size*vec[0]*6.0*vdc::camera_zoom << "; // camera location" << std::endl;
	//pfile << "#declare CY = " << size*vec[1]*6.0*vdc::camera_zoom << ";" << std::endl;
	//pfile << "#declare CZ = " << size*vec[2]*6.0*vdc::camera_zoom << ";" << std::endl;
   pfile << "Set_Camera(<CX,CY,CZ>, <LX,LY,LZ>, 15)" << std::endl;
	pfile << "Set_Camera_Aspect(4,3)"  << std::endl;
	pfile << "Set_Camera_Sky(<0,0,1>)" << std::endl;

   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
   pfile << "// Global constants and appearance" << std::endl;
   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
	pfile << "global_settings { assumed_gamma 2.0 }" << std::endl;
   // background colour
	pfile << "background { color " << vdc::background_colour << " } // background  colour" << std::endl;


	pfile << "light_source { <2*CX, 2*CY, 2*CZ> color White} // lights" << std::endl;
   pfile << "#declare Initial_Frame = " << vdc::start_file_id << ";" << std::endl;
   pfile << "#declare Final_Frame = "   << vdc::final_file_id << ";\n" << std::endl;

   pfile << "//------------------------------------------" << std::endl;
   pfile << "// Global settings for appearance" << std::endl;
   pfile << "//------------------------------------------\n" << std::endl;

   pfile << "#declare ref = 0.05; // reflection level of objects" << std::endl;
   pfile << "#declare dif = 1.0; // diffuse level of objects" << std::endl;
   pfile << "#declare amb = 0.1; // ambient level of objects\n" << std::endl;

   pfile << "// spin scale (default 2.0)" << std::endl;
   pfile << "#declare global_spin_scale = true;" << std::endl;
   pfile << "#declare sscale = 2.0;\n" << std::endl;

   pfile << "// radius scale (default 1.2)" << std::endl;
   pfile << "#declare global_radius_scale = true;" << std::endl;
   pfile << "#declare rscale = 1.2;\n" << std::endl;

   pfile << "// cube scale (default 1.2)" << std::endl;
   pfile << "#declare global_cube_scale = true;" << std::endl;
   pfile << "#declare cscale = 3.54;\n" << std::endl;

   pfile << "// global objects" << std::endl;
   pfile << "#declare global_cones   = false;" << std::endl;
   pfile << "#declare global_arrows  = true;" << std::endl;
   pfile << "#declare global_spheres = true;" << std::endl;
   pfile << "#declare global_cubes   = false;\n" << std::endl;

   pfile << "// non-magnetic element colours" << std::endl;
   pfile << "#declare global_non_magnet_colour = true;" << std::endl;
   pfile << "#declare nmcolour = pigment {color rgb < 0.2 0.2 0.2 >};\n" << std::endl;

   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
   pfile << "// Colour macro" << std::endl;
   pfile << "// The default is vdc generated and takes rgb values, but can be overridden to" << std::endl;
   pfile << "// apply non-linear colour scales, colour maps etc." << std::endl;
   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
   pfile << "#macro spinrgb(sx, sy, sz, cr, cg, cb)" << std::endl;
   pfile << "   pigment {color rgb <cr cg cb>}" << std::endl;
   pfile << "#end\n" << std::endl;

   pfile << "// alternative colour schemes" << std::endl;
   pfile << "//#include \"util/povray_colours/jet.inc\"" << std::endl;
   pfile << "//#include \"util/povray_colours/purple_white.inc\"" << std::endl;
   pfile << "//#include \"util/povray_colours/blue_gold.inc\"" << std::endl;
   pfile << "//#include \"util/povray_colours/color_wheel.inc\"\n" << std::endl;

   //---------------------------------------------------------------------------
   // Determine non-magnetic materials looping over all non-magnetic atoms
   //---------------------------------------------------------------------------
   std::vector<bool> is_nm_mat(vdc::materials.size(),false);
   for(size_t i=0; i < vdc::sliced_nm_atoms_list.size(); i++){

      // get atom ID
      unsigned int atom = vdc::sliced_nm_atoms_list[i];

      const int mat = vdc::nm_type[atom];
      is_nm_mat[mat] = true;

   }

   // check there are shape sizes defined for all materials
   while (vdc::atom_sizes.size()  < vdc::materials.size()){ vdc::atom_sizes.push_back(1.2);  }
   while (vdc::arrow_sizes.size() < vdc::materials.size()){ vdc::arrow_sizes.push_back(2.0); }

   //for( auto i : remove_materials){
   for(int i=0; i< remove_materials.size(); i++){
      std::cout << i << "\t" << remove_materials[i] << std::endl;
      //std::cout << std::find(remove_materials.begin(), remove_materials.end(), ) << std::endl;
   }

   //-----------------------------------------------------
   // optionally output sticks macro if requested in vdc
   //-----------------------------------------------------
   /*if(vdc::povsticks){ // version with truncated cylinders
      pfile << "\n" << std::endl;
      pfile << "//----------------------------------------------------------------" << std::endl;
      pfile << "// Sticks macro" << std::endl;
      pfile << "//----------------------------------------------------------------" << std::endl;
      pfile << "#macro stick(sx,sy,sz,ex,ey,ez,r1,r2)" << std::endl;
      pfile << "   #declare dx = ex-sx;" << std::endl;
      pfile << "   #declare dy = ey-sy;" << std::endl;
      pfile << "   #declare dz = ez-sz;" << std::endl;
      pfile << "   #declare r = sqrt(dx*dx + dy*dy + dz*dz);" << std::endl;
      pfile << "   #declare xh = 0.5*dx/r;" << std::endl;
      pfile << "   #declare yh = 0.5*dy/r;" << std::endl;
      pfile << "   #declare zh = 0.5*dz/r;" << std::endl;
      pfile << "   difference{" << std::endl;
      pfile << "      cylinder {" << std::endl;
      pfile << "         <sx+xh*r1,sy+yh*r1,sz+zh*r1>," << std::endl;
      pfile << "         <ex-xh*r2,ey-yh*r2,ez-zh*r2>, 0.2" << std::endl;
      pfile << "      }" << std::endl;
      pfile << "      sphere{<sx,sy,sz>, r1*0.5}" << std::endl;
      pfile << "      sphere{<ex,ey,ez>, r2*0.5}" << std::endl;
      pfile << "      texture { pigment {color rgb < 0.1 0.1 0.1 >} finish {reflection ref diffuse dif ambient amb } }" << std::endl;
      pfile << "   }" << std::endl;
      pfile << "#end\n" << std::endl;
   }*/

   if(vdc::povsticks){
      pfile << "\n" << std::endl;
      pfile << "//----------------------------------------------------------------" << std::endl;
      pfile << "// Sticks macro" << std::endl;
      pfile << "//----------------------------------------------------------------" << std::endl;
      pfile << "#macro stick(sx,sy,sz,ex,ey,ez,r1,r2)" << std::endl;
      pfile << "   difference{" << std::endl;
      pfile << "      cylinder {" << std::endl;
      pfile << "         <sx,sy,sz>," << std::endl;
      pfile << "         <ex,ey,ez>, 0.2" << std::endl;
      pfile << "      }" << std::endl;
      pfile << "      sphere{<sx,sy,sz>, r1*0.5}" << std::endl;
      pfile << "      sphere{<ex,ey,ez>, r2*0.5}" << std::endl;
      pfile << "      texture { pigment {color rgb < 0.1 0.1 0.1 >} finish {reflection ref diffuse dif ambient amb } }" << std::endl;
      pfile << "   }" << std::endl;
      pfile << "#end\n" << std::endl;
   }

   // Output material specific macros
	for(unsigned int indx=0; indx < vdc::materials.size(); indx++){

      const int imat = indx+1; // correct for standard material numbering as in vampire material file

      // check for removed material
      if (std::find(remove_materials.begin(), remove_materials.end(), imat) == remove_materials.end() ){

         pfile << "//-----------------------------------------------------------------------------------" << std::endl;
         pfile << "// Material " << vdc::materials[indx].id+1 << "\tName: " << vdc::materials[indx].name << "\tElement: " << vdc::materials[indx].element << std::endl;
         pfile << "//-----------------------------------------------------------------------------------" << std::endl;

         // if material is magnetic
         if(is_nm_mat[indx] == false){


            pfile << "#if(global_spin_scale) #declare sscale"  << imat << " = sscale;" << std::endl;
            pfile << "#else #declare sscale" << imat << " = " << vdc::arrow_sizes[indx] << "; #end" << std::endl;
            pfile << "#if(global_radius_scale) #declare rscale"  << imat << " = rscale;" << std::endl;
            pfile << "#else #declare rscale" << imat << " = " << vdc::atom_sizes[indx]  << "; #end" << std::endl;
            pfile << "#if(global_cube_scale) #declare cscale"  << imat << " = cscale;" << std::endl;
            pfile << "#else #declare cscale" << imat << " = 3.54; #end" << std::endl;

            pfile << "#if(global_cones)   #declare cones"   << imat << " = global_cones;   #else #declare cones"   << imat << " = false; #end" << std::endl;
            pfile << "#if(global_cubes)   #declare arrows"  << imat << " = global_arrows;  #else #declare arrows"  << imat << " = true;  #end" << std::endl;
            pfile << "#if(global_spheres) #declare spheres" << imat << " = global_spheres; #else #declare spheres" << imat << " = true;  #end" << std::endl;
            pfile << "#if(global_cubes)   #declare cubes"   << imat << " = global_cubes;   #else #declare cubes"   << imat << " = false; #end" << std::endl;
            pfile << "#declare spincolors" << imat << " = true; // enable colours defined in vdc" << std::endl;
            pfile << "#declare spincolor" << imat << "  = pigment {color rgb < 0.1 0.1 0.1 >};" << std::endl;
            pfile << "//-------------------------------------" << std::endl;
            pfile << "#macro spinm"<< imat << "(cx,cy,cz,sx,sy,sz, cr,cg,cb)" << std::endl;
            pfile << "union{"      << std::endl;
            pfile << "   #if(spheres" << imat << ") sphere {<cx,cy,cz>,0.5*rscale"<< imat << "} #end" << std::endl;
            pfile << "   #if(cubes"   << imat << ") box {<cx-cscale"<< imat << "*0.5,cy-cscale" << imat << "*0.5,cz-cscale"<< imat << "*0.5>,<cx+cscale"<< imat << "*0.5,cy+cscale" << imat << "*0.5,cz+cscale"<< imat << "*0.5>} #end" << std::endl;
            pfile << "   #if(cones"   << imat << ") cone {<cx+0.5*sx*sscale" << imat << ",cy+0.5*sy*sscale"<< imat << ",cz+0.5*sz*sscale"<< imat << ">,0.0 <cx-0.5*sx*sscale"<< imat << ",cy-0.5*sy*sscale"<< imat << ",cz-0.5*sz*sscale"<< imat << ">,sscale" << imat << "*0.5} #end" << std::endl;
            pfile << "   #if(arrows"  << imat << ")" << std::endl;
            pfile << "      cylinder {<cx+sx*0.5*sscale"<< imat <<",    cy+sy*0.5*sscale"<< imat <<",    cz+sz*0.5*sscale"<< imat << ">" << std::endl;
            pfile << "                <cx-sx*0.5*sscale"<< imat <<",    cy-sy*0.5*sscale"<< imat <<",    cz-sz*0.5*sscale"<< imat <<">,sscale"<< imat <<"*0.12}" << std::endl;
            pfile << "      cone     {<cx+sx*0.5*1.6*sscale"<< imat <<",cy+sy*0.5*1.6*sscale"<< imat <<",cz+sz*0.5*1.6*sscale"<< imat <<">,sscale"<< imat <<"*0.0" << std::endl;
            pfile << "                <cx+sx*0.5*sscale"<< imat <<",    cy+sy*0.5*sscale"<< imat <<",    cz+sz*0.5*sscale"<< imat <<"    >,sscale"<< imat <<"*0.2}" << std::endl;
            pfile << "   #end" << std::endl;
            pfile << "   #if(spincolors"<< imat << ") texture { spinrgb(sx,sy,sz,cr,cg,cb) finish {reflection ref diffuse dif ambient amb } }" << std::endl;
            pfile << "   #else texture { spincolor"<< imat << " finish {reflection ref diffuse dif ambient amb }} #end" << std::endl;
            pfile << "}"    << std::endl;
            pfile << "#end\n" << std::endl;
         }
         else{
            pfile << "#if(global_spheres) #declare spheres" << imat << " = global_spheres; #else #declare spheres" << imat << " = true;  #end" << std::endl;
            pfile << "#if(global_cubes)   #declare cubes"   << imat << " = global_cubes;   #else #declare cubes"   << imat << " = false; #end" << std::endl;
            pfile << "#if(global_radius_scale) #declare rscale"  << imat << " = rscale;" << std::endl;
            pfile << "#else #declare rscale" << imat << " = " << vdc::atom_sizes[indx]  << "; #end" << std::endl;
            pfile << "#if(global_cube_scale) #declare cscale"  << imat << " = cscale;" << std::endl;
            pfile << "#else #declare cscale" << imat << " = 3.54; #end" << std::endl;
            pfile << "#declare spincolors"<< imat << " = false;" << std::endl;
            pfile << "#declare spincolor" << imat << " = pigment {color rgb < 0.2 0.2 0.2 >};" << std::endl;
            pfile << "//-------------------------------------" << std::endl;
            pfile << "#macro spinm"  << imat << "(cx,cy,cz,sx,sy,sz,cr,cg,cb)" << std::endl;
            pfile << "union{" << std::endl;
            pfile << "   #if(spheres"   << imat << ") sphere {<cx,cy,cz>,0.5*rscale"<< imat << "} #end" << std::endl;
            pfile << "   #if(cubes"     << imat << ") box {<cx-cscale"<< imat << "*0.5,cy-cscale" << imat << "*0.5,cz-cscale"<< imat << "*0.5>," << std::endl;
            pfile << "                        <cx+cscale"<< imat << "*0.5,cy+cscale" << imat << "*0.5,cz+cscale"<< imat << "*0.5>} #end" << std::endl;
            pfile << "   #if(spincolors"<< imat << ") texture { pigment {color rgb <cr cg cb>}finish {reflection {ref} diffuse 1 ambient 0}}" << std::endl;
            pfile << "   #else texture { spincolor"<< imat << " finish {reflection {ref} diffuse 1 ambient 0}} #end" << std::endl;
            pfile << "}"    << std::endl;
            pfile << "#end" << std::endl;
         }
      } // end of check if material is to be removed
	}

   // frame specific povray output
   pfile << "//----------------------------------------------------------------" << std::endl;
   pfile << "// Include spin and sticks data" << std::endl;
   pfile << "//----------------------------------------------------------------" << std::endl;
   pfile << "#include concat(\"spins-\", str(frame_number, -8, 0) \".inc\")" << std::endl;

   // optionally include sticks
   if(vdc::povsticks) pfile << "#include \"sticks.inc\"" << std::endl;

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
