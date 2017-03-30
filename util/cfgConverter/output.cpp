/// Output to povray

#include "cfgConverter.hpp"

void outpov()
{

   //Ouput povray file header
   double size, mag_vec;
   double vec[3];

   size = sqrt(dim[0] * dim[0] + dim[1] * dim[1] + dim[2] * dim[2]);
   vec[0] = (1.0 / dim[0]);
   vec[1] = (1.0 / dim[1]);
   vec[2] = (1.0 / dim[2]);
   mag_vec = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
   vec[0] /= mag_vec;
   vec[1] /= mag_vec;
   vec[2] /= mag_vec;

   pfile << "#include \"colors.inc\"" << std::endl;
   pfile << "#include \"metals.inc\"" << std::endl;
   pfile << "#include \"screen.inc\"" << std::endl;
   pfile << "#declare LX=" << dim[0] * 0.5 << ";" << std::endl;
   pfile << "#declare LY=" << dim[1] * 0.5 << ";" << std::endl;
   pfile << "#declare LZ=" << dim[2] * 0.5 << ";" << std::endl;
   pfile << "#declare CX=" << size * vec[0] * 6.0 << ";" << std::endl;
   pfile << "#declare CY=" << size * vec[1] * 6.0 << ";" << std::endl;
   pfile << "#declare CZ=" << size * vec[2] * 6.0 << ";" << std::endl;
   pfile << "#declare ref=0.05;" << std::endl;
   pfile << "global_settings { assumed_gamma 2.0 }" << std::endl;
   pfile << "background { color Gray30 }" << std::endl;

   pfile << "Set_Camera(<CX,CY,CZ>, <LX,LY,LZ>, 15)" << std::endl;
   pfile << "Set_Camera_Aspect(4,3)" << std::endl;
   pfile << "Set_Camera_Sky(<0,0,1>)" << std::endl;
   pfile << "light_source { <2*CX, 2*CY, 2*CZ> color White}" << std::endl;

   for (int imat = 0; imat < n_mat; imat++)
   {
      pfile << "#declare sscale" << imat << "=2.0;" << std::endl;
      pfile << "#declare rscale" << imat << "=1.2;" << std::endl;
      pfile << "#declare cscale" << imat << "=3.54;" << std::endl;
      pfile << "#declare cones" << imat << "=0;" << std::endl;
      pfile << "#declare arrows" << imat << "=1;" << std::endl;
      pfile << "#declare spheres" << imat << "=1;" << std::endl;
      pfile << "#declare cubes" << imat << "=0;" << std::endl;
      pfile << "#declare spincolors" << imat << "=1;" << std::endl;
      pfile << "#declare spincolor" << imat << "=pigment {color rgb < 0.1 0.1 0.1 >};" << std::endl;
      pfile << "#macro spinm" << imat << "(cx,cy,cz,sx,sy,sz, cr,cg,cb)" << std::endl;
      pfile << "union{" << std::endl;
      pfile << "#if(spheres" << imat << ") sphere {<cx,cy,cz>,0.5*rscale" << imat << "} #end" << std::endl;
      pfile << "#if(cubes" << imat << ") box {<cx-cscale" << imat << "*0.5,cy-cscale" << imat << "*0.5,cz-cscale" << imat << "*0.5>,<cx+cscale" << imat << "*0.5,cy+cscale" << imat << "*0.5,cz+cscale" << imat << "*0.5>} #end" << std::endl;
      pfile << "#if(cones" << imat << ") cone {<cx+0.5*sx*sscale0,cy+0.5*sy*sscale" << imat << ",cz+0.5*sz*sscale" << imat << ">,0.0 <cx-0.5*sx*sscale" << imat << ",cy-0.5*sy*sscale" << imat << ",cz-0.5*sz*sscale" << imat << ">,sscale0*0.5} #end" << std::endl;
      pfile << "#if(arrows" << imat << ") cylinder {<cx+sx*0.5*sscale" << imat << ",cy+sy*0.5*sscale" << imat << ",cz+sz*0.5*sscale" << imat << ">,<cx-sx*0.5*sscale" << imat << ",cy-sy*0.5*sscale" << imat << ",cz-sz*0.5*sscale" << imat << ">,sscale" << imat << "*0.12}";
      pfile << "cone {<cx+sx*0.5*1.6*sscale" << imat << ",cy+sy*0.5*1.6*sscale" << imat << ",cz+sz*0.5*1.6*sscale" << imat << ">,sscale" << imat << "*0.0 <cx+sx*0.5*sscale" << imat << ",cy+sy*0.5*sscale" << imat << ",cz+sz*0.5*sscale" << imat << ">,sscale" << imat << "*0.2} #end" << std::endl;
      pfile << "#if(spincolors" << imat << ") texture { pigment {color rgb <cr cg cb>}finish {reflection {ref} diffuse 1 ambient 0}}" << std::endl;
      pfile << "#else texture { spincolor" << imat << " finish {reflection {ref} diffuse 1 ambient 0}} #end" << std::endl;
      pfile << "}" << std::endl;
      pfile << "#end" << std::endl;
   }
   pfile << "#include \"" << incpov_file_sstr.str() << "\"" << std::endl;

   pfile.close();

}

void outputinc(){

   std::ofstream incpfile;
   incpfile.open(incpov_file.c_str());

   double sx, sy, sz, red, green, blue, ireal;
   unsigned int si = 0;

   // Read in spin coordinates and output to povray file
   for (int i = 0; i < n_local_atoms; i++)
   {
      spinfile >> sx >> sy >> sz;
      rgb(sz, red, green, blue);
      incpfile << "spinm" << mat[si] << "(" << cx[si] << "," << cy[si] << "," << cz[si] << ","
               << sx << "," << sy << "," << sz << "," << red << "," << green << "," << blue << ")" << std::endl;
      si++;
   }

   // close master file
   spinfile.close();

   // now read subsidiary files
   for (int file = 0; file < n_files; file++)
   {
      std::ifstream infile;
      infile.open(filenames[file].c_str());

      // read number of atoms in this file
      getline(infile, dummy);
      n_local_atoms = atoi(dummy.c_str());
      for (int i = 0; i < n_local_atoms; i++)
      {
         infile >> sx >> sy >> sz;
         rgb(sz, red, green, blue);
         incpfile << "spinm" << mat[si] << "(" << cx[si] << "," << cy[si] << "," << cz[si] << ","
                  << sx << "," << sy << "," << sz << "," << red << "," << green << "," << blue << ")" << std::endl;
         si++;
      }
      // close subsidiary file
      infile.close();
   }
   	// close povray inc file
	incpfile.close();
}