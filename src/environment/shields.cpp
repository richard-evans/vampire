//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// Vampire headers
#include "environment.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "create.hpp"
#include "vio.hpp"
#include "errors.hpp"
#include <algorithm>
// C++ headers
#include <math.h>

namespace env = environment::internal;

// typesafe sign function
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}


namespace environment{

namespace internal{

bool  in_x(double x, double z){

   const double xr = x - 500.0; // shift to zero
   const double zr = z - 300.0; // shift to base of side shield

   const double zmin = 200.0 * exp(-((fabs(xr)-190)*0.01));
   //std::cout << x << '\t' << z << '\t'  << xr << '\t' << zr << '\t' << zmin << std::endl;
   if(zr > zmin && z < 500.0) return true;
   return false;

}

//------------------------------------------------------------------------------
// Function to calculate basic shield geometry for reader
//------------------------------------------------------------------------------
int in_shield(double x, double y, double z,int shield){

//
//
//   else if (expoential_shields){
//
//    // height of inner sensor region
//    const double stack_height = 200; // Angstroms
//
//    const double xr = x;
//    const double yr = y;
//    const double zr = z; // reduced height
//
// //   Bottom shield
  //  if(z <= 302 && z >= -2)    return 1;
  // //
  // //if (zr < 310  && zr > 290) return true;
  // //
  // // // Top shield (++z) 31-51 nm
  // if(z >= 518 && z <= 722.9) return 2;
  // //
  // //  //Top shield (+z) 52-72 nm
//  std::cout << env::shield_shape[shield] << '\t' << x << '\t' << y << '\t' << z << '\t' << env::shield_min_x[shield] << '\t' << env::shield_max_x[shield] << '\t' << env::shield_min_y[shield] << '\t' << env::shield_max_y[shield] << '\t' << env::shield_min_z[shield] << '\t' << env::shield_max_z[shield] <<std::endl;
  if(env::shield_shape[shield] == "cube" && x >= env::shield_min_x[shield] && x <= env::shield_max_x[shield] &&
      y >= env::shield_min_y[shield] && y <= env::shield_max_y[shield] &&
      z >= env::shield_min_z[shield] && z <= env::shield_max_z[shield]) {
        return 1;
      }
  else if (env::shield_shape[shield] == "exponential"){

  //  std::cout << "EXPONENTIAL CELLS!!!" << std::endl;
    double xmin = env::shield_min_x[shield];
    double xmax = env::shield_max_x[shield];
    double zmin = env::shield_min_z[shield];
    double zmax = env::shield_max_z[shield];

    //double f = exp((x-xmax)*0.01);
    //double f2 = exp((-x+(xmax-xmin)-xmax)*0.01);
    //double fmin = exp((xmin-xmax)*0.01);
    //double fmax = exp((xmax-xmax)*0.01);
    double g;// = zmin+(zmax-zmin)*(f-fmin)/(fmax-fmin);
    double g2;// = zmin+(zmax-zmin)*(f2-fmin)/(fmax-fmin);
   if (x !=xmax){
      g2 = -10*(zmax-zmin)/(x-xmax) + zmin;
    }
    else g2 = 100000;

    if (x !=xmin){
    g =  10*(zmax-zmin)/(x-xmin) + zmin;
    }
    else g = 1000000;


   // std::cout <<  g2 <<"\t" <<   z << "\t" << zmax << '\t' <<  zmin << "\t" << xmax << "\t" << x << std::endl;

     if(g < z && env::pos_or_neg[shield] == "pos" && x >= env::shield_min_x[shield] && x <= env::shield_max_x[shield] &&
         y >= env::shield_min_y[shield] && y <= env::shield_max_y[shield] &&
         z >= env::shield_min_z[shield] && z <= env::shield_max_z[shield]){//}&& z <= zmax && x <= xmax && x >= xmin && z >= zmin && env::pos_or_neg[shield] == "pos") {
    //  std::cout << "neg" <<std::endl;
      //     std::cout << x << '\t'  << z << '\t' <<g2 << '\t' << g << '\t' <<  std::endl;
      return 1;
     }
     else if( g2 < z && env::pos_or_neg[shield] == "neg" && x >= env::shield_min_x[shield] && x <= env::shield_max_x[shield] &&
         y >= env::shield_min_y[shield] && y <= env::shield_max_y[shield] &&
         z >= env::shield_min_z[shield] && z <= env::shield_max_z[shield]){ //z <= zmax && x <= xmax && x >= xmin && z >= zmin) {
   //    std::cout << "neg" <<std::endl;
  //    std::cout << x << '\t'  << z << '\t' <<g2 << '\t' << g << '\t' <<  std::endl;

      // std::cout << zmin << '\t' << zmax << "\t" << xmin << '\t' << xmax << '\t' << g2 << '\t' << g << '\t' << x << '\t'  << z << '\t' << std::endl;
       return 1;
     }
     else{
    //          std::cout << "not in" <<std::endl;
              return 0;
     }


  }
  // // //  // side shields
  // if(in_x(x, z)) return 4;


  //if (sqrt(x*x + y*y + z*z)  < 90) return 1;
  //else
else   return 0;

}

int bias_shields(){

  double shield_Ms = 1;
  double x_size = dim[1];
  double y_size = 1000000;
  double z_size = dim[2];

  double x_pos = x_size/2.0;
  double y_pos;
  double z_pos = dim[2]/2.0;

  double y_pos_1 = -y_size/2.0;
  double y_pos_2 =  y_size/2.0 +dim[0];


   double prefactor = shield_Ms;///(4.0*M_PI);
  //save this new m as the initial value, so it can be saved and used in the final equation.
    for (int cell = 0; cell < num_cells; cell ++){

    std::vector <double > B(3,0.0);
    bias_field_x[cell] = 0;
    bias_field_y[cell] = 0;
    bias_field_z[cell] = 0;

     //cell position in Angstrom
     double x_cell = cell_coords_array_x[cell];
     double y_cell = cell_coords_array_y[cell];
     double z_cell = cell_coords_array_z[cell];

     const double xb = x_size * 0.5;
     const double yb = y_size * 0.5;
     const double zb = z_size * 0.5;

     for (int shield = 0; shield < 2; shield++){

       if (shield == 0) y_pos = y_pos_1;
       if (shield == 1) y_pos = y_pos_2;
       //calculates the vector in A from the cell to the shields
       double x = sqrt((x_cell - x_pos)*(x_cell - x_pos));
       double y = sqrt((y_cell - y_pos)*(y_cell - y_pos));
       double z = sqrt((z_cell - z_pos)*(z_cell - z_pos));

       double Bx = 0.0;
       double By = 0.0;
       double Bz = 0.0;

       for(int k=1; k<3; k++){

           // predefine power as fixed for loop iteration
           const double m1k = pow(-1,k);

           for(int l=1; l<3; l++){

              // predefine power as fixed for loop iteration
              const double m1l = pow(-1,l);

              for(int m=1; m<3; m++){

                 const double m1m = pow(-1,m);
                 const double m1klm = pow(-1,k+l+m);

                 const double xp = x + xb*m1k;
                 const double yp = y + yb*m1l;
                 const double zp = z + zb*m1m;

                 const double xabs = fabs(xp);
                 const double yabs = fabs(yp);

                 double r = sqrt(xp*xp + yp*yp + zp*zp);

                 Bx = Bx + m1klm * log(zp + r);
                 By = By - m1klm * sign(yp) * sign(xp) * atan(xabs * zp / (yabs * r));
                 Bz = Bz + m1klm * log(xp + r);


              }
           }
       }
       bias_field_x[cell] = bias_field_x[cell] + By*prefactor;
       bias_field_y[cell] = bias_field_y[cell] + Bx*prefactor;
       bias_field_z[cell] = bias_field_z[cell] + Bz*prefactor;

     }
   // / std::cout << bias_field_x[cell] << '\t' << bias_field_y[cell] << '\t' << bias_field_z[cell] << std::endl;

  }
//std::cin.get();





  return 0;

}

int read_in_shield_info(){

  std::ifstream ifile;
  ifile.open("shield_geom");

  //int shield_number;
  string shield_type;

  if (ifile.good()){
    std::cout << "Creating shields" << std::endl;
  }

  //-------------------------------------------------------
  // Material 0
  //-------------------------------------------------------


  // Print informative message to zlog file
  zlog << zTs() << "Creating shields " << std::endl;

  int line_counter=0;
  // Loop over all lines and pass keyword to matching function
  while (! ifile.eof() ){
      line_counter++;
      // read in whole line
      std::string line;
      getline(ifile,line);

      // save a copy of the line before stripping characters in case of error
      std::string original_line = line;

      // Clear whitespace, quotes and tabs
      line.erase(remove(line.begin(), line.end(), '\t'), line.end());
      line.erase(remove(line.begin(), line.end(), ' '), line.end());
      line.erase(remove(line.begin(), line.end(), '\"'), line.end());

      // remove carriage returns for dos formatted files
                  line.erase(remove(line.begin(), line.end(), '\r'), line.end());

      // strip key,word,unit,value
      std::string key="";
      std::string word="";
      std::string value="";
      std::string unit="";
      std::string index="";
      int super_index=1; // Inital values *as would be read from input file*
      int sub_index=1;

      // get size of string
      int linelength = line.length();
      int last=0;

      // set character triggers
      const char* colon=":";	// Word identifier
      const char* eq="=";		// Value identifier
      const char* exc="!";		// Unit identifier
      const char* hash="#";	// Comment identifier
      const char* si="[";		// Index identifier
      const char* ei="]";		// End index identifier

      // Determine key and super index by looping over characters in line
      for(int i=0;i<linelength;i++){
          char c=line.at(i);
          last=i;



          // if character is not ":" or "=" or "!" or "#" interpret as key
          if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash) && (c != *si) && (c != *ei)){
              key.push_back(c);
          }
          // Check for number of materials statement
          else if(c == *eq){
              // break to read in value
              break;
          }
          // Check for super_index
          else if(c ==*si){
              const int old=last;
              // Get super index
              for(int j=old+1;j<linelength;j++){
                  c=line.at(j);
                  if(c != *ei){
                      index.push_back(c);
                  }
                  else{
                      break;
                  }
                  last=j;
              }

              // check for valid index
              super_index = atoi(index.c_str());
              if((super_index>=1) && (super_index<env::num_shields+1)){
                  break;
              }
              else{
                  std::cerr << "Invalid index number " << index << " on line " << line_counter << " in material input file" << std::endl;
                  std::cerr << "Causes could be invalid character or outside of range, ie less than 1 or greater than max_materials=" << mp::max_materials << ", exiting" << std::endl;
                  err::vexit();
              }

          }
          // For anything else
          else break;
      }
      const int end_key=last;

      //
      //err::vexit();
      // Determine the rest
      for(int i=end_key;i<linelength;i++){


          char c=line.at(i);
          //std::cout << c << std::endl;
           // colon found - interpret as word
          if(c== *colon){
              for(int j=i+1;j<linelength;j++){

                  // if character is not special add to value
                  char c=line.at(j);
                  //std::cout << c << std::endl;
                  if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
                      // check for sub-index
                      if(c == *si){
                          index="";
                          while(line.at(j+1) != *ei){
                              j++;
                              index.push_back(line.at(j));
                          }
                          sub_index=atoi(index.c_str());
                          // Check for valid index
                          if((sub_index<1) || (sub_index>=mp::max_materials)){
                              std::cerr << "Invalid sub-index number " << index << " on line " << line_counter << " in material input file" << std::endl;
                              std::cerr << "Causes could be invalid character or outside of range, ie less than 1 or greater than max_materials=" << mp::max_materials << ", exiting" << std::endl;
                              err::vexit();
                          }
                          // end of word
                          break;
                      }
                    else  word.push_back(c);

                  }
                  // if character is special then go back to main loop
                  else{

                      i=j-1;
                    //  std::cout << i << '\t' << line.at(i) <<std::endl;
                      break;
                  }
              }

          }
          // equals found - interpret as value
          else if(c== *eq){

              for(int j=i+1;j<linelength;j++){
                  // if character is not special add to value
                  char c=line.at(j);
                  if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
                      value.push_back(c);

                  }
                  //if character is special then go back to main loop
                  else{
                      i=j-1;
                      break;
                  }
              }
          }
          // exclaimation mark found - interpret as unit
          else if(c== *exc){
              for(int j=i+1;j<linelength;j++){
                  // if character is not special add to value
                  char c=line.at(j);
                  if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
                      unit.push_back(c);
                  }
                  // if character is special then go back to main loop
                  else{
                      i=j-1;
                      break;
                  }
              }
          }
          // hash found - interpret as comment
          else if(c== *hash){
              break;
          }


      }

      std::string test="shield-type";
      if(word==test){
        if (value == "cube"){
          env::shield_shape[super_index-1] = value;
        }
        else if (value == "exponential"){
          env::shield_shape[super_index-1] = value;
        }
      }
      test="which-side-would-you-like-exp";
      if(word==test){
        env::pos_or_neg[super_index-1] = value;

      }


      test="minimum-x";
      if(word==test){
        double g=atof(value.c_str());
        vin::check_for_valid_value(g, word, 1, "environment", unit, "length", -1e10, 1e10,"shield_geom","-100 - 100 cms");
        env::shield_min_x[super_index-1] = g;
      //  std::cout << word << '\t' << env::shield_min_x[super_index-1] << '\t' << super_index - 1 <<std::endl;

      }

      test="minimum-y";
      if(word==test){
        double g=atof(value.c_str());
        vin::check_for_valid_value(g, word, 1, "environment", unit, "length", -1e10, 1e10,"shield_geom","-100 - 100 cms");
        env::shield_min_y[super_index-1] = g;

      }
      test="minimum-z";
      if(word==test){
        double g=atof(value.c_str());
        vin::check_for_valid_value(g, word, 1, "environment", unit, "length", -1e10, 1e10,"shield_geom","-100 - 100 cms");
        env::shield_min_z[super_index-1] = g;

      }

      test="maximum-x";
      if(word==test){
        double g=atof(value.c_str());
         vin::check_for_valid_value(g, word, 1, "environment", unit, "length", -1e10, 1e10,"shield_geom","-100 - 100 cms");
        env::shield_max_x[super_index-1] = g;

      }

      test="maximum-y";
      if(word==test){
        double g=atof(value.c_str());
        vin::check_for_valid_value(g, word, 1, "environment", unit, "length", -1e10, 1e10,"shield_geom","-100 - 100 cms");
        env::shield_max_y[super_index-1] = g;

      }
      test="maximum-z";
      if(word==test){
        double g=atof(value.c_str());
        vin::check_for_valid_value(g, word, 1, "environment", unit, "length", -1e10, 1e10,"shield_geom","-100 - 100 cms");
        env::shield_max_z[super_index-1] = g;

      }

      test="maximum-cell-size";
      if(word==test){
        double g=atof(value.c_str());
        vin::check_for_valid_value(g, word, 1, "environment", unit, "length", 20, 2000,"shield_geom","2 nm to 100 nm");
        env::shield_max_cell_size[super_index-1] = g;

      }

      test="minimum-cell-size";
      if(word==test){
        double g=atof(value.c_str());
         vin::check_for_valid_value(g, word, 1, "environment", unit, "length", 10, 10,"shield_geom","1 nm to 10 nm");
        env::shield_min_cell_size[super_index-1] = g;

      }

      test="Ms";
      if(word==test){
        double g=atof(value.c_str());
        vin::check_for_valid_positive_value(g, word, 1, "environment", unit, "magnetisation", 1e-35,1e-15 ,"shield_geom","1e-32 - 1e-18");
        env::shield_ms[super_index-1] = g;

      }

      test="Ku";
      if(word==test){
        double g=atof(value.c_str());
        vin::check_for_valid_value(g, word, 1, "environment", unit, "anisotropy", 0,1e-15 ,"shield_geom","0 - 1e-18");
        env::shield_ku[super_index-1] = g;
      }

      test="Tc";
      if(word==test){
        double g=atof(value.c_str());
        vin::check_for_valid_positive_value(g, word, 1, "environment", unit, "none", 10, 30000 ,"shield_geom","0 -30000");
        env::shield_Tc[super_index-1] = g;
      }

      test="A";
      if(word==test){
        double g=atof(value.c_str());
        vin::check_for_valid_value(g, word, 1, "environment", unit, "exchange", 0, 1e19 ,"shield_geom","0 - 1e-19");
        env::shield_A[super_index-1][sub_index-1] = g;
        //std::cout << super_index << '\t' << sub_index << '\t' << g << std::endl;
      }
      test="applied-field-strength";
      if(word==test){
         double K=atof(value.c_str());
         env::H_strength[super_index] = K;
      }
      test="applied-field-unit-vector";
      if(word==test){
         // temporary storage container
         std::vector<double> u(3);

         // read values from string
         u=vin::doubles_from_string(value);
         // Copy sanitised unit vector to material
         env::shield_Hext_x[super_index] =u.at(0);
         env::shield_Hext_y[super_index] =u.at(1);
         env::shield_Hext_z[super_index] =u.at(2);
      }
      test="initial-spin-direction";
      if(word==test){
      //    // first test for random spins
      //    test="random";
      //    if(value==test){
      //       env::random_spins[super_index]=true;
      //    }
      //    else{
            // temporary storage container
            std::vector<double> u(3);
            // read values from string
            u=vin::doubles_from_string(value);
            std::cout << super_index << '\t' << u.at(0) <<std::endl;
          // vin::check_for_valid_value(u.at(0), word, line, prefix, unit, "none", 0,1,"input","0- 1");
          // vin::check_for_valid_value(u.at(1), word, line, prefix, unit, "none", 0,1,"input","0- 1");
          // vin::check_for_valid_value(u.at(2), word, line, prefix, unit, "none", 0,1,"input","0- 1");

            // Copy sanitised unit vector to material
            env::initial_spin_x[super_index-1]=u.at(0);
            env::initial_spin_y[super_index-1]=u.at(1);
            env::initial_spin_z[super_index-1]=u.at(2);

            // ensure random spins is unset
            env::random_spins[super_index-1]=false;
         }
      //    // return
      //    return true;
      // }



    //  string empty="";
    //  std::cout << key << std::endl;
    //  if(key!=empty){
          //std::cout << key << "[" << super_index << "]:" << word << "[" << sub_index << "]=" << value << " !" << unit << std::endl;
          //std::cout << "\t" << "key:  " << key << std::endl;
          //std::cout << "\t" << "word: " << word << std::endl;
          //std::cout << "\t" << "value:" << value << std::endl;
          //std::cout << "\t" << "unit: " << unit << std::endl;
    //  int matchcheck = vin::match_material(word, value, unit, line_counter, super_index-1, sub_index-1, original_line, ifile);
    //      if(matchcheck==EXIT_FAILURE){
    //          err::vexit();
    //       }
     }

     return 0;

  }



    //std::string line; // declare a string to hold line of text
    // while(getline(ifile,line) ){
    //   std::stringstream line_stream(line);
    //
    //   line_stream >> shield_number >> shield_type;
    //   if (shield_type == "cube"){
    //
    //   }
    //
    // }


  }
}
