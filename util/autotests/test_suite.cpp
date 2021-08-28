#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <sys/param.h>
#include <unistd.h>

// global scope command variable for executing vampire
std::string cmd;

bool exchange_test(std::string dir, double result){

   // get root directory
   std::string path = std::filesystem::current_path();

   std::cout << "Testing exchange energy for " << dir << "    : " << std::flush;

   // change directory
   int test = chdir(dir.c_str());
   if(test !=0){
      std::cerr << "Error changing into directory " << dir << ". Returning as failed test." << std::endl;
      return false;
   }

   // run vampire
   system(cmd.c_str());

   std::string line;

   // open output file
   std::ifstream ifile;
   ifile.open("output");

   // read value after header
   for(int i=0; i<9; i++) getline(ifile, line);
   std::stringstream liness(line);
   double value = 0.0;
   liness >> value;

   // return to parent directory
   chdir(path.c_str());

   // now test value obtained from code
   const double ratio = value/result;

   if(ratio >0.99999 && ratio < 1.00001){
      std::cout << "OK" << std::endl;
      return true;
   }
   else{
      std::cout << "FAIL | expected: " << result << "\tobtained:  " << value << "\t" << line << std::endl;
      return false;
   }

   return false;

}

bool integrator_test(std::string dir, double rx, double ry, double rz){

   // get root directory
   std::string path = std::filesystem::current_path();

   std::cout << "Testing integrator for " << dir << "   \t: " << std::flush;

   // change directory
   int test = chdir(dir.c_str());
   if(test !=0){
      std::cerr << "Error changing into directory " << dir << ". Returning as failed test." << std::endl;
      return false;
   }

   // run vampire
   system(cmd.c_str());

   std::string line;

   // open output file
   std::ifstream ifile;
   ifile.open("output");

   // read value after header
   for(int i=0; i<982; i++) getline(ifile, line);
   std::stringstream liness(line);
   double v1 = 0.0;
   double vx = 0.0;
   double vy = 0.0;
   double vz = 0.0;
   double vm = 0.0;

   liness >> v1 >> vx >> vy >> vz >> vm;

   // return to parent directory
   chdir(path.c_str());

   // now test value obtained from code
   const double ratiox = vx/rx;
   const double ratioy = vy/ry;
   const double ratioz = vz/rz;

   if(ratiox >0.99999 && ratiox < 1.00001 && ratioy >0.99999 && ratioy < 1.00001 && ratioz >0.99999 && ratioz < 1.00001){
      std::cout << "OK" << std::endl;
      return true;
   }
   else{
      std::cout << "FAIL | expected: " << rx << "\t" << ry << "\t" << rz << "\t" << "\tobtained:  " << vx << "\t" << vy << "\t" << vz << "\t" << "\t" << line << std::endl;
      return false;
   }

   return false;

}

int main(){

   std::cout << "--------------------------------------------------" << std::endl;
   std::cout << "    Running system test suite for vampire code" << std::endl;
   std::cout << "--------------------------------------------------" << std::endl;

   // Get vampire path and command
   std::string path = std::filesystem::current_path();
   cmd = path+"/../vampire-serial 1>/dev/null";

   unsigned int fail = 0;

   // Exchange energy tests
   if( !exchange_test("crystals/SC" , -3.0e-17 ) ) fail += 1;
   if( !exchange_test("crystals/FCC", -2.4e-16 ) ) fail += 1;

   // Integrator tests
   if( !integrator_test("dynamics/Heun",-0.106813,-0.337996,0.935067) ) fail += 1;

   // Summary
   std::cout << "--------------------------------------------------" << std::endl;
   if(fail >0){
      std::cout << "Failed " << fail << " tests : OVERALL FAIL" << std::endl;
   }
   else{
      std::cout << "Failed " << fail << " tests : OVERALL PASS" <<	std::endl;
   }
   std::cout << "--------------------------------------------------" << std::endl;

   return 0;

}
