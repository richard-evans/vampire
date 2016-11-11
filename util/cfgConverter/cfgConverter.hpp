

#include <string>
#include <vector>

struct material_t{

   double mx,my,mz,mu_s,magm;

};


class data
{
public:
   data(int argc, char* argv[]);
   ~data();
private:
   std::vector <int> mat;
   std::vector <int> cat; 
   std::vector <double> cx;
   std::vector <double> cy;
   std::vector <double> cz;
   std::vector <std::string> type;
};

class cfg
{
public:
   cfg(int argc, char* argv[]);
   ~cfg();
   void read();
   void output();

private:
   void outputinc();
   void outputpv();
   void gencells();
   data atoms;
   data cells;

};
