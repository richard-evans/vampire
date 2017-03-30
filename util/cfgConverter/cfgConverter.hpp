

#include <string>
#include <vector>

struct material_t{

   double mx,my,mz,mu_s,magm;

};


class cells
{
public:
   cells(int numx, int numy, int numz, double maxx, double maxy, double maxz);
   ~cells();

   void addspin(double cx, double cy, double cz, double sx, double sy, double sz);

   std::vector <double> outcx();
   std::vector <double> outcy();
   std::vector <double> outcz();
   std::vector <double> outsx();
   std::vector <double> outsy();
   std::vector <double> outsz();

private:
   void binning(double cx, double cy, double cz, double sx, double sy, double sz);
   std::vector <int> cellxEdge;
   std::vector <int> cellyEdge;
   std::vector <int> cellzEdge;

   std::vector < std::vector < std::vector < double > > > cellsx;
   std::vector < std::vector < std::vector < double > > > cellsy;
   std::vector < std::vector < std::vector < double > > > cellsz;

   std::vector < std::vector < std::vector < int > > > bincount;
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
