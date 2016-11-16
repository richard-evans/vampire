/// Output to povray

#include "cfgConverter.hpp"

//--------------------------Public---------------------------------

cells::cells(int numx, int numy, int numz, double maxx, double maxy, double maxz){

   bincount.resize(numx);
   cellsx.resize(numx);
   cellsy.resize(numx);
   cellsz.resize(numx);
   for (int i = 0; i < numx; ++i) {
      bincount[i].resize(numy);
      cellsx[i].resize(numy);
      cellsy[i].resize(numy);
      cellsz[i].resize(numy);

      for (int j = 0; j < numy; ++j){
         bincount[i][j].resize(numz);
         cellsx[i][j].resize(numz);
         cellsy[i][j].resize(numz);
         cellsz[i][j].resize(numz);
      }
   }

   cellxEdge.resize(numx);
   for (int i = 0; i < numx; i++)
      cellxEdge[i] = (maxx/numx) * i;

   cellyEdge.resize(numy);
   for (int i = 0; i < numy; i++)
      cellyEdge[i] = (maxy/numy) * i;

   cellzEdge.resize(numz);
   for (int i = 0; i < numz; i++)
      cellzEdge[i] = (maxz/numz) * i;


}

void cells::addspin(double cx, double cy, double cz, double sx, double sy, double sz){
   binning(cx, cy, cz, sx, sy, sz);
}


std::vector <double> cells::outcx(){
   std::vector <double> Output;

   for (int i = 0; i < cellsx.size(); i++){
      for (int j = 0; j < cellsx[0].size(); j++){
         for (int k = 0; k < cellsx[0][0].size(); k++){

            Output.push_back((cellxEdge[i]+cellxEdge[i+1])/2);

         }
      }
   }

   return Output;
}
std::vector <double> cells::outcy(){
   std::vector <double> Output;

   for (int i = 0; i < cellsy.size(); i++){
      for (int j = 0; j < cellsy[0].size(); j++){
         for (int k = 0; k < cellsy[0][0].size(); k++){

            Output.push_back((cellyEdge[j]+cellyEdge[j+1])/2);

         }
      }
   }

   return Output;
}
std::vector <double> cells::outcz(){
   std::vector <double> Output;

   for (int i = 0; i < cellsz.size(); i++){
      for (int j = 0; j < cellsz[0].size(); j++){
         for (int k = 0; k < cellsz[0][0].size(); k++){

            Output.push_back((cellzEdge[k]+cellzEdge[k+1])/2);

         }
      }
   }

   return Output;
}
std::vector <double> cells::outsx(){
   std::vector <double> Output;

   for (int i = 0; i < cellsx.size(); i++){
      for (int j = 0; j < cellsx[0].size(); j++){
         for (int k = 0; k < cellsx[0][0].size(); k++){

            Output.push_back(cellsx[i][j][k]);

         }
      }
   }

   return Output;
}
std::vector <double> cells::outsy(){
   std::vector <double> Output;

   for (int i = 0; i < cellsy.size(); i++){
      for (int j = 0; j < cellsy[0].size(); j++){
         for (int k = 0; k < cellsy[0][0].size(); k++){

            Output.push_back(cellsy[i][j][k]);
            
         }
      }
   }

   return Output;
}
std::vector <double> cells::outsz(){
   std::vector <double> Output;

   for (int i = 0; i < cellsz.size(); i++){
      for (int j = 0; j < cellsz[0].size(); j++){
         for (int k = 0; k < cellsz[0][0].size(); k++){

            Output.push_back(cellsz[i][j][k]);
            
         }
      }
   }

   return Output;
}

//--------------------------Private---------------------------------

void cells::binning(double cx, double cy, double cz, double sx, double sy, double sz){
   int cellSize = 100;

   int xbin, ybin, zbin;
   for (int i = 0; i < cellSize; i++)
   {
      if (cx >= cellxEdge[i] && cx < cellxEdge[i+1] )
         xbin=i;
      if (cy >= cellyEdge[i] && cy < cellyEdge[i+1] )
         ybin=i;
      if (cz >= cellzEdge[i] && cx < cellxEdge[i+1] )
         zbin=i;
   }

   cellsx[xbin][ybin][zbin] =+ (1.0/(bincount[xbin][ybin][zbin]))*(sx-cellsx[xbin][ybin][zbin]);
   cellsy[xbin][ybin][zbin] =+ (1.0/(bincount[xbin][ybin][zbin]))*(sy-cellsy[xbin][ybin][zbin]);
   cellsz[xbin][ybin][zbin] =+ (1.0/(bincount[xbin][ybin][zbin]))*(sz-cellsz[xbin][ybin][zbin]);

}
