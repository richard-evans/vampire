//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "create.hpp"
#include "vmath.hpp"
#include "voronoi.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace create{
namespace internal{

void voronoi_grain_rounding(std::vector <std::vector <double> > & grain_coord_array,
                            std::vector <std::vector <std::vector <double> > > &  grain_vertices_array){

	const int max_vertices=50;
   double tmp_grain_pointx_array[max_vertices];
	double tmp_grain_pointy_array[max_vertices];

	// calculate grain rounding
	for(unsigned int grain=0;grain<grain_coord_array.size();grain++){

		// get number of vertices for each grain
		const int nv = grain_vertices_array[grain].size();

		// Exclude grains with zero vertices
		if(nv!=0){

			// allocate temporary array for area calculation
			std::vector<std::vector <double> > rnd;
			rnd.resize(48);
			for(int idx=0; idx<48;idx++){
				rnd[idx].resize(2,0.0);
			}

			double area_frac=create_voronoi::area_cutoff;
			double deltar=0.5; // Angstroms
			double radius=0.0; // Strating radius
			double area=0.0; // Starting area

			// Set temporary vertex coordinates
			int num_vertices = grain_vertices_array[grain].size();
				for(int vertex=0;vertex<num_vertices;vertex++){
					tmp_grain_pointx_array[vertex]=grain_vertices_array[grain][vertex][0]; //-grain_coord_array[grain][0];
					tmp_grain_pointy_array[vertex]=grain_vertices_array[grain][vertex][1]; //-grain_coord_array[grain][1];
				}

			// calculate voronoi area
			double varea=0.0;
			for(int r=0;r<1000;r++){
				radius+=deltar;

				for(int i=0;i<48;i++){
					double theta = 2.0*M_PI*double(i)/48.0;
					double x = radius*cos(theta);
					double y = radius*sin(theta);

					// Check to see if site is within polygon
					if(vmath::point_in_polygon(x,y,tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)==true){
						rnd[i][0]=x;
						rnd[i][1]=y;
					}
				}

				//update area
				varea=0.0;
				for(int i=0;i<48;i++){
					int nvi = i+1;
					if(nvi>=48) nvi=0;
					varea+=0.5*sqrt((rnd[nvi][0]-rnd[i][0])*(rnd[nvi][0]-rnd[i][0]))*sqrt((rnd[nvi][1]-rnd[i][1])*(rnd[nvi][1]-rnd[i][1]));
				}
			}

			// reset polygon positions and radius
			for(int idx=0; idx<48;idx++){
				rnd[idx][0]=0.0;
				rnd[idx][1]=0.0;
			}
			radius=0.0;
			// expand polygon
			for(int r=0;r<100;r++){
			if(area<area_frac*varea){
				radius+=deltar;

				//loop over coordinates
				for(int i=0;i<48;i++){
					double theta = 2.0*M_PI*double(i)/48.0;
					double x = radius*cos(theta);
					double y = radius*sin(theta);

					// Check to see if site is within polygon
					if(vmath::point_in_polygon(x,y,tmp_grain_pointx_array,tmp_grain_pointy_array,num_vertices)==true){
						rnd[i][0]=x;
						rnd[i][1]=y;
					}
				}

				//update area
				area=0.0;
				for(int i=0;i<48;i++){
					int nvi = i+1;
					if(nvi>=48) nvi=0;
					area+=0.5*sqrt((rnd[nvi][0]-rnd[i][0])*(rnd[nvi][0]-rnd[i][0]))*sqrt((rnd[nvi][1]-rnd[i][1])*(rnd[nvi][1]-rnd[i][1]));
				}
			}
		}

		// set new polygon points
		grain_vertices_array[grain]=rnd;

		} // end of check for nv=0

	} // end of grain loop

   return;

}

} // end of namespace internal
} // end of namespace create
