//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016, 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "grains.hpp"

// grains module headers
#include "internal.hpp"

// Array-format adapter over round_polygon(): trims the sharp corners off
// each grain's polygon to area_fraction of its own area. Vertices are in
// the grain's own centre-relative frame, so centre is the origin here.
namespace grains{
namespace internal{

void voronoi_grain_rounding(std::vector <std::vector <double> > & grain_coord_array,
                            std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                            double area_fraction){

   for(size_t grain=0; grain<grain_vertices_array.size(); grain++){

      const size_t nv = grain_vertices_array[grain].size();
      if(nv == 0) continue;

      polygon_t cell(nv);
      for(size_t v=0; v<nv; v++){
         cell[v] = point2_t(grain_vertices_array[grain][v][0], grain_vertices_array[grain][v][1]);
      }
      cell = ensure_ccw(cell);

      const polygon_t rounded = round_polygon(cell, point2_t(0.0, 0.0), area_fraction);

      grain_vertices_array[grain].assign(rounded.size(), std::vector<double>(2));
      for(size_t v=0; v<rounded.size(); v++){
         grain_vertices_array[grain][v][0] = rounded[v].x;
         grain_vertices_array[grain][v][1] = rounded[v].y;
      }

   }

   return;

}

// Organic-rounding overload: each grain gets an independent random blob outline
// instead of a uniformly rounded polygon.
void voronoi_grain_rounding(std::vector <std::vector <double> > & grain_coord_array,
                            std::vector <std::vector <std::vector <double> > > &  grain_vertices_array,
                            double area_fraction,
                            double roughness_amplitude, int roughness_modes, MTRand& rng){

   for(size_t grain=0; grain<grain_vertices_array.size(); grain++){

      const size_t nv = grain_vertices_array[grain].size();
      if(nv == 0) continue;

      polygon_t cell(nv);
      for(size_t v=0; v<nv; v++){
         cell[v] = point2_t(grain_vertices_array[grain][v][0], grain_vertices_array[grain][v][1]);
      }
      cell = ensure_ccw(cell);

      const polygon_t rounded = round_polygon(cell, point2_t(0.0, 0.0), area_fraction,
                                               roughness_amplitude, roughness_modes, rng);

      grain_vertices_array[grain].assign(rounded.size(), std::vector<double>(2));
      for(size_t v=0; v<rounded.size(); v++){
         grain_vertices_array[grain][v][0] = rounded[v].x;
         grain_vertices_array[grain][v][1] = rounded[v].y;
      }

   }

   return;

}

} // end of namespace internal
} // end of namespace grains
