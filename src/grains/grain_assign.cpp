//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>

// grains module headers
#include "internal.hpp"

namespace grains{
   namespace internal{

      //--------------------------------------------------------------------
      // Bins every atom into a 2D array of supercells indexed by
      // (x/unit_cell_dim_x, y/unit_cell_dim_y), so that a grain only ever
      // needs to scan the atoms in the supercells its bounding box covers
      // rather than every atom in the system.
      //--------------------------------------------------------------------
      supercell_array_t build_supercell_array(const std::vector<cs::catom_t>& catom_array,
                                               double unit_cell_dim_x, double unit_cell_dim_y,
                                               int num_unit_cells_x, int num_unit_cells_y){

         supercell_array_t supercell_array;
         supercell_array.resize(num_unit_cells_x);
         for(int i=0; i<num_unit_cells_x; i++) supercell_array[i].resize(num_unit_cells_y);

         for(unsigned int atom=0; atom<catom_array.size(); atom++){
            const int cx = int(catom_array[atom].x/unit_cell_dim_x);
            const int cy = int(catom_array[atom].y/unit_cell_dim_y);
            supercell_array.at(cx).at(cy).push_back(atom);
         }

         return supercell_array;

      }

      //--------------------------------------------------------------------
      // Computes the range of supercells a grain's vertex polygon touches,
      // and flattens the vertex coordinates (relative to the grain centre
      // x0,y0) to the double* layout vmath::point_in_polygon*() expects.
      //
      // Under periodicity the raw range can span more supercells than the
      // domain has (an unwrapped, boundary-straddling cell), so it is capped
      // to at most one full domain width/height.
      //--------------------------------------------------------------------
      grain_footprint_t compute_grain_footprint(const std::vector<std::vector<double> >& vertices,
                                                 double x0, double y0,
                                                 double unit_cell_dim_x, double unit_cell_dim_y,
                                                 bool periodic,
                                                 int num_unit_cells_x,
                                                 int num_unit_cells_y){

         grain_footprint_t fp;
         fp.minx = 10000000;
         fp.maxx = 0;
         fp.miny = 10000000;
         fp.maxy = 0;
         fp.periodic = periodic;
         fp.num_unit_cells_x = num_unit_cells_x;
         fp.num_unit_cells_y = num_unit_cells_y;

         const int num_vertices = vertices.size();
         fp.px.resize(num_vertices);
         fp.py.resize(num_vertices);

         for(int vertex=0; vertex<num_vertices; vertex++){
            fp.px[vertex] = vertices[vertex][0];
            fp.py[vertex] = vertices[vertex][1];
            const int x = int((fp.px[vertex]+x0)/unit_cell_dim_x);
            const int y = int((fp.py[vertex]+y0)/unit_cell_dim_y);
            if(x < fp.minx) fp.minx = x;
            if(x > fp.maxx) fp.maxx = x;
            if(y < fp.miny) fp.miny = y;
            if(y > fp.maxy) fp.maxy = y;
         }

         if(periodic){
            if(num_unit_cells_x > 0 && fp.maxx - fp.minx + 1 > num_unit_cells_x) fp.maxx = fp.minx + num_unit_cells_x - 1;
            if(num_unit_cells_y > 0 && fp.maxy - fp.miny + 1 > num_unit_cells_y) fp.maxy = fp.miny + num_unit_cells_y - 1;
         }

         return fp;

      }

      //--------------------------------------------------------------------
      // Calls visit(atom_index) for every atom binned into a supercell
      // within fp's bounding box. Under periodicity (fp.periodic) indices
      // are wrapped modulo fp.num_unit_cells_x/y, since fp.minx/maxx/miny/maxy
      // may be negative or exceed the domain for a cell built in unwrapped
      // coordinates.
      //--------------------------------------------------------------------
      void assign_atoms_in_footprint(const supercell_array_t& supercell_array,
                                      const grain_footprint_t& fp,
                                      const std::function<void(int atom)>& visit){

         for(int i=fp.minx; i<=fp.maxx; i++){
            const int wi = fp.periodic ? (((i % fp.num_unit_cells_x) + fp.num_unit_cells_x) % fp.num_unit_cells_x) : i;
            for(int j=fp.miny; j<=fp.maxy; j++){
               const int wj = fp.periodic ? (((j % fp.num_unit_cells_y) + fp.num_unit_cells_y) % fp.num_unit_cells_y) : j;
               for(unsigned int id=0; id<supercell_array[wi][wj].size(); id++){
                  visit(supercell_array[wi][wj][id]);
               }
            }
         }

         return;

      }

      //--------------------------------------------------------------------
      // See internal.hpp.
      //--------------------------------------------------------------------
      double wrap_minimum_image(double d, double domain){

         if(!(domain > 0.0)) return d;
         return d - domain*std::round(d/domain);

      }

   } // end of internal namespace
} // end of grains namespace
