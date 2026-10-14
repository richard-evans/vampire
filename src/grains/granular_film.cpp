//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <random>
#include <cmath>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>

// Vampire headers
#include "create.hpp"
#include "errors.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "random.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vmath.hpp"
#include <algorithm>

// grains module headers
#include "internal.hpp"

namespace grains{

int voronoi_film(std::vector<cs::catom_t> & catom_array){

	// check calling of routine if error checking is activated
	if(err::check==true){
      terminaltextcolor(RED);
      std::cerr << "grains::voronoi_film has been called" << std::endl;
      terminaltextcolor(WHITE);
   }

	// grain_coord_array[g] = (x,y) of grain g's site, grain_vertices_array[g] = its polygon's vertex list
	std::vector <std::vector <double> > grain_coord_array;
	std::vector <std::vector <std::vector <double> > > grain_vertices_array;

	// Generate the grain seeds: sample target diameters, choose how many
	// grains are needed, place sites, fit per-seed weights w_i = r_i^2 for
	// the tessellation below. create:grain-size falls back to
	// dimensions:particle-size when unset.
	const double effective_grain_size = (internal::grain_size > 0.0) ?
		internal::grain_size : cs::particle_scale;

	internal::grain_size_distribution_params_t size_params;
	size_params.distribution         = internal::grain_size_distribution;
	size_params.mean_diameter        = effective_grain_size;
	size_params.sd                   = internal::grain_size_sd;
	size_params.second_mean_diameter = internal::grain_size_second;
	size_params.second_sd            = internal::grain_size_second_sd;
	size_params.second_fraction      = internal::grain_size_second_fraction;
	size_params.regularity           = internal::grain_regularity;
	size_params.periodic             = internal::grain_periodic_boundaries;
	if(internal::grain_size_distribution == internal::grain_size_file){
		size_params.file_diameters = internal::read_grain_size_distribution_file(internal::grain_size_distribution_file);
	}

	// Seed the grains module's own random generator so grain-structure
	// generation is reproducible independently of the global random stream.
	internal::grnd.seed(internal::grain_structure_seed);

	const internal::grain_seed_set_t seeds = internal::generate_grain_seeds(
		cs::system_dimensions[0], cs::system_dimensions[1], size_params, internal::grnd);

	const size_t num_seeds = seeds.sites.size();
	grain_coord_array.resize(num_seeds);
	grain_vertices_array.resize(num_seeds);
	for(size_t s=0; s<num_seeds; s++){
		grain_coord_array[s].resize(2);
		grain_coord_array[s][0] = seeds.sites[s].x;
		grain_coord_array[s][1] = seeds.sites[s].y;
	}

	int grain = int(num_seeds);

	//-----------------------
	// Check for grains >=1
	//-----------------------
	if(grain<1){
		terminaltextcolor(RED);
		std::cerr << "Error! - No grains found in structure - Increase system dimensions" << std::endl;
		terminaltextcolor(WHITE);
		zlog << zTs() << "Error! - No grains found in structure - Increase system dimensions" << std::endl;
		err::vexit();
	}

	// create:grain-tessellation selects plain Voronoi vs weighted Laguerre
	internal::populate_vertex_points_power(grain_coord_array, grain_vertices_array, seeds.weights,
		cs::system_dimensions[0], cs::system_dimensions[1],
		internal::grain_tessellation == internal::tessellation_laguerre,
		internal::grain_periodic_boundaries);

	// non-periodic films only: remove every grain cut by the film edge
	const int num_boundary_removed = internal::grain_periodic_boundaries ? 0 :
		internal::remove_boundary_grains(grain_coord_array, grain_vertices_array,
			0.0, 0.0, cs::system_dimensions[0], cs::system_dimensions[1]);

	// grain-size statistics, measured before spacing/rounding/faceting shrink the cells
	if(internal::output_grain_statistics_file && vmpi::master){
		internal::write_grain_statistics_file("grain_statistics.txt",
			grain_coord_array, grain_vertices_array, seeds.population, seeds.target_diameter,
			0.0, 0.0, cs::system_dimensions[0], cs::system_dimensions[1]);
	}

	{
		const internal::grain_statistics_summary_t stats =
			internal::compute_grain_size_statistics(grain_coord_array, grain_vertices_array, seeds.population);

		const bool is_bimodal = (internal::grain_size_distribution == internal::grain_size_bimodal);
		const int placed = stats.placed_primary + stats.placed_second;

		std::ostringstream msg;
		msg << "Granular film: requested mean diameter = " << effective_grain_size/10.0 << " nm";
		if(is_bimodal){
			msg << ", second population = " << internal::grain_size_second/10.0 <<
				" nm, requested number fraction = " << internal::grain_size_second_fraction;
		}
		msg << std::endl;
		msg << "  Placed " << placed << " grains";
		if(is_bimodal){
			msg << " (" << stats.placed_primary << " primary, " << stats.placed_second << " second population)";
		}
		// exclude boundary-removed cells to leave only the weight-fitting count
		msg << ", " << (stats.num_vanished - num_boundary_removed) << " vanished during weight fitting." << std::endl;
		if(num_boundary_removed > 0){
			msg << "  Removed " << num_boundary_removed << " grain(s) touching the film boundary "
				"(create:grain-periodic-boundaries = false), leaving a void margin at the film edges." << std::endl;
		}
		msg << "  Measured mean equivalent diameter: primary = " << stats.mean_diameter_primary/10.0 <<
			" nm (sd " << stats.sd_diameter_primary/10.0 << " nm)";
		if(is_bimodal){
			msg << ", second = " << stats.mean_diameter_second/10.0 << " nm (sd " << stats.sd_diameter_second/10.0 << " nm)";
		}
		msg << std::endl;
		if(is_bimodal){
			msg << "  Realised number fraction (second population) = " << stats.realised_number_fraction <<
				", realised area fraction = " << stats.realised_area_fraction << std::endl;
		}
		if(!seeds.pack.converged){
			msg << "  Warning: grain seed packing did not fully relax (jammed with max overlap fraction " <<
				seeds.pack.achieved_max_overlap << " after " << seeds.pack.iterations << " iterations)." << std::endl;
		}
		if(!seeds.fit.converged){
			msg << "  Warning: grain weight fitting did not fully converge (achieved tolerance " <<
				seeds.fit.achieved_tol << " after " << seeds.fit.iterations << " iterations)." << std::endl;
		}
		// warn if a grain-size-sd/second-* keyword is set but ignored (distribution == delta)
		if(internal::grain_size_distribution == internal::grain_size_delta &&
			internal::grain_size_shape_param_explicitly_set){
			msg << "  Warning: create:grain-size-distribution = delta (monodisperse), so grain-size-sd and any "
				"grain-size-second-* settings have no effect on the realised grain sizes." << std::endl;
		}

		std::cout << msg.str();
		zlog << zTs() << msg.str();
	}

	// Wulff faceting: clip each tessellated cell toward a faceted equilibrium
	// shape, before the spacing offset and rounding below.
	const std::vector<double> grain_facet_orientation_array = internal::generate_grain_facet_orientations(
		int(grain_coord_array.size()), internal::grain_facet_symmetry,
		internal::grain_facet_orientation, internal::grain_facet_orientation_angle,
		internal::grain_facet_orientation_spread, internal::grnd);

	// record each grain's facet orientation for later per-grain property lookup
	internal::set_raw_orientation(grain_facet_orientation_array);

	internal::apply_grain_faceting(grain_coord_array, grain_vertices_array, grain_facet_orientation_array,
		internal::grain_facet_symmetry, internal::grain_facet_anisotropy, internal::grain_facet_strength);

	// pull each grain boundary inward by half the requested spacing
	const double effective_grain_spacing = (internal::grain_spacing >= 0.0) ?
		internal::grain_spacing : cs::particle_spacing;

	// create:grain-spacing-sd jitters each grain's own half-gap independently
	const std::vector<double> grain_spacing_delta = internal::generate_grain_spacing_jitter(
		int(grain_coord_array.size()), 0.5*effective_grain_spacing, internal::grain_spacing_sd, internal::grnd);
	internal::apply_grain_spacing(grain_coord_array, grain_vertices_array, grain_spacing_delta);

   // round grain corners if requested (create:grain-rounding; 0.0 = off)
	if(internal::grain_rounding > 0.0){
		// create:grain-boundary-roughness switches to the organic-rounding overload
		if(internal::grain_boundary_roughness > 0.0){
			internal::voronoi_grain_rounding(grain_coord_array, grain_vertices_array, internal::grain_rounding,
				internal::grain_boundary_roughness, internal::grain_boundary_roughness_modes, internal::grnd);
		}
		else{
			internal::voronoi_grain_rounding(grain_coord_array, grain_vertices_array, internal::grain_rounding);
		}
	}

   // Bin atoms into supercells once, shared by every grain's atom-assignment pass
   const internal::supercell_array_t supercell_array =
      internal::build_supercell_array(catom_array, cs::unit_cell.dimensions[0], cs::unit_cell.dimensions[1],
                                               cs::total_num_unit_cells[0], cs::total_num_unit_cells[1]);

   // Determine order for core-shell grains
   std::list<create::core_radius_t> material_order = create::sorted_core_shell_materials();

	// determine whether elliptical grain rounding is active
	const bool elliptical_rounding = internal::voronoi_elliptical_rounding > 0.0;

	std::cout <<"Generating Voronoi Grains" << std::flush;
	zlog << zTs() << "Generating Voronoi Grains" << std::flush;
	if(elliptical_rounding){
		zlog << zTs() << "Applying elliptical grain rounding with rounding factor " <<
			internal::voronoi_elliptical_rounding << " about height fraction " <<
			internal::voronoi_elliptical_rounding_height << std::endl;
	}

	// optionally output grain vertices to file
	std::ofstream gvfile;
	if(internal::output_gv_file && vmpi::master){
		gvfile.open("grain_shapes.txt");
	}

	// loop over all grains with vertices
	for(unsigned int grain=0;grain<grain_coord_array.size();grain++){
		// Exclude grains with zero vertices

		internal::print_grain_progress(grain, grain_coord_array.size());
		if(grain_vertices_array[grain].size()!=0){

			if(internal::output_gv_file){
				const double dx = grain_coord_array[grain][0];
				const double dy = grain_coord_array[grain][1];
				internal::write_grain_vertices(grain, dx, dy, gvfile, grain_vertices_array[grain]);
			}

         // determine coordinate offset for grains
         const double x0 = grain_coord_array[grain][0];
         const double y0 = grain_coord_array[grain][1];

         internal::grain_footprint_t fp =
            internal::compute_grain_footprint(grain_vertices_array[grain], x0, y0,
                                                       cs::unit_cell.dimensions[0], cs::unit_cell.dimensions[1],
                                                       internal::grain_periodic_boundaries,
                                                       cs::total_num_unit_cells[0], cs::total_num_unit_cells[1]);
         const int num_vertices = fp.px.size();

         internal::assign_atoms_in_footprint(supercell_array, fp, [&](int atom){

            // Get atomic position
            double x = catom_array[atom].x;
            double y = catom_array[atom].y;

            // compute in-plane scale factor for elliptical grain rounding
            const double erf = internal::elliptical_rounding_factor(catom_array[atom].z);

            // skip atoms where the grain cross-section has closed completely
            if(erf <= 0.0) return;

            // atom position relative to grain centre, minimum-image wrapped if periodic
            double dx = x-x0;
            double dy = y-y0;
            if(internal::grain_periodic_boundaries){
               dx = internal::wrap_minimum_image(dx, cs::system_dimensions[0]);
               dy = internal::wrap_minimum_image(dy, cs::system_dimensions[1]);
            }

            if(mp::material[catom_array[atom].material].core_shell_size>0.0){
               // core-shell grains: assign largest shell first, smallest core last (overwrites)
               for(std::list<create::core_radius_t>::iterator it = material_order.begin(); it !=  material_order.end(); it++){
                  int mat = (it)->mat;
                  double factor = mp::material[mat].core_shell_size;
                  double maxz=create::get_material_height_max(mat)*cs::system_dimensions[2];
                  double minz=create::get_material_height_min(mat)*cs::system_dimensions[2];
                  double cz=catom_array[atom].z;
                  const int atom_uc_cat = catom_array[atom].uc_category;
                  const int mat_uc_cat = create::get_material_unit_cell_category(mat);
                  // check for within core shell range
                  const bool in_grain = vmath::point_in_polygon_scaled(dx, dy, factor*erf, fp.px.data(), fp.py.data(), num_vertices);
                  if(in_grain){
                     if((cz>=minz) && (cz<maxz) && (atom_uc_cat == mat_uc_cat) ){
                        catom_array[atom].include=true;
                        catom_array[atom].material=mat;
                        catom_array[atom].grain=grain;
                     }
                     // if set to clear atoms then remove atoms within radius
                     else if(cs::fill_core_shell==false){
                        catom_array[atom].include=false;
                     }
                  }
               }
            }
            // Check to see if site is within the (optionally elliptically rounded) grain
            else if(vmath::point_in_polygon_scaled(dx,dy,erf,fp.px.data(),fp.py.data(),num_vertices)==true){
               catom_array[atom].include=true;
               catom_array[atom].grain=grain;
            }

         });
		}
	}
	terminaltextcolor(GREEN);
	std::cout << "done!" << std::endl;
	terminaltextcolor(WHITE);
	zlog << "done!" << std::endl;

	// add final grain for continuous layer
	grain_coord_array.push_back(std::vector <double>());
	grain_coord_array[grain_coord_array.size()-1].push_back(0.0); // x
	grain_coord_array[grain_coord_array.size()-1].push_back(0.0); // y

	// check for continuous layer
	for(unsigned int atom=0; atom < catom_array.size(); atom++){
	  if(mp::material[catom_array[atom].material].continuous==true && catom_array[atom].include == false ){
	    catom_array[atom].include=true;
	    catom_array[atom].grain=int(grain_coord_array.size()-1);
	  }
	}

	// set number of grains
	grains::num_grains = int(grain_coord_array.size());

	// sort atoms by grain number
	create::sort_atoms_by_grain(catom_array);

	return EXIT_SUCCESS;
}

} // End of grains namespace
