//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "create.hpp"
#include "vio.hpp"
#include "vmath.hpp"

// create module headers
#include "internal.hpp"

//------------------------------------------------------------------------------
//	   Function to decide if test point is within polygon defined by points
//------------------------------------------------------------------------------
bool create::internal::point_in_polygon(const create::internal::points_t test,             // pair of test points
													 std::vector<create::internal::points_t>& points){  // list of points defining a polygon

	// Add tiny amount to atomic coordinates to include atoms at 0,0,0
	test.x += 1e-10;
	test.y += 1e-10;

	const int poly_sides = points.size();
	int j = poly_sides - 1 ;
	bool odd_nodes=false;

	for ( int i = 0 ; i < poly_sides; i++) {
		if ( ( points[i].y < test.y && points[j].y >= test.y ) || ( points[j].y < test.y && points[i].y >= test.y) ) {
			if (points[i].x+(test.y-points[i].y) / ( points[j].y - points[i].y ) * ( points[j].x - points[i].x) < test.x) {
				odd_nodes = !odd_nodes;
			}
		}
		j=i;
	}

  return odd_nodes;

}
