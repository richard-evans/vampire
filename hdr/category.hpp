#include <vector>

namespace cat{

		class category_t {
		public:
			
			// list of atoms in each category
			std::vector <int> category_atoms;

			// magnetisation variables
			double magm;
			double mean_magm;
			double mx;
			double my;
			double mz;
			
			// energy and toque variables
			double torque;
			double mean_torque;
			double energy;
			
			// constraint variables
			double theta;
			double phi;
			
			// member functions
			category_t(); // constructor
			
		};
		
		extern std::vector <category_t> category;


} // End of namespace cat