//----------------------------------------------------------------------------
// Optionally calculate random atomic anisotropy axes for spherical harmonics
//----------------------------------------------------------------------------
if(sim::spherical_harmonics && sim::random_anisotropy){
   // Resize arrays
   atoms::uniaxial_anisotropy_vector_x.resize(atoms::num_atoms,0.0);
   atoms::uniaxial_anisotropy_vector_y.resize(atoms::num_atoms,0.0);
   atoms::uniaxial_anisotropy_vector_z.resize(atoms::num_atoms,0.0);

   std::vector<double> grain_anisotropy_directions(0);
   // check for grain level random anisotropy
   if(grains::random_anisotropy){
      // resize array storing grain anisotropy vectors
      grain_anisotropy_directions.resize(3*grains::num_grains);

      // calculate anisotropy directions for all grains on root process
      if(vmpi::my_rank == 0){
         for(int g=0; g<grains::num_grains; g++){

            double x = mtrandom::gaussian();
            double y = mtrandom::gaussian();
            double z = mtrandom::gaussian();

            // Calculate vector length
            const double r = 1.0/sqrt (x*x + y*y + z*z);

            grain_anisotropy_directions[3*g + 0] = x*r;
            grain_anisotropy_directions[3*g + 1] = y*r;
            grain_anisotropy_directions[3*g + 2] = z*r;

         }
      }

      #ifdef MPICF
         // Broadcast calculated anisotropy directions to all nodes
         MPI_Bcast(&grain_anisotropy_directions[0], grain_anisotropy_directions.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      #endif

   }

   // Unroll anisotropy directions
   for(int atom = 0; atom < atoms::num_atoms; ++atom){

      // Determine material number
      int imaterial = atoms::type_array[atom];

      // Calculate random anisotropy directions on unit sphere
      if(mp::material[imaterial].random_anisotropy){

         double x = mtrandom::gaussian();
         double y = mtrandom::gaussian();
         double z = mtrandom::gaussian();

         // Calculate vector length
         const double r = 1.0/sqrt (x*x + y*y + z*z);

         // Save direction
         atoms::uniaxial_anisotropy_vector_x[atom] = x*r;
         atoms::uniaxial_anisotropy_vector_y[atom] = y*r;
         atoms::uniaxial_anisotropy_vector_z[atom] = z*r;

      }
      // If random grain anisotropy defined set local anisotropy to the grain anisotropy direction
      else if(mp::material[imaterial].random_grain_anisotropy){

         const int grain = atoms::grain_array[atom];

         atoms::uniaxial_anisotropy_vector_x[atom] = grain_anisotropy_directions[3*grain + 0];
         atoms::uniaxial_anisotropy_vector_y[atom] = grain_anisotropy_directions[3*grain + 1];
         atoms::uniaxial_anisotropy_vector_z[atom] = grain_anisotropy_directions[3*grain + 2];

      }
      // Otherwise unroll anisotropy directions
      else{
         atoms::uniaxial_anisotropy_vector_x[atom] = mp::material[imaterial].UniaxialAnisotropyUnitVector[0];
         atoms::uniaxial_anisotropy_vector_y[atom] = mp::material[imaterial].UniaxialAnisotropyUnitVector[1];
         atoms::uniaxial_anisotropy_vector_z[atom] = mp::material[imaterial].UniaxialAnisotropyUnitVector[2];
      }

   }
}

///--------------------------------------------------------------------------------------------------------------
///  Function to calculate random spherical harmonic anisotropy fields
///
///  (c) R F L Evans 2015
///
///  In this function uniaxial anisotropy is calculated using spherical harmonics,
///  except each atom is allowed a locally defined anisotropy axis. This comes with
///  a performance cost, and so this version is only caled if needed (defined by the
///  sim::random_anisotropy flag).
///
///--------------------------------------------------------------------------------------------------------------
void calculate_random_spherical_harmonic_fields(const int start_index,const int end_index){

  // rescaling prefactor
  const double scale = 2.0/3.0; // Factor to rescale anisotropies to usual scale

  // constant factors
  const double oneo8 = 1.0/8.0;
  const double oneo16 = 1.0/16.0;

  // loop over all atoms
  for(int atom=start_index; atom<end_index; atom++){

    // Determine atom type
    const int imaterial=atoms::type_array[atom];

    // determine harmonic constants for material
    const double k2 = mp::material_spherical_harmonic_constants_array[3*imaterial + 0];
    const double k4 = mp::material_spherical_harmonic_constants_array[3*imaterial + 1];
    const double k6 = mp::material_spherical_harmonic_constants_array[3*imaterial + 2];

    // determine anisotropy direction and dot product
    const double ex = atoms::uniaxial_anisotropy_vector_x[atom];
    const double ey = atoms::uniaxial_anisotropy_vector_y[atom];
    const double ez = atoms::uniaxial_anisotropy_vector_z[atom];
    const double sx = atoms::x_spin_array[atom];
    const double sy = atoms::y_spin_array[atom];
    const double sz = atoms::z_spin_array[atom];

    const double sdote = (sx*ex + sy*ey + sz*ez);
    const double sdote3 = sdote*sdote*sdote;
    const double sdote5 = sdote3*sdote*sdote;

    // calculate field (double negative from scale factor and negative derivative)
    atoms::x_total_spin_field_array[atom] += scale*ex*(k2*3.0*sdote + k4*oneo8*(140.0*sdote3 - 60.0*sdote) + k6*oneo16*(1386.0*sdote5 - 1260.0*sdote3 + 210.0*sdote));
    atoms::y_total_spin_field_array[atom] += scale*ey*(k2*3.0*sdote + k4*oneo8*(140.0*sdote3 - 60.0*sdote) + k6*oneo16*(1386.0*sdote5 - 1260.0*sdote3 + 210.0*sdote));
    atoms::z_total_spin_field_array[atom] += scale*ez*(k2*3.0*sdote + k4*oneo8*(140.0*sdote3 - 60.0*sdote) + k6*oneo16*(1386.0*sdote5 - 1260.0*sdote3 + 210.0*sdote));

  }

  return;

}
///--------------------------------------------------------------------------------------------------------------
///  Function to calculate spherical harmonic anisotropy energy
///
///  (c) R F L Evans 2015
///
///  In this function uniaxial anisotropy is calculated using spherical harmonics,
///  except each atom is allowed a locally defined anisotropy axis. This comes with
///  a performance cost, and so this version is only called if needed (defined by the
///  sim::random_anisotropy flag).
///
///--------------------------------------------------------------------------------------------------------------
double spin_spherical_harmonic_random_aniostropy_energy(const int atom, const int imaterial, const double sx, const double sy, const double sz){

   // rescaling prefactor
   const double scale = -2.0/3.0; // Factor to rescale anisotropies to usual scale

   // constant factors
   const double oneo2 = 0.5;
   const double oneo8 = 1.0/8.0;
   const double oneo16 = 1.0/16.0;

   // determine harmonic constants for material
   const double k2 = mp::material_spherical_harmonic_constants_array[3*imaterial + 0];
   const double k4 = mp::material_spherical_harmonic_constants_array[3*imaterial + 1];
   const double k6 = mp::material_spherical_harmonic_constants_array[3*imaterial + 2];

   // determine anisotropy direction and dot product
	const double ex = atoms::uniaxial_anisotropy_vector_x[atom];
	const double ey = atoms::uniaxial_anisotropy_vector_y[atom];
	const double ez = atoms::uniaxial_anisotropy_vector_z[atom];

   const double sdote2 = (sx*ex + sy*ey + sz*ez)*(sx*ex + sy*ey + sz*ez);
   const double sdote4 = sdote2*sdote2;
   const double sdote6 = sdote4*sdote2;

   // calculate field (double negative from scale factor and negative derivative)
   const double energy = scale*(k2*oneo2*(3.0*sdote2 - 1.0) +
                                k4*oneo8*(35.0*sdote4 - 30.0*sdote2 + 3.0) +
                                k6*oneo16*(231.0*sdote6 - 315.0*sdote4 + 105.0*sdote2 - 5.0));

   return energy;

}
