/*
void calculate_spherical_harmonic_fields(const int start_index,const int end_index){

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
      const double ex = mp::material[imaterial].UniaxialAnisotropyUnitVector[0];
      const double ey = mp::material[imaterial].UniaxialAnisotropyUnitVector[1];
      const double ez = mp::material[imaterial].UniaxialAnisotropyUnitVector[2];
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

}*/






///--------------------------------------------------------------------------------------------------------------
///  Function to calculate spherical harmonic anisotropy energy
///
///  (c) R F L Evans 2015
///
///  Higher order anisotropies generally need to be described using spherical harmonics. The usual form (a
///  series in S leads to cross pollution of terms, giving strange temperature dependencies.
///
///  The harmonics are described with Legendre polynomials with even order, which for 2nd, 4th and 6th are:
///  ( http://en.wikipedia.org/wiki/Legendre_polynomials )
///
///  k_2(sz) = (1/2) *(3sz^2 - 1)
///  k_4(sz) = (1/8) *(35sz^4 - 30sz^2 + 3)
///  k_6(sz) = (1/16)*(231sz^6 - 315*sz^4 + 105sz^2 - 5)
///
///  The harmonics feature an arbritrary -2/3 factor compared with the usual form, and so in VAMPIRE these are
///  renormalised to maintain consistency for the 2nd order terms. This can be projected onto
///  any arbritrary direction ex,ey,ez allowing higher order anisotropy terms along any direction. This
///  direction is shared with the other uniaxial anisotropy coefficients since they should not be used
///  simultaneously.
///
///--------------------------------------------------------------------------------------------------------------
double spin_spherical_harmonic_aniostropy_energy(const int imaterial, const double sx, const double sy, const double sz){

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
   const double ex = mp::material[imaterial].UniaxialAnisotropyUnitVector[0];
   const double ey = mp::material[imaterial].UniaxialAnisotropyUnitVector[1];
   const double ez = mp::material[imaterial].UniaxialAnisotropyUnitVector[2];

   const double sdote2 = (sx*ex + sy*ey + sz*ez)*(sx*ex + sy*ey + sz*ez);
   const double sdote4 = sdote2*sdote2;
   const double sdote6 = sdote4*sdote2;

   // calculate field (double negative from scale factor and negative derivative)
   const double energy = scale*(k2*oneo2*(3.0*sdote2 - 1.0) +
                                k4*oneo8*(35.0*sdote4 - 30.0*sdote2 + 3.0) +
                                k6*oneo16*(231.0*sdote6 - 315.0*sdote4 + 105.0*sdote2 - 5.0));

   return energy;

}
