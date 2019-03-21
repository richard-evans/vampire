//----------------------------------------------------------------------------------
// Slonczewski spin torque field
//----------------------------------------------------------------------------------

// save polarization to temporary constant
const double stpx = slonczewski_spin_polarization_unit_vector[0];
const double stpy = slonczewski_spin_polarization_unit_vector[1];
const double stpz = slonczewski_spin_polarization_unit_vector[2];

const double staj = slonczewski_aj[material];
const double stbj = slonczewski_bj[material];

// calculate field
hx += staj*(sy*stpz - sz*stpy) + stbj*stpx;
hy += staj*(sz*stpx - sx*stpz) + stbj*stpy;
hz += staj*(sx*stpy - sy*stpx) + stbj*stpz;

// save field to spin field array
atoms::x_total_spin_field_array[atom]+=hx;
atoms::y_total_spin_field_array[atom]+=hy;
atoms::z_total_spin_field_array[atom]+=hz;
