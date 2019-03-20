// calculate resistance

// parallel resitsors (one from each stack)

R = R_P + R_AP
I = V/R

for(all cell in stack)

calculate mi
if(magnetic){
   calculate RP and RAP
   calculate mj
   R += R_P +( R_AP - R_P) (1 - m1 . m2) / 2
   RP = 0
   RAP = 0
   mi = mj
}
else{
   // calculate sum of resistances until next moment reached
   RP+=RP
   RAP+=RAP
}

if(end) R += (RP + RAP)/2 // assume non magnetic transport

I = V/R
calculate spin torque in stack -> assume lambda < cell size - each mi/mj pair form a polariser/free layer
spin torque propto mi x mj (no effect for parallel spins, effects at DW and MTJ layer)

STT_factor = mi x mj

//----------------------------------------------------------------------------------
// Slonczewski spin torque field
//----------------------------------------------------------------------------------

// save polarization to temporary constant
const double stpx = slonczewski_spin_polarization_unit_vector[0]; // = mi
const double stpy = slonczewski_spin_polarization_unit_vector[1];
const double stpz = slonczewski_spin_polarization_unit_vector[2];

const double staj = slonczewski_aj[material];
const double stbj = slonczewski_bj[material];

// calculate field
hx += staj*(sy*stpz - sz*stpy) + stbj*stpx;
hy += staj*(sz*stpx - sx*stpz) + stbj*stpy;
hz += staj*(sx*stpy - sy*stpx) + stbj*stpz;
