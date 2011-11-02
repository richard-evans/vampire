#ifndef PROGRAM_H_
#define PROGRAM_H_

//==========================================================
// Namespace program
//==========================================================
namespace program
{
	// program functions
	extern int bmark();
	extern void time_series();
	extern int hysteresis();
	extern int static_hysteresis();
	extern int curie_temperature();
	extern void field_cool();
	extern void two_temperature_pulse();
	extern void hamr();
	extern void cmc_anisotropy();
	extern void hybrid_cmc();

	// Sundry programs and diagnostics not under general release
	extern int LLB_Boltzmann();
	extern int timestep_scaling();
	extern void boltzmann_dist();
	
}

#endif /*PROGRAM_H_*/
