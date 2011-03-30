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
	extern int two_temperature_pulse();

	// Sundry programs and diagnostics not under general release
	extern int LLB_Boltzmann();
	extern int timestep_scaling();
	
}

#endif /*PROGRAM_H_*/
