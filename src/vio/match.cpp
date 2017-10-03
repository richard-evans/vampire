//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and Rory Pond 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk and rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <sstream>

// Vampire headers
// Headers
#include "anisotropy.hpp"
#include "vio.hpp"
#include "sim.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "material.hpp"
#include "micromagnetic.hpp"
#include "environment.hpp"
#include "grains.hpp"
#include "stats.hpp"
#include "units.hpp"
#include "config.hpp"
#include "demag.hpp"
#include "cells.hpp"
#include "voronoi.hpp"
#include "ltmp.hpp"
#include "random.hpp"
#include "spintorque.hpp"
#include "unitcell.hpp"

// vio module headers
#include "internal.hpp"

namespace vin{

    /// @brief Function to match keywords, variables and units to an initialisation variable.
    ///
    /// @section License
    /// Use of this code, either in source or compiled form, is subject to license from the authors.
    /// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
    ///
    /// @section Information
    /// @author  Richard Evans, rfle500@york.ac.uk
    /// @version 1.1
    /// @date    18/01/2010
    ///
    /// @param[in] keyword Unique string variable linked to an initialisation variable
    /// @param[in] value Value of keyword linked to initialisation variable
    /// @return EXIT_SUCCESS
    ///
    /// @internal
    ///	Created:		15/01/2010
    ///	Revision:	  ---
    ///=====================================================================================
    ///
    int match(string const key, string const word, string const value, string const unit, int const line){

        std::string test;

        //-------------------------------------------------------------------
        // Call module input parameters
        //-------------------------------------------------------------------
        if(ltmp::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        else if(anisotropy::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        else if(cells::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        else if(create::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        else if(dipole::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        else if(exchange::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        else if(sim::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        else if(st::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        else if(unitcell::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        else if(micromagnetic::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        else if(environment::match_input_parameter(key, word, value, unit, line)) return EXIT_SUCCESS;
        //===================================================================
        // Test for create variables
        //===================================================================
        else
        test="create";
        if(key==test){
            int frs=vin::match_create(word, value, unit, line);
            return frs;
        }
        //===================================================================
        // Test for dimension variables
        //===================================================================
        else
        test="dimensions";
        if(key==test){
            int frs=vin::match_dimension(word, value, unit, line);
            return frs;
        }
        //===================================================================
        // Test for simulation variables
        //===================================================================
        else
        test="sim";
        if(key==test){
            int frs=vin::match_sim(word, value, unit, line);
            return frs;
        }
        //===================================================================
        // Test for data file output
        //===================================================================
        else
        test="output";
        if(key==test){
            int frs=vin::match_vout_list(word, value, line, vout::file_output_list);
            return frs;
        }
        //===================================================================
        // Test for screen output
        //===================================================================
        else
        test="screen";
        if(key==test){
            int frs=vin::match_vout_list(word, value, line, vout::screen_output_list);
            return frs;
        }
        //===================================================================
        // Test for grain output
        //===================================================================
        else
        test="grain";
        if(key==test){
            int frs=vin::match_vout_grain_list(word, value, line, vout::grain_output_list);
            return frs;
        }
        //===================================================================
        // Test for config output
        //===================================================================
        else
        test="config";
        if(key==test){
            int frs=config::match_input_parameter(word, value, unit, line);
            return frs;
        }
        //-------------------------------------------------------------------
        // Get material filename
        //-------------------------------------------------------------------
        else
        test="material";
        if(key==test){
            test="file";
            if(word==test){
                std::string matfile=value;
                // strip quotes
                matfile.erase(remove(matfile.begin(), matfile.end(), '\"'), matfile.end());
                test="";
                if(matfile!=test){
                    //std::cout << matfile << std::endl;
                read_mat_file(matfile,line);
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error - empty filename in control statement \'material:" << word << "\' on line " << line << " of input file" << std::endl;
                    terminaltextcolor(WHITE);
                    return EXIT_FAILURE;
                }
            }
        }
        else
            terminaltextcolor(RED);
            std::cerr << "Error - Unknown control statement \'" << key <<":"<< word << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
            return EXIT_FAILURE;

    } // end of match function

    int match_create(string const word, string const value, string const unit, int const line){
        ///-------------------------------------------------------------------
        /// system_creation_flags[1] - Set system particle shape
        ///-------------------------------------------------------------------

        std::string prefix="create:";

        // cs::system_creation_flags needs refactoring for readability and bug resistance
        std::string test="full";
        if(word==test){
            cs::system_creation_flags[1]=0;
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="cube";
        if(word==test){
            cs::system_creation_flags[1]=1;
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="cylinder";
        if(word==test){
            cs::system_creation_flags[1]=2;
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="ellipsoid";
        if(word==test){
            cs::system_creation_flags[1]=3;
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="sphere";
        if(word==test){
            cs::system_creation_flags[1]=4;
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="truncated-octahedron";
        if(word==test){
            cs::system_creation_flags[1]=5;
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="tear-drop";
        if(word==test){
            cs::system_creation_flags[1]=6;
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        // system_creation_flags[2] - Set system type
        //-------------------------------------------------------------------
        test="particle";
        if(word==test){
            cs::system_creation_flags[2]=0;
            return EXIT_SUCCESS;
        }
        else
        test="particle-array";
        if(word==test){
            cs::system_creation_flags[2]=1;
            return EXIT_SUCCESS;
        }
        else
        test="hexagonal-particle-array";
        if(word==test){
            cs::system_creation_flags[2]=2;
            return EXIT_SUCCESS;
        }
        else
        test="voronoi-film";
        if(word==test){
            cs::system_creation_flags[2]=3;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        else
        test="voronoi-size-variance";
        if(word==test){
            double vsd=atof(value.c_str());
            check_for_valid_value(vsd, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
            create_voronoi::voronoi_sd=vsd;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        else
        test="voronoi-row-offset";
        if(word==test){
            create_voronoi::parity=1;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        else
        test="voronoi-random-seed";
        if(word==test){
            int vs=atoi(value.c_str());
            check_for_valid_int(vs, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
            mtrandom::voronoi_seed=vs;
                return EXIT_SUCCESS;
        }
        else
        test="voronoi-rounded-grains";
        if(word==test){
            create_voronoi::rounded=true;
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="voronoi-rounded-grains-area";
        if(word==test){
            double vsd=atof(value.c_str());
            check_for_valid_value(vsd, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
            create_voronoi::area_cutoff=vsd;
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="particle-centre-offset"; //parity
        if(word==test){
            cs::particle_creation_parity=1;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        else
        test="single-spin";
        if(word==test){
            cs::single_spin=true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        else
        test="periodic-boundaries-x";
        if(word==test){
            cs::pbc[0]=true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        else
        test="periodic-boundaries-y";
        if(word==test){
            cs::pbc[1]=true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        else
        test="periodic-boundaries-z";
        if(word==test){
            cs::pbc[2]=true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        else
        test="select-material-by-geometry";
        if(word==test){
            cs::SelectMaterialByGeometry=true; // default
            // also check for value
            std::string VFalse="false";
            if(value==VFalse){
                cs::SelectMaterialByGeometry=false;
            }
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="fill-core-shell-particles";
        if(word==test){
            cs::fill_core_shell=true; // default
            // also check for value
            std::string VFalse="false";
            if(value==VFalse){
                cs::fill_core_shell=false;
            }
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="interfacial-roughness";
        if(word==test){
            cs::interfacial_roughness=true; // default
            // also check for value
            std::string VFalse="false";
            if(value==VFalse){
                cs::interfacial_roughness=false;
            }
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="material-interfacial-roughness";
        if(word==test){
            cs::interfacial_roughness_local_height_field=true; // default
            // also check for value
            std::string VFalse="false";
            if(value==VFalse){
                cs::interfacial_roughness_local_height_field=false;
            }
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="interfacial-roughness-random-seed";
        if(word==test){
            unsigned int vs=atoi(value.c_str());
            cs::interfacial_roughness_random_seed=vs;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="interfacial-roughness-number-of-seed-points";
        if(word==test){
            int sc=atoi(value.c_str());
            check_for_valid_int(sc, word, line, prefix, 0, 100000,"input","0 - 100,000");
            cs::interfacial_roughness_seed_count=sc;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="interfacial-roughness-type";
        if(word==test){
            std::string loctest="peaks";
            if(value==loctest){
                cs::interfacial_roughness_type=1;
                return EXIT_SUCCESS;
            }
            else
            loctest="troughs";
            if(value==loctest){
                cs::interfacial_roughness_type=-1;
                return EXIT_SUCCESS;
            }
            else{
                cs::interfacial_roughness_type=0;
                return EXIT_SUCCESS;
            }
        }
        //--------------------------------------------------------------------
        test="interfacial-roughness-seed-radius";
        if(word==test){
            double irsr=atof(value.c_str());
            // Test for valid range
            check_for_valid_value(irsr, word, line, prefix, unit, "length", 0.0, 10000.0,"input","0.0 - 1 micrometre");
            cs::interfacial_roughness_mean_seed_radius=irsr;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="interfacial-roughness-seed-radius-variance";
        if(word==test){
            double irsrv=atof(value.c_str());
            // Test for valid range
            check_for_valid_value(irsrv, word, line, prefix, unit, "none", 0.0, 1.0,"input","0.0 - 1.0");
            cs::interfacial_roughness_seed_radius_variance=irsrv;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="interfacial-roughness-mean-height";
        if(word==test){
            double irmh=atof(value.c_str());
            // Test for valid range
            check_for_valid_value(irmh, word, line, prefix, unit, "length", 0.1, 100.0,"input","0.1 Angstroms - 10 nanometres");
            cs::interfacial_roughness_mean_seed_height=irmh;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="interfacial-roughness-maximum-height";
        if(word==test){
            double shm=atof(value.c_str());
            // Test for valid range
            check_for_valid_value(shm, word, line, prefix, unit, "length", 0.1, 100.0,"input","0.1 Angstroms - 10 nanometres");
            cs::interfacial_roughness_seed_height_max=shm;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="interfacial-roughness-height-field-resolution";
        if(word==test){
            double irhfr=atof(value.c_str());
            // Test for valid range
            check_for_valid_value(irhfr, word, line, prefix, unit, "length", 0.1, 100.0,"input","0.1 Angstroms - 10 nanometres");
            cs::interfacial_roughness_height_field_resolution=irhfr;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="multilayers";
        if(word==test){
            int nmul=atoi(value.c_str());
            // Test for valid range
            check_for_valid_int(nmul, word, line, prefix, 1, 100,"input","1 - 100, specifying the number of multilayers to be generated");
            cs::multilayers = true;
            cs::num_multilayers = nmul;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="height-categorization";
        if(word==test){
            // Test for different options
            test="default";
            if(value==test){
                // do nothing
                return EXIT_SUCCESS;
            }
            test="multilayers";
            if(value==test){
                cs::multilayer_height_category = true;
                return EXIT_SUCCESS;
            }
            else{
                terminaltextcolor(RED);
                std::cerr << "Error - value for \'create:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"default\"" << std::endl;
                std::cerr << "\t\"multilayers\"" << std::endl;
                zlog << zTs() << "Error - value for \'create:" << word << "\' must be one of:" << std::endl;
                zlog << zTs() << "\t\"default\"" << std::endl;
                zlog << zTs() << "\t\"multilayers\"" << std::endl;
                terminaltextcolor(WHITE);
                err::vexit();
            }
        }
        //--------------------------------------------------------------------
        // keyword not found
        //--------------------------------------------------------------------
        else{
        terminaltextcolor(RED);
            std::cerr << "Error - Unknown control statement \'create:" << word << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
        return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    }

    int match_dimension(string const word, string const value, string const unit, int const line){
        //-------------------------------------------------------------------
        // System dimension variables
        //-------------------------------------------------------------------
        std::string prefix="dimensions:";
        //--------------------------------------------------------------------
        std::string test="system-size";
        if(word==test){
            double d=atof(value.c_str());
            check_for_valid_value(d, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
            cs::system_dimensions[0]=d;
            cs::system_dimensions[1]=d;
            cs::system_dimensions[2]=d;
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="system-size-x";
        if(word==test){
            double dx=atof(value.c_str());
            check_for_valid_value(dx, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
            cs::system_dimensions[0]=dx;
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="system-size-y";
        if(word==test){
            double dy=atof(value.c_str());
            check_for_valid_value(dy, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
            cs::system_dimensions[1]=dy;
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="system-size-z";
        if(word==test){
            double dz=atof(value.c_str());
            check_for_valid_value(dz, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
            cs::system_dimensions[2]=dz;
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="particle-size";
        if(word==test){
            double psize=atof(value.c_str());
            check_for_valid_value(psize, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
            cs::particle_scale=psize;
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="particle-spacing";
        if(word==test){
            double pspacing=atof(value.c_str());
            check_for_valid_value(pspacing, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
            cs::particle_spacing=pspacing;
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="particle-shape-factor-x";
        if(word==test){
            double sfx=atof(value.c_str());
            check_for_valid_value(sfx, word, line, prefix, unit, "none", 0.001, 2.0,"input","0.001 - 2.0");
            cs::particle_shape_factor_x=sfx;
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="particle-shape-factor-y";
        if(word==test){
            double sfy=atof(value.c_str());
            check_for_valid_value(sfy, word, line, prefix, unit, "none", 0.001, 2.0,"input","0.001 - 2.0");
            cs::particle_shape_factor_y=sfy;
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="particle-shape-factor-z";
        if(word==test){
            double sfz=atof(value.c_str());
            check_for_valid_value(sfz, word, line, prefix, unit, "none", 0.001, 2.0,"input","0.001 - 2.0");
            cs::particle_shape_factor_z=sfz;
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="particle-array-offset-x";
        if(word==test){
            double paox=atof(value.c_str());
            check_for_valid_value(paox, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.0 - 1.0 millimetre");
            cs::particle_array_offset_x=paox;
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="particle-array-offset-y";
        if(word==test){
            double paoy=atof(value.c_str());
            check_for_valid_value(paoy, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.0 - 1.0 millimetre");
            cs::particle_array_offset_y=paoy;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        else{
        terminaltextcolor(RED);
            std::cerr << "Error - Unknown control statement \'dimensions:"<< word << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
        return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    }

    int match_sim(string const word, string const value, string const unit, int const line){

        std::string prefix="sim:";

        //-------------------------------------------------------------------
        // System simulation variables
        //-------------------------------------------------------------------
        std::string test="integrator";
        if(word==test){
            test="llg-heun";
            if(value==test){
                sim::integrator=0;
                return EXIT_SUCCESS;
            }
            test="monte-carlo";
            if(value==test){
                sim::integrator=1;
                return EXIT_SUCCESS;
            }
            test="llg-midpoint";
            if(value==test){
                sim::integrator=2;
                return EXIT_SUCCESS;
            }
            test="constrained-monte-carlo";
            if(value==test){
                sim::integrator=3;
                return EXIT_SUCCESS;
            }
            test="hybrid-constrained-monte-carlo";
            if(value==test){
                sim::integrator=4;
                return EXIT_SUCCESS;
            }
            else{
            terminaltextcolor(RED);
                std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"llg-heun\"" << std::endl;
                std::cerr << "\t\"llg-midpoint\"" << std::endl;
                std::cerr << "\t\"monte-carlo\"" << std::endl;
                std::cerr << "\t\"constrained-monte-carlo\"" << std::endl;
            terminaltextcolor(WHITE);
                err::vexit();
            }
        }
        //-------------------------------------------------------------------
        test="program";
        if(word==test){
            test="benchmark";
            if(value==test){
                sim::program=0;
                return EXIT_SUCCESS;
            }
            test="time-series";
            if(value==test){
                sim::program=1;
                return EXIT_SUCCESS;
            }
            test="hysteresis-loop";
            if(value==test){
                sim::program=2;
                return EXIT_SUCCESS;
            }
            test="static-hysteresis-loop";
            if(value==test){
                sim::program=3;
                return EXIT_SUCCESS;
            }
            test="curie-temperature";
            if(value==test){
                sim::program=4;
                return EXIT_SUCCESS;
            }
            test="field-cool";
            if(value==test){
                sim::program=5;
                return EXIT_SUCCESS;
            }
            test="localised-field-cool";
            if(value==test){
                sim::program=16;
                return EXIT_SUCCESS;
            }
            test="laser-pulse";
            if(value==test){
                sim::program=6;
                return EXIT_SUCCESS;
            }
            test="hamr-simulation";
            if(value==test){
                sim::program=7;
                return EXIT_SUCCESS;
            }
            test="cmc-anisotropy";
            if(value==test){
                sim::program=8;
                return EXIT_SUCCESS;
            }
            test="hybrid-cmc";
            if(value==test){
                sim::program=9;
                return EXIT_SUCCESS;
            }
            test="reverse-hybrid-cmc";
            if(value==test){
                sim::program=10;
                return EXIT_SUCCESS;
            }
            test="LaGrange-Multiplier";
            if(value==test){
                sim::program=11;
                return EXIT_SUCCESS;
            }
            test="partial-hysteresis-loop";
            if(value==test){
                sim::program=12;
                return EXIT_SUCCESS;
            }
            test="localised-temperature-pulse";
            if(value==test){
                sim::program=13;
                return EXIT_SUCCESS;
            }
            test="effective-damping";
            if(value==test){
                sim::program=14;
                return EXIT_SUCCESS;
            }
            test="fmr";
            if(value==test){
                sim::program=15;
                return EXIT_SUCCESS;
            }
            test="diagnostic-boltzmann";
            if(value==test){
                sim::program=50;
                return EXIT_SUCCESS;
            }
            test="setting";
            if(value==test){
                sim::program=51;
                return EXIT_SUCCESS;
            }
            test="disk-tracks";
            if(value==test){
                sim::program=52;
                return EXIT_SUCCESS;
            }
            else{
            terminaltextcolor(RED);
                std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"benchmark\"" << std::endl;
                std::cerr << "\t\"time-series\"" << std::endl;
                std::cerr << "\t\"hysteresis-loop\"" << std::endl;
                std::cerr << "\t\"static-hysteresis-loop\"" << std::endl;
                std::cerr << "\t\"curie-temperature\"" << std::endl;
                std::cerr << "\t\"field-cool\"" << std::endl;
                std::cerr << "\t\"localised-field-cool\"" << std::endl;
                std::cerr << "\t\"laser-pulse\"" << std::endl;
                std::cerr << "\t\"cmc-anisotropy\"" << std::endl;
                std::cerr << "\t\"hybrid-cmc\"" << std::endl;
                std::cerr << "\t\"reverse-hybrid-cmc\"" << std::endl;
                std::cerr << "\t\"localised-temperature-pulse\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
            }
        }
        //-------------------------------------------------------------------
        test="enable-dipole-fields";
        if(word==test){
            sim::hamiltonian_simulation_flags[4]=1;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="enable-fmr-field";
        if(word==test){
            sim::hamiltonian_simulation_flags[5]=1;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="time-step";
        if(word==test){
            double dt=atof(value.c_str());
            check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
            mp::dt_SI=dt;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="total-time-steps";
        if(word==test){
            int tt=atoi(value.c_str());
            check_for_valid_int(tt, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
            sim::total_time=tt;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="loop-time-steps";
        if(word==test){
            int tt=atoi(value.c_str());
            check_for_valid_int(tt, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
            sim::loop_time=tt;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="partial-time-steps";
        if(word==test){
            int tt=atoi(value.c_str());
            check_for_valid_int(tt, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
            terminaltextcolor(YELLOW);
            std::cout << "Warning: Keyword \'partial-time-steps\' is deprecated and may be removed in a future release. Please use \'time-steps-increment\' instead." << std::endl;
            terminaltextcolor(WHITE);
        sim::partial_time=tt;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="time-steps-increment";
        if(word==test){
            int tt=atoi(value.c_str());
            check_for_valid_int(tt, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
            sim::partial_time=tt;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="equilibration-time-steps";
        if(word==test){
            int tt=atoi(value.c_str());
            check_for_valid_int(tt, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
            sim::equilibration_time=tt;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="simulation-cycles";
        if(word==test){
            int r=atoi(value.c_str());
            check_for_valid_int(r, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
            sim::runs=r;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="maximum-temperature";
        if(word==test){
            double T=atof(value.c_str());
            check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
            sim::Tmax=T;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="minimum-temperature";
        if(word==test){
            double T=atof(value.c_str());
            check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
            sim::Tmin=T;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="equilibration-temperature";
        if(word==test){
            double T=atof(value.c_str());
            check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
            sim::Teq=T;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="temperature";
        if(word==test){
            double T=atof(value.c_str());
            check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
            sim::temperature=T;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="temperature-increment";
        if(word==test){
            double T=atof(value.c_str());
            check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
            sim::delta_temperature=T;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="cooling-time";
        if(word==test){
            double T=atof(value.c_str());
            check_for_valid_value(T, word, line, prefix, unit, "time", 1.0e-18, 1.0,"input","1 attosecond - 1 s");
            sim::cooling_time=T;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="laser-pulse-temporal-profile";
        if(word==test){
            test="square";
            if(value==test){
                sim::pump_function=square;
                return EXIT_SUCCESS;
            }
            test="two-temperature";
            if(value==test){
                sim::pump_function=two_temperature;
                return EXIT_SUCCESS;
            }
            test="double-pulse-two-temperature";
            if(value==test){
                sim::pump_function=double_pump_two_temperature;
                return EXIT_SUCCESS;
            }
            test="double-pulse-square";
            if(value==test){
                sim::pump_function=double_pump_square;
                return EXIT_SUCCESS;
            }
            else{
            terminaltextcolor(RED);
                std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"square\"" << std::endl;
                std::cerr << "\t\"double-pulse-square\"" << std::endl;
                std::cerr << "\t\"two-temperature\"" << std::endl;
                std::cerr << "\t\"double-pulse-two-temperature\"" << std::endl;
            terminaltextcolor(WHITE);
                err::vexit();
            }
        }
        //--------------------------------------------------------------------
        test="laser-pulse-time";
        if(word==test){
            double pt=atof(value.c_str());
            check_for_valid_value(pt, word, line, prefix, unit, "time", 1.0e-18, 1.0,"input","1 attosecond - 1 s");
            sim::pump_time=pt;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="laser-pulse-power";
        if(word==test){
            double pp=atof(value.c_str());
            check_for_valid_value(pp, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
            sim::pump_power=pp;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="second-laser-pulse-time";
        if(word==test){
            double pt=atof(value.c_str());
            check_for_valid_value(pt, word, line, prefix, unit, "time", 1.0e-18, 1.0,"input","1 attosecond - 1 s");
            sim::double_pump_time=pt;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="second-laser-pulse-power";
        if(word==test){
            double pp=atof(value.c_str());
            check_for_valid_value(pp, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
            sim::double_pump_power=pp;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="second-laser-pulse-maximum-temperature";
        if(word==test){
            double T=atof(value.c_str());
            check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
            sim::double_pump_Tmax=T;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="second-laser-pulse-delay-time";
        if(word==test){
            double pd=atof(value.c_str());
            check_for_valid_value(pd, word, line, prefix, unit, "time", 0.0, 1.0,"input","0 - 1 s");
            sim::double_pump_delay=pd;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="two-temperature-heat-sink-coupling";
        if(word==test){
            double hscc=atof(value.c_str());
            check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
            sim::HeatSinkCouplingConstant=hscc;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="two-temperature-electron-heat-capacity";
        if(word==test){
            double hscc=atof(value.c_str());
            check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
            sim::TTCe=hscc;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="two-temperature-phonon-heat-capacity";
        if(word==test){
            double hscc=atof(value.c_str());
            check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
            sim::TTCl=hscc;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="two-temperature-electron-phonon-coupling";
        if(word==test){
            double hscc=atof(value.c_str());
            check_for_valid_value(hscc, word, line, prefix, unit, "none", 0.0, 1.0e40,"input","0.0 - 1.0E40");
            sim::TTG=hscc;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="cooling-function";
        if(word==test){
            test="exponential";
            if(value==test){
                sim::cooling_function_flag=0;
                return EXIT_SUCCESS;
            }
            test="gaussian";
            if(value==test){
                sim::cooling_function_flag=1;
                return EXIT_SUCCESS;
            }
            test="double-gaussian";
            if(value==test){
                sim::cooling_function_flag=2;
                return EXIT_SUCCESS;
            }
            test="linear";
            if(value==test){
                sim::cooling_function_flag=3;
                return EXIT_SUCCESS;
            }
            else{
            terminaltextcolor(RED);
                std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"exponential\"" << std::endl;
                std::cerr << "\t\"gaussian\"" << std::endl;
                std::cerr << "\t\"double-gaussian\"" << std::endl;
                std::cerr << "\t\"linear\"" << std::endl;
            terminaltextcolor(WHITE);
                err::vexit();
            }
        }
        //--------------------------------------------------------------------
        test="applied-field-strength";
        if(word==test){
            double H=atof(value.c_str());
            check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
            sim::H_applied=H;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="minimum-applied-field-strength";
        if(word==test){
            double H=atof(value.c_str());
            check_for_valid_value(H, word, line, prefix, unit, "field", 0.0, 1.0e3,"input","0 - 1,000 T");
            sim::Hmin=H;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="maximum-applied-field-strength";
        if(word==test){
            double H=atof(value.c_str());
            check_for_valid_value(H, word, line, prefix, unit, "field", 0.0, 1.0e3,"input","0 - 1,000 T");
            sim::Hmax=H;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="equilibration-applied-field-strength";
        if(word==test){
            double H=atof(value.c_str());
            check_for_valid_value(H, word, line, prefix, unit, "field", 0.0, 1.0e3,"input","0 - 1,000 T");
            sim::Heq=H;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="applied-field-strength-increment";
        if(word==test){
            double H=atof(value.c_str());
            check_for_valid_value(H, word, line, prefix, unit, "field", 1.0e-6, 1.0e3,"input","1 uT - 1,000 T");
            sim::Hinc=H;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="applied-field-angle-theta";
        if(word==test){
            double angle=atof(value.c_str());
            check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
            sim::applied_field_angle_theta=angle;
            sim::applied_field_set_by_angle=true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="applied-field-angle-phi";
        if(word==test){
            double angle=atof(value.c_str());
            check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
            sim::applied_field_angle_phi=angle;
            sim::applied_field_set_by_angle=true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="applied-field-unit-vector";
        if(word==test){
            std::vector<double> u(3);
            u=doubles_from_string(value);
            check_for_valid_unit_vector(u, word, line, prefix, "input");
            sim::H_vec[0]=u.at(0);
            sim::H_vec[1]=u.at(1);
            sim::H_vec[2]=u.at(2);
            sim::applied_field_set_by_angle=false;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="demagnetisation-factor";
        if(word==test){
            std::vector<double> u(3);
            u=doubles_from_string(value);
            vin::check_for_valid_three_vector(u, word, line, prefix, "input");
            // Extra check for demagnetisation-factor Nx+Ny+Nz=1
            double sum=u.at(0)+u.at(1)+u.at(2);
            if(fabs(1.0-sum)>1.e-4){
            terminaltextcolor(RED);
                std::cerr << "Error: sum of all elements of variable " << prefix << word << " on line " << line << " of input file must equal 1." << std::endl;
                terminaltextcolor(WHITE);
            zlog << zTs() << "Error: sum of all elements of variable " << prefix << word << " on line " << line << " of input file must equal 1." << std::endl;
                err::vexit();
            }
            sim::demag_factor[0]=u.at(0);
            sim::demag_factor[1]=u.at(1);
            sim::demag_factor[2]=u.at(2);
            sim::ext_demag=true;
            // force calculation of system magnetization
            stats::calculate_system_magnetization=true;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="mpi-mode";
        if(word==test){
            test="geometric-decomposition";
            if(value==test){
                vmpi::mpi_mode=0;
                return EXIT_SUCCESS;
            }
            test="replicated-data";
            if(value==test){
                vmpi::mpi_mode=1;
                vmpi::replicated_data_staged=false;
                return EXIT_SUCCESS;
            }
            test="replicated-data-staged";
            if(value==test){
                vmpi::mpi_mode=1;
                vmpi::replicated_data_staged=true;
                return EXIT_SUCCESS;
            }
            else{
            terminaltextcolor(RED);
                std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"geometric-decomposition\"" << std::endl;
                std::cerr << "\t\"replicated-data\"" << std::endl;
                std::cerr << "\t\"replicated-data-staged\"" << std::endl;
            terminaltextcolor(WHITE);
                err::vexit();
            }
        }
        //--------------------------------------------------------------------
        test="mpi-ppn";
        if(word==test){
            int ppn=atoi(value.c_str());
            check_for_valid_int(ppn, word, line, prefix, 1, 1024,"input","1 - 1024");
            vmpi::ppn=ppn;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="integrator-random-seed";
        if(word==test){
            int is=atoi(value.c_str());
            check_for_valid_int(is, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
            mtrandom::integration_seed=is;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="constraint-rotation-update";
        if(word==test){
            sim::constraint_rotation=true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="constraint-angle-theta";
        if(word==test){
            double angle=atof(value.c_str());
            check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
            sim::constraint_theta=angle;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="constraint-angle-theta-minimum";
        if(word==test){
            double angle=atof(value.c_str());
            check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
            sim::constraint_theta_min=angle;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="constraint-angle-theta-maximum";
        if(word==test){
            double angle=atof(value.c_str());
            check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
            sim::constraint_theta_max=angle;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="constraint-angle-theta-increment";
        if(word==test){
            double angle=atof(value.c_str());
            check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
            sim::constraint_theta_delta=angle;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="constraint-angle-phi";
        if(word==test){
            double angle=atof(value.c_str());
            check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
            sim::constraint_phi=angle;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="constraint-angle-phi-minimum";
        if(word==test){
            double angle=atof(value.c_str());
            check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
            sim::constraint_phi_min=angle;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="constraint-angle-phi-maximum";
        if(word==test){
            double angle=atof(value.c_str());
            check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
            sim::constraint_phi_max=angle;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="constraint-angle-phi-increment";
        if(word==test){
            double angle=atof(value.c_str());
            check_for_valid_value(angle, word, line, prefix, unit, "none", 0.0, 360.0,"input","0.0 - 360.0 degrees");
            sim::constraint_phi_delta=angle;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="monte-carlo-algorithm";
        if(word==test){
            // include namesapce here to access enum values
            using namespace sim;
            test="spin-flip";
            if(value==test){
                sim::mc_algorithm=spin_flip;
                return EXIT_SUCCESS;
            }
            test="uniform";
            if(value==test){
                sim::mc_algorithm=uniform;
                return EXIT_SUCCESS;
            }
            test="angle";
            if(value==test){
                sim::mc_algorithm=angle;
                return EXIT_SUCCESS;
            }
            test="hinzke-nowak";
            if(value==test){
                sim::mc_algorithm=hinzke_nowak;
                return EXIT_SUCCESS;
            }
            else{
            terminaltextcolor(RED);
                std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"spin-flip\"" << std::endl;
                std::cerr << "\t\"uniform\"" << std::endl;
                std::cerr << "\t\"angle\"" << std::endl;
                std::cerr << "\t\"hinzke-nowak\"" << std::endl;
            terminaltextcolor(WHITE);
                err::vexit();
            }
        }
        //-------------------------------------------------------------------
        test="save-checkpoint";
        if(word==test){
            test="end";
            if(value==test){
                sim::save_checkpoint_flag=true; // Save checkpoint
                sim::save_checkpoint_continuous_flag=false; // do not save checkpoints during simulation
                return EXIT_SUCCESS;
            }
            test="continuous";
            if(value==test){
                sim::save_checkpoint_flag=true; // Save checkpoint
                sim::save_checkpoint_continuous_flag=true; // save checkpoints during simulation
                return EXIT_SUCCESS;
            }
            else{
                terminaltextcolor(RED);
                std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"end\"" << std::endl;
                std::cerr << "\t\"continuous\"" << std::endl;
                terminaltextcolor(WHITE);
                err::vexit();
            }
        }
        //--------------------------------------------------------------------
        test="save-checkpoint-rate";
        if(word==test){
            int scr=atoi(value.c_str());
            check_for_valid_int(scr, word, line, prefix, 1, 2000000000,"input","1 - 2,000,000,000");
            sim::save_checkpoint_rate=scr;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="load-checkpoint";
        if(word==test){
            test="restart";
            if(value==test){
                sim::load_checkpoint_flag=true; // Load spin configurations
                sim::load_checkpoint_continue_flag=false; // Restart simulation with checkpoint configuration
                return EXIT_SUCCESS;
            }
            test="continue";
            if(value==test){
                sim::load_checkpoint_flag=true; // Load spin configurations
                sim::load_checkpoint_continue_flag=true; // Continue simulation from saved time with checkpoint configuration
                return EXIT_SUCCESS;
            }
            else{
                terminaltextcolor(RED);
                std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"restart\"" << std::endl;
                std::cerr << "\t\"continue\"" << std::endl;
                terminaltextcolor(WHITE);
                err::vexit();
            }
        }
        //--------------------------------------------------------------------
        test="fmr-field-strength";
        if(word==test){
            double H=atof(value.c_str());
            check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
            sim::fmr_field_strength=H;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="fmr-field-frequency";
        if(word==test){
            double w = atof(value.c_str());
            check_for_valid_value(w, word, line, prefix, unit, "none", 0.0, 1.0e4,"input","0 - 10,000 GHz");
            sim::fmr_field_frequency = w;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="fmr-field-unit-vector";
        if(word==test){
            std::vector<double> u(3);
            u=doubles_from_string(value);
            check_for_valid_unit_vector(u, word, line, prefix, "input");
            sim::fmr_field_unit_vector = u;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        else{
        terminaltextcolor(RED);
            std::cerr << "Error - Unknown control statement \'sim:"<< word << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
        return EXIT_FAILURE;
        }


        return EXIT_SUCCESS;
    }

    int match_vout_list(string const word, string const value, int const line, std::vector<unsigned int> & output_list){

        std::string prefix="output:";

        std::string test="time-steps";
        if(word==test){
            output_list.push_back(0);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="real-time";
        if(word==test){
            output_list.push_back(1);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="temperature";
        if(word==test){
            output_list.push_back(2);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="applied-field-strength";
        if(word==test){
            output_list.push_back(3);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="applied-field-unit-vector";
        if(word==test){
            output_list.push_back(4);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="applied-field-alignment";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(12);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="magnetisation";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(5);
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="magnetisation-length";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(6);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="mean-magnetisation-length";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(7);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="mean-magnetisation";
        if(word==test){
            stats::calculate_system_magnetization=true;
            output_list.push_back(48);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="material-magnetisation";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(8);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="material-mean-magnetisation-length";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(9);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="material-mean-magnetisation";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(49);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="total-torque";
        if(word==test){
            stats::calculate_torque=true;
            output_list.push_back(14);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="mean-total-torque";
        if(word==test){
            stats::calculate_torque=true;
            output_list.push_back(15);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="constraint-phi";
        if(word==test){
            stats::calculate_torque=true;
            output_list.push_back(16);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="constraint-theta";
        if(word==test){
            stats::calculate_torque=true;
            output_list.push_back(17);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="material-constraint-phi";
        if(word==test){
            stats::calculate_torque=true;
            output_list.push_back(18);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="material-constraint-theta";
        if(word==test){
            stats::calculate_torque=true;
            output_list.push_back(19);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="material-mean-torque";
        if(word==test){
            stats::calculate_torque=true;
            output_list.push_back(20);
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="mean-susceptibility";
        if(word==test){
            // Set flags for calculations of susceptibility and magnetization
            stats::calculate_system_susceptibility=true;
            stats::calculate_system_magnetization=true;
            output_list.push_back(21);
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="material-mean-susceptibility";
        if(word==test){
            // Set flags for calculations of susceptibility and magnetization
            stats::calculate_material_susceptibility=true;
            stats::calculate_material_magnetization=true;
            output_list.push_back(50);
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="electron-temperature"; // identical to temperature
        if(word==test){
            output_list.push_back(2);
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="phonon-temperature";
        if(word==test){
            output_list.push_back(22);
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="material-temperature";
        if(word==test){
            output_list.push_back(23);
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="material-applied-field-strength";
        if(word==test){
            output_list.push_back(24);
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="material-fmr-field-strength";
        if(word==test){
            output_list.push_back(25);
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="material-applied-field-alignment";
        if(word==test){
            stats::calculate_material_magnetization=true;
            output_list.push_back(26);
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="total-energy";
        if(word==test){
            output_list.push_back(27);
            stats::calculate_energy=true;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="mean-total-energy";
        if(word==test){
            output_list.push_back(28);
            stats::calculate_energy=true;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="anisotropy-energy";
        if(word==test){
            output_list.push_back(29);
            stats::calculate_energy=true;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="mean-anisotropy-energy";
        if(word==test){
            output_list.push_back(30);
            stats::calculate_energy=true;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="exchange-energy";
        if(word==test){
            output_list.push_back(35);
            stats::calculate_energy=true;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="mean-exchange-energy";
        if(word==test){
            output_list.push_back(36);
            stats::calculate_energy=true;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="applied-field-energy";
        if(word==test){
            output_list.push_back(37);
            stats::calculate_energy=true;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="mean-applied-field-energy";
        if(word==test){
            output_list.push_back(38);
            stats::calculate_energy=true;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="magnetostatic-energy";
        if(word==test){
            output_list.push_back(39);
            stats::calculate_energy=true;
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="mean-magnetostatic-energy";
        if(word==test){
            output_list.push_back(40);
            stats::calculate_energy=true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="height-magnetisation-normalised";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(43);
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="material-height-magnetisation-normalised";
        if(word==test){
            stats::calculate_material_height_magnetization=true;
            output_list.push_back(44);
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="height-magnetisation";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(45);
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="material-height-magnetisation";
        if(word==test){
            stats::calculate_material_height_magnetization=true;
            output_list.push_back(46);
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="fmr-field-strength";
        if(word==test){
            output_list.push_back(47);
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="mean-height-magnetisation-length";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(51);
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="mean-height-magnetisation";
        if(word==test){
            stats::calculate_height_magnetization=true;
            output_list.push_back(52);
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="mpi-timings";
        if(word==test){
            vmpi::DetailedMPITiming=true;
            output_list.push_back(60);
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="gnuplot-array-format";
        if(word==test){
            vout::gnuplot_array_format=true;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        test="output-rate";
        if(word==test){
            int r=atoi(value.c_str());
            check_for_valid_int(r, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
            vout::output_rate=r;
            return EXIT_SUCCESS;
        }

        //--------------------------------------------------------------------
        // keyword not found
        //--------------------------------------------------------------------
        else{
            terminaltextcolor(RED);
            std::cerr << "Error - Unknown control statement "<< prefix << word << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
        return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    }
    int match_vout_grain_list(string const word, string const value, int const line, std::vector<unsigned int> & output_list){

        std::string prefix="grain:";

        std::string test="time-steps";
        if(word==test){
            output_list.push_back(0);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="real-time";
        if(word==test){
            output_list.push_back(1);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="temperature";
        if(word==test){
            output_list.push_back(2);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="applied-field-strength";
        if(word==test){
            output_list.push_back(3);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="applied-field-unit-vector";
        if(word==test){
            output_list.push_back(4);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="magnetisation";
        if(word==test){
            output_list.push_back(10);
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="mag-m";
        if(word==test){
            output_list.push_back(11);
            return EXIT_SUCCESS;
        }
        else
        //--------------------------------------------------------------------
        test="material-magnetisation";
        if(word==test){
            output_list.push_back(13);
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="electron-temperature"; // identical to temperature
        if(word==test){
            output_list.push_back(2);
            return EXIT_SUCCESS;
        }
        //-------------------------------------------------------------------
        test="phonon-temperature";
        if(word==test){
            output_list.push_back(22);
            return EXIT_SUCCESS;
        }
        else
        //-------------------------------------------------------------------
        test="output-rate";
        if(word==test){
            int r=atoi(value.c_str());
            check_for_valid_int(r, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
            vout::output_grain_rate=r;
            return EXIT_SUCCESS;
        }
        //--------------------------------------------------------------------
        // keyword not found
        //--------------------------------------------------------------------
        else{
        terminaltextcolor(RED);
            std::cerr << "Error - Unknown control statement \'grain:" << word << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
        return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    }

    // temporary array of materials for reading in material data
    std::vector<mp::materials_t> read_material(0);

    int read_mat_file(std::string const matfile, int const line_number){


        // resize temporary materials array for storage of variables
        read_material.resize(mp::max_materials);
        cmc::cmc_mat.resize(mp::max_materials);

            // Print informative message to zlog file
        zlog << zTs() << "Opening material file \"" << matfile << "\"." << std::endl;

        // Open file read only
        std::stringstream inputfile;
        inputfile.str( vin::get_string(matfile.c_str(), "material", line_number) );
        //-------------------------------------------------------
        // Material 0
        //-------------------------------------------------------


        // Print informative message to zlog file
        zlog << zTs() << "Parsing material file for parameters." << std::endl;

        int line_counter=0;
        // Loop over all lines and pass keyword to matching function
        while (! inputfile.eof() ){
            line_counter++;
            // read in whole line
            std::string line;
            getline(inputfile,line);

            // save a copy of the line before stripping characters in case of error
            std::string original_line = line;

            // Clear whitespace, quotes and tabs
            line.erase(remove(line.begin(), line.end(), '\t'), line.end());
            line.erase(remove(line.begin(), line.end(), ' '), line.end());
            line.erase(remove(line.begin(), line.end(), '\"'), line.end());

            // remove carriage returns for dos formatted files
                        line.erase(remove(line.begin(), line.end(), '\r'), line.end());

            // strip key,word,unit,value
            std::string key="";
            std::string word="";
            std::string value="";
            std::string unit="";
            std::string index="";
            int super_index=1; // Inital values *as would be read from input file*
            int sub_index=1;

            // get size of string
            int linelength = line.length();
            int last=0;

            // set character triggers
            const char* colon=":";	// Word identifier
            const char* eq="=";		// Value identifier
            const char* exc="!";		// Unit identifier
            const char* hash="#";	// Comment identifier
            const char* si="[";		// Index identifier
            const char* ei="]";		// End index identifier

            // Determine key and super index by looping over characters in line
            for(int i=0;i<linelength;i++){
                char c=line.at(i);
                last=i;

                // if character is not ":" or "=" or "!" or "#" interpret as key
                if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash) && (c != *si) && (c != *ei)){
                    key.push_back(c);
                }
                // Check for number of materials statement
                else if(c == *eq){
                    // break to read in value
                    break;
                }
                // Check for superindex
                else if(c ==*si){
                    const int old=last;
                    // Get super index
                    for(int j=old+1;j<linelength;j++){
                        c=line.at(j);
                        if(c != *ei){
                            index.push_back(c);
                        }
                        else{
                            break;
                        }
                        last=j;
                    }

                    // check for valid index
                    super_index = atoi(index.c_str());
                    if((super_index>=1) && (super_index<mp::max_materials+1)){
                        break;
                    }
                    else{
                        std::cerr << "Invalid index number " << index << " on line " << line_counter << " in material input file" << std::endl;
                        std::cerr << "Causes could be invalid character or outside of range, ie less than 1 or greater than max_materials=" << mp::max_materials << ", exiting" << std::endl;
                        err::vexit();
                    }

                }
                // For anything else
                else break;
            }
            const int end_key=last;

            //
            //err::vexit();
            // Determine the rest
            for(int i=end_key;i<linelength;i++){

                char c=line.at(i);
                // colon found - interpret as word
                if(c== *colon){
                    for(int j=i+1;j<linelength;j++){
                        // if character is not special add to value
                        char c=line.at(j);
                        if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
                            // check for sub-index
                            if(c == *si){
                                index="";
                                while(line.at(j+1) != *ei){
                                    j++;
                                    index.push_back(line.at(j));
                                }
                                sub_index=atoi(index.c_str());
                                // Check for valid index
                                if((sub_index<1) || (sub_index>=mp::max_materials)){
                                    std::cerr << "Invalid sub-index number " << index << " on line " << line_counter << " in material input file" << std::endl;
                                    std::cerr << "Causes could be invalid character or outside of range, ie less than 1 or greater than max_materials=" << mp::max_materials << ", exiting" << std::endl;
                                    err::vexit();
                                }
                                // end of word
                                break;
                            }
                            else word.push_back(c);
                        }
                        // if character is special then go back to main loop
                        else{
                            i=j-1;
                            break;
                        }
                    }
                }
                // equals found - interpret as value
                else if(c== *eq){
                    for(int j=i+1;j<linelength;j++){
                        // if character is not special add to value
                        char c=line.at(j);
                        if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
                            value.push_back(c);
                        }
                        // if character is special then go back to main loop
                        else{
                            i=j-1;
                            break;
                        }
                    }
                }
                // exclaimation mark found - interpret as unit
                else if(c== *exc){
                    for(int j=i+1;j<linelength;j++){
                        // if character is not special add to value
                        char c=line.at(j);
                        if((c != *colon) && (c != *eq) && (c != *exc) && (c != *hash)){
                            unit.push_back(c);
                        }
                        // if character is special then go back to main loop
                        else{
                            i=j-1;
                            break;
                        }
                    }
                }
                // hash found - interpret as comment
                else if(c== *hash){
                    break;
                }
                //break;
            }
            string empty="";
            if(key!=empty){
                //std::cout << key << "[" << super_index << "]:" << word << "[" << sub_index << "]=" << value << " !" << unit << std::endl;
                //std::cout << "\t" << "key:  " << key << std::endl;
                //std::cout << "\t" << "word: " << word << std::endl;
                //std::cout << "\t" << "value:" << value << std::endl;
                //std::cout << "\t" << "unit: " << unit << std::endl;
            int matchcheck = vin::match_material(word, value, unit, line_counter, super_index-1, sub_index-1, original_line, matfile);
                if(matchcheck==EXIT_FAILURE){
                    err::vexit();
                }
            }
        }

        // resize global material array
        mp::material.resize(mp::num_materials);

        // Copy data to global material array
        for(int mat=0;mat<mp::num_materials;mat++){
            mp::material[mat]=read_material[mat];
        }

        // Resize read array to zero
        read_material.resize(0);


        // Close file
        //inputfile.close();

        return EXIT_SUCCESS;

    }
    ///-------------------------------------------------------------------
    /// Function to match material key words
    ///-------------------------------------------------------------------
    int match_material(string const word,
                string const value,
                string const unit,
                int const line,
                int const super_index,
                int const sub_index,
                std::string const line_string,
                std::string const filename_string)
    {
            std::string prefix="material:";
            //------------------------------------------------------------
            std::string test="num-materials";
            if(word==test){
                unsigned int nm = atoi(value.c_str());
                check_for_valid_int(nm, word, line, prefix, 1, 100,"material","1 - 100");
                mp::num_materials=nm;
                return EXIT_SUCCESS;
            }
            //------------------------------------------------------------
            else
            test="material-name";
            if(word==test){
                read_material[super_index].name=value;
                return EXIT_SUCCESS;
            }
            //------------------------------------------------------------
            else
            test="damping-constant";
            if(word==test){
                double damping=atof(value.c_str());
                check_for_valid_positive_value(damping, word, line, prefix, unit, "none", 0.0, 10.0,"material","0.0 - 10.0");
                read_material[super_index].alpha=damping;
                return EXIT_SUCCESS;
            }
            //------------------------------------------------------------
            else
            test="atomic-spin-moment";
            if(word==test){
                double mu_s=atof(value.c_str());
                check_for_valid_positive_value(mu_s, word, line, prefix, unit, "moment", 0.1*9.24e-24, 1e8*9.24e-24,"material","0.1 - 1e8 mu_B");
                read_material[super_index].moment_flag=true;
                read_material[super_index].mu_s_SI=mu_s;
                return EXIT_SUCCESS;
            }
            //------------------------------------------------------------
            else
            test="relative-gamma";
            if(word==test){
                double gr = atof(value.c_str());
                // Test for valid range
                check_for_valid_positive_value(gr, word, line, prefix, unit, "none", 0.01, 100.0,"material"," 0.01 - 100.0");
                read_material[super_index].gamma_rel=gr;
                return EXIT_SUCCESS;
            }
            //------------------------------------------------------------
            else
            test="initial-spin-direction";
            if(word==test){
                // first test for random spins
                test="random";
                if(value==test){
                    read_material[super_index].random_spins=true;
                }
                else{
                    // temporary storage container
                    std::vector<double> u(3);

                    // read values from string
                    u=doubles_from_string(value);

                    // check for sane input and normalise if necessary
                    check_for_valid_unit_vector(u, word, line, prefix, "material");

                    // Copy sanitised unit vector to material
                    read_material[super_index].initial_spin[0]=u.at(0);
                    read_material[super_index].initial_spin[1]=u.at(1);
                    read_material[super_index].initial_spin[2]=u.at(2);

                    // ensure random spins is unset
                    read_material[super_index].random_spins=false;
                }
                // return
                return EXIT_SUCCESS;
            }
            //------------------------------------------------------------
            else
            test="material-element";
            if(word==test){
                // Test for 3 characters
                if(value.length()>3){
                terminaltextcolor(RED);
                    std::cerr << "Error - element identifier on line "<<  line << " of material file must be a maximum of three characters long" << std::endl;
                    terminaltextcolor(WHITE);
            }
                else{
                    // pad value to be equal to 3 characters
                    string tmp="   ";
                    for(unsigned int i=0;i<3;i++){
                        if(i<value.length()){
                            tmp.at(i)=value.at(i);
                        }
                    }
                    read_material[super_index].element=tmp;
                    return EXIT_SUCCESS;
                }
            }
            //--------------------------------------------------------------------
            else
            test="geometry-file";
            if(word==test){

               // Open geometry file
               std::stringstream gfile;
          		gfile.str( vin::get_string(value.c_str(), "material", line) );

                gfile >> read_material[super_index].geometry;
                if((read_material[super_index].geometry<3) || (read_material[super_index].geometry>100)){
                    terminaltextcolor(RED);
                std::cerr << "Error in geometry input file " << value.c_str() << " - first number must be non zero integer in the range 3-100"<< std::endl;
                    terminaltextcolor(WHITE);
                return EXIT_FAILURE;
                }
                //std::cout << "ngp " << read_material[super_index].geometry << std::endl;
                for(int c=0;c<read_material[super_index].geometry;c++){
                    for(int xy=0;xy<2;xy++){
                        double var;
                        gfile >> var;
                        if(gfile.eof()){
                    terminaltextcolor(RED);
                            std::cerr << "Error in geometry input file " << value.c_str() << " end of file reached before reading all coordinates" << std::endl;
                            terminaltextcolor(WHITE);
                    return EXIT_FAILURE;
                        }
                        if((var<0.0) || (var > 1.0)){
                    terminaltextcolor(RED);
                            std::cerr << "Error in geometry input file " << value.c_str() << " value is outside of valid range (0.0-1.0)" << std::endl;
                            terminaltextcolor(WHITE);
                    return EXIT_FAILURE;
                        }
                        else read_material[super_index].geometry_coords[c][xy]=var;
                    }
                    //std::cout << read_material[super_index].geometry_coords[c][0] << "\t" << read_material[super_index].geometry_coords[c][1] << std::endl;
                }
                //double min=atof(value.c_str());
                //if((min<-0.11) || (min > 1.11)){
                //	std::cerr << "Error in input file - material[" << super_index << "]:min is outside of valid range (0.0-1.0)" << std::endl;
                //	return EXIT_FAILURE;}
                //else{
                //	read_material[super_index].min=min;
                    return EXIT_SUCCESS;
                //}
            }
            //--------------------------------------------------------------------
            test="core-shell-size";
            if(word==test){
                double css=atof(value.c_str());
                check_for_valid_positive_value(css, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
                read_material[super_index].core_shell_size=css;
                cs::core_shell_particles = true;
                return EXIT_SUCCESS;
            }
            //-------------------------------------------------------------------
            else
            test="interface-roughness";
            if(word==test){
                double ir=atof(value.c_str());
                check_for_valid_value(ir, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
                read_material[super_index].interface_roughness=ir;
                return EXIT_SUCCESS;
            }
            else
            //-------------------------------------------------------------------
            test="density";
            if(word==test){
                double d=atof(value.c_str());
                check_for_valid_positive_value(d, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
                read_material[super_index].density=d;
                return EXIT_SUCCESS;
            }
            else
            test="continuous";
            if(word==test){
                read_material[super_index].continuous=true;
                return EXIT_SUCCESS;
            }
            else
            //-------------------------------------------------------------------
            test="intermixing";
            if(word==test){
                double i=atof(value.c_str());
                //check_for_valid_value(i, word, line, prefix, unit, "none", 0.0, 1.0,"material"," 0.0 - 1.0");
                if((i<0.0) || (i > 1.0)){
                terminaltextcolor(RED);
                    std::cerr << "Error in input file - material[" << super_index+1 << "]:intermixing[" << sub_index+1 <<"] is outside of valid range (0.0-1.0)" << std::endl;
                terminaltextcolor(WHITE);
                    return EXIT_FAILURE;}
                else{
                    read_material[super_index].intermixing[sub_index]=i;
                    return EXIT_SUCCESS;
                }
            }
            //--------------------------------------------------------------------
            else
            test="constrained"; // determines use of alternate integrator
            if(word==test){
                string t="true";
                string f="false";
                if(value==t){
                    read_material[super_index].constrained=true;
                    return EXIT_SUCCESS;
                }
                else if(value==f){
                    read_material[super_index].constrained=false;
                    return EXIT_SUCCESS;
                }
                else {
                    terminaltextcolor(RED);
                    std::cerr << "Error in input file - material[" << super_index+1 << "]:constrained must be either true or false" << std::endl;
                    terminaltextcolor(WHITE);
                    return EXIT_FAILURE;
                }
            }
            //--------------------------------------------------------------------
            test="constraint-angle-theta";
            if(word==test){
                double angle=atof(value.c_str());
                // Test for valid range
                if((angle>=0.0) && (angle<=360.0)){
                    cmc::cmc_mat[super_index].constraint_theta=angle;
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 360.0" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //--------------------------------------------------------------------
            test="constraint-angle-theta-minimum";
            if(word==test){
                double angle=atof(value.c_str());
                // Test for valid range
                if((angle>=0.0) && (angle<=360.0)){
                    cmc::cmc_mat[super_index].constraint_theta_min=angle;
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 360.0" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //--------------------------------------------------------------------
            test="constraint-angle-theta-maximum";
            if(word==test){
                double angle=atof(value.c_str());
                // Test for valid range
                if((angle>=0.0) && (angle<=360.0)){
                    cmc::cmc_mat[super_index].constraint_theta_max=angle;
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 360.0" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //--------------------------------------------------------------------
            test="constraint-angle-theta-increment";
            if(word==test){
                double angle=atof(value.c_str());
                // Test for valid range
                if((angle>=0.0) && (angle<=360.0)){
                    cmc::cmc_mat[super_index].constraint_theta_delta=angle;
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 360.0" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //--------------------------------------------------------------------
            test="constraint-angle-phi-minimum";
            if(word==test){
                double angle=atof(value.c_str());
                // Test for valid range
                if((angle>=0.0) && (angle<=180.0)){
                    cmc::cmc_mat[super_index].constraint_phi_min=angle;
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 180.0" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //--------------------------------------------------------------------
            test="constraint-angle-phi";
            if(word==test){
                double angle=atof(value.c_str());
                // Test for valid range
                if((angle>=0.0) && (angle<=180.0)){
                    cmc::cmc_mat[super_index].constraint_phi=angle;
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 180.0" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //--------------------------------------------------------------------
            test="constraint-angle-phi-maximum";
            if(word==test){
                double angle=atof(value.c_str());
                // Test for valid range
                if((angle>=0.0) && (angle<=180.0)){
                    cmc::cmc_mat[super_index].constraint_phi_max=angle;
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 180.0" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //--------------------------------------------------------------------
            test="constraint-angle-phi-increment";
            if(word==test){
                double angle=atof(value.c_str());
                // Test for valid range
                if((angle>=0.0) && (angle<=180.0)){
                    cmc::cmc_mat[super_index].constraint_phi_delta=angle;
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 180.0" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //--------------------------------------------------------------------
            test="temperature";
            if(word==test){
                double T=atof(value.c_str());
                // Test for valid range
                if((T>=0.0) && (T<1.0E5)){
                    read_material[super_index].temperature=T;
                    // set local temperature flag
                    sim::local_temperature=true;
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 1.0E5" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //--------------------------------------------------------------------
            test="maximum-temperature";
            if(word==test){
               double T=atof(value.c_str());
               check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
               read_material[super_index].maximum_temperature = T;
               sim::local_temperature=true; // set local temperature flag
               return EXIT_SUCCESS;
            }
            //--------------------------------------------------------------------
            test="minimum-temperature";
            if(word==test){
               double T=atof(value.c_str());
               check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e6,"input","0.0 - 1,000,000 K");
               read_material[super_index].minimum_temperature = T;
               sim::local_temperature=true; // set local temperature flag
               return EXIT_SUCCESS;
            }
            //--------------------------------------------------------------------
            test="use-phonon-temperature";
            /*
            logical use-phonon-temperature
                This flag enables specific materials to couple to the phonon temperature
                of the system for simulations using the two temperature model. The default
                is for all materials to use the electron temperature. Valid values are true,
                false or (blank) [same as true].
            */
            if(word==test){
                // Test for sane input
                //bool sanitised_bool=check_for_valid_bool(value, word, line, prefix,"material");

                // set flag
                read_material[super_index].couple_to_phonon_temperature=true; //sanitised_bool;

                // enable local temperature flag
                sim::local_temperature=true;
                return EXIT_SUCCESS;
            }
            //--------------------------------------------------------------------
            test="fill-space";
            /*
            logical fill-space [false]
                This flag causes the material to be a fill material, instead of
                have atoms generated as part of the usual structure, for example
                particle shpaes, particle arrays and voronoi films. The default
                value is false for all materials. Valid values are true,
                false or (blank) [same as true].
            */
            if(word==test){
                // Test for sane input
                //bool sanitised_bool=check_for_valid_bool(value, word, line, prefix,"material");

                // set flag
                read_material[super_index].fill=true; //sanitised_bool;

                return EXIT_SUCCESS;
            }
            //--------------------------------------------------------------------
            test="applied-field-strength";
            if(word==test){
                double H=atof(value.c_str());
                // test for unit
                string unit_type="field";
                // if no unit given, assume internal
                if(unit.size() != 0){
                    units::convert(unit,H,unit_type);
                }
                string str="field";
                if(unit_type==str){
                    // Test for valid range
                    if((H>=0.0) && (H<1.0E5)){
                        read_material[super_index].applied_field_strength=H;
                        // set local applied field flag
                        sim::local_applied_field=true;
                        return EXIT_SUCCESS;
                    }
                    else{
                        terminaltextcolor(RED);
                        std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 1.0E5" << std::endl;
                        terminaltextcolor(WHITE);
                        err::vexit();
                    }
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - unit type \'" << unit_type << "\' is invalid for parameter material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 1.0E5" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //------------------------------------------------------------
            test="applied-field-unit-vector";
            if(word==test){
                // temporary storage container
                std::vector<double> u(3);

                // read values from string
                u=doubles_from_string(value);

                // check size
                if(u.size()!=3){
                    terminaltextcolor(RED);
                    std::cerr << "Error in input file - material[" << super_index+1 << "]:"<< word << " must have three values." << std::endl;
                    terminaltextcolor(WHITE);
                    zlog << zTs() << "Error in input file - material[" << super_index+1 << "]:"<< word << " must have three values." << std::endl;
                    return EXIT_FAILURE;
                }

                // Normalise
                double ULength=sqrt(u.at(0)*u.at(0)+u.at(1)*u.at(1)+u.at(2)*u.at(2));

                // Check for correct length unit vector
                if(ULength < 1.0e-9){
                    terminaltextcolor(RED);
                    std::cerr << "Error in input file - material[" << super_index+1 << "]:"<< word << " must be normalisable (possibly all zero)." << std::endl;
                    terminaltextcolor(WHITE);
                    zlog << zTs() << "Error in input file - material[" << super_index+1 << "]:"<< word << " must be normalisable (possibly all zero)." << std::endl;
                    return EXIT_FAILURE;
                }
                u.at(0)/=ULength;
                u.at(1)/=ULength;
                u.at(2)/=ULength;

                // Copy applied field direction to material
                read_material[super_index].applied_field_unit_vector=u;

                // set local applied field flag
                sim::local_applied_field=true;

                return EXIT_SUCCESS;

            }
            //--------------------------------------------------------------------
            test="fmr-field-strength";
            if(word==test){
                double H=atof(value.c_str());
                // test for unit
                string unit_type="field";
                // if no unit given, assume internal
                if(unit.size() != 0){
                    units::convert(unit,H,unit_type);
                }
                string str="field";
                if(unit_type==str){
                    // Test for valid range
                    if((H>=0.0) && (H<1.0E5)){
                        read_material[super_index].applied_field_strength=H;
                        // set local fmr flag
                        sim::local_fmr_field=true;
                        return EXIT_SUCCESS;
                    }
                    else{
                        terminaltextcolor(RED);
                        std::cerr << "Error - sim:" << word << " on line " << line << " of input file must be in the range 0 - 1.0E5" << std::endl;
                        terminaltextcolor(WHITE);
                        err::vexit();
                    }
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - unit type \'" << unit_type << "\' is invalid for parameter material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 1.0E5" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //--------------------------------------------------------------------
            test="fmr-field-frequency";
            if(word==test){
                double f=atof(value.c_str());
                // Test for valid range
                if((f>=0.0) && (f<1.0E20)){
                    read_material[super_index].fmr_field_frequency=f;
                    // set local fmr flag
                    sim::local_fmr_field=true;
                    return EXIT_SUCCESS;
                }
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " is outside of valid range 0.0 - 1.0E20" << std::endl;
                    terminaltextcolor(WHITE);
                    err::vexit();
                }
            }
            //------------------------------------------------------------
            test="fmr-field-unit-vector";
            if(word==test){
                // temporary storage container
                std::vector<double> u(3);

                // read values from string
                u=doubles_from_string(value);

                // check size
                if(u.size()!=3){
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " must have three values." << std::endl;
                    terminaltextcolor(WHITE);
                    zlog << zTs() << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " must have three values." << std::endl;
                    return EXIT_FAILURE;
                }

                // Normalise
                double ULength=sqrt(u.at(0)*u.at(0)+u.at(1)*u.at(1)+u.at(2)*u.at(2));

                // Check for correct length unit vector
                if(ULength < 1.0e-9){
                    terminaltextcolor(RED);
                    std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " must be normalisable (possibly all zero)." << std::endl;
                    terminaltextcolor(WHITE);
                    zlog << zTs() << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " must be normalisable (possibly all zero)." << std::endl;
                    return EXIT_FAILURE;
                }
                u.at(0)/=ULength;
                u.at(1)/=ULength;
                u.at(2)/=ULength;

                // Copy field direction to material
                read_material[super_index].fmr_field_unit_vector=u;

                // set local applied field flag
                sim::local_fmr_field=true;

                return EXIT_SUCCESS;

            }
            //--------------------------------------------------------------------
            else
            test="temperature-rescaling-exponent";
            if(word==test){
                double alpha=atof(value.c_str());
                check_for_valid_value(alpha, word, line, prefix, unit, "none", 0.0, 10.0,"material"," 0.0 - 10.0");
                read_material[super_index].temperature_rescaling_alpha=alpha;
                return EXIT_SUCCESS;
            }
            //--------------------------------------------------------------------
            test="temperature-rescaling-curie-temperature";
            if(word==test){
                double Tc=atof(value.c_str());
                check_for_valid_value(Tc, word, line, prefix, unit, "none", 0.0, 10000.0,"material"," 0 - 10000 K");
                read_material[super_index].temperature_rescaling_Tc=Tc;
                return EXIT_SUCCESS;
            }
            //--------------------------------------------------------------------
            test="non-magnetic";
            /*
            integer non-magnetic [0]
            The default value is 0 for all materials. Valid values are
            remove, (blank) [same as remove] and keep.
            Value = keep causes the material to be identified as non magnetic
            with all atoms of this type KEPT in the simulation.
            Value = 1 remove/blank causes the material to be identified as non magnetic
            with all atoms of this type REMOVED from the simulation.
            The atomic positions of non-magnetic atoms are saved separately
            with the usual atomic spin configuration for post processing.
            The default value is false for all materials. Valid values are
            true, false or (blank) [same as true].
            */
            if(word==test){
               test="keep";
               // keep all atoms in simulation (for efficient parallelization)
               if(value==test){
                  // set flag
                  read_material[super_index].non_magnetic = 2;
                  return EXIT_SUCCESS;
               }
               // delete all atoms in simulation (default)
               test="remove";
               if(value==test){
                  // set flag
                  read_material[super_index].non_magnetic = 1;
                  return EXIT_SUCCESS;
               }
               test="";
               if(value==test){
                  // set flag
                  read_material[super_index].non_magnetic = 1;
                  return EXIT_SUCCESS;
               }
               else{
                  terminaltextcolor(RED);
                  std::cerr << "Error on line " << line << " of material file - material[" << super_index+1 << "]:"<< word << " = " << value <<" is not a valid option: remove, keep." << std::endl;
                  terminaltextcolor(WHITE);
                  err::vexit();
               }
            }
            //-------------------------------------------------------------------
            // Call module input parameters
            //-------------------------------------------------------------------
            else if(anisotropy::match_material_parameter(word, value, unit, line, super_index, sub_index, mp::max_materials)) return EXIT_SUCCESS;
            else if(create::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
            else if(dipole::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
            else if(exchange::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
            else if(sim::match_material_parameter(word, value, unit, line, super_index)) return EXIT_SUCCESS;
            else if(st::match_material(word, value, unit, line, super_index)) return EXIT_SUCCESS;
            else if(unitcell::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
            else if(micromagnetic::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
            else if(environment::match_material_parameter(word, value, unit, line, super_index, sub_index)) return EXIT_SUCCESS;
            //--------------------------------------------------------------------
            // keyword not found
            //--------------------------------------------------------------------
            else{
                terminaltextcolor(RED);
                std::cerr << "Error - Unknown control statement '" << line_string << "' on line " << line << " of material file '" << filename_string << "'" << std::endl;
                terminaltextcolor(WHITE);
                zlog << zTs() << "Error - Unknown control statement '" << line_string << " on line " << line << " of material file '" << filename_string << "'" << std::endl;
                return EXIT_FAILURE;
            }

            return EXIT_FAILURE;
    }

}
