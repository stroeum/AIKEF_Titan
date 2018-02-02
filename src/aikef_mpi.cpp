

#include "defines.h"
#include "CHybrid.h"
#include "utils.h"

#include "parameters.h"
#include "absolute_Globals.h"

#include "SpiceZfc.h"

#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>



using namespace std;

CHybrid* Hybrid;



//!-------------------------------------------------------------//
//! main:
//! - Acces routine of programm.
//! - Hybrid Class is instanced and Hybrid Cycle is performed
//!   within this routine
//!-------------------------------------------------------------//
int main(int argc, char *argv[])
{

#if defined use_rotating_frame
  if (use_Star) {
    D_REAL Position_Obstacle[3] = {0.,0.,0.}; // The obstacle is always at (0,0,0)
    D_REAL r_StarObstacle[3] = {0.,0.,0.};
    D_REAL l = 0., w = 0.;
    for (int m=0; m<3; m++) r_StarObstacle[m] = Position_Obstacle[m]-Position_Star[m];
    l = vec_len(r_StarObstacle);
    w = sqrt((G*M_Star*SI_m0)/pow(l*SI_x0,3.0))*SI_t0;
    w_Frame[0] = 0.;
    w_Frame[1] = 0.;
    w_Frame[2] = w;
  }
#endif

#if defined useSPICE
        //! Load SPICE Kernel at the biginning
	furnsh_c(kernel_name);
#endif
	 
	//! Allocate Instance of Hybrid Class
	Hybrid = new CHybrid;

	//! init MPI, get myRnak and num_processes
   	Hybrid->init_MPI(argc, argv);

	//! read parameter, alloc memory, restore state
   	Hybrid->init();

	//! --------------------------------------------
	//! --- Initialize and collect particle --------
	//! - (some initialization steps are left out)--
	//! --------------------------------------------
	Hybrid->fill_Particle_in_Box();
// 	if(TL==0)
// 	for( int i=0;i<1000; i++ )
// 	Hybrid->inject_obstacle_ions();

	Hybrid->collect_RHOnp1_UIplus_LAM_GAM();
        
	Hybrid->copy_Field(id_rho_n, id_rho_np1);
        
	
#if defined nonadiabatic_gradPE_TERM
        if(TL==0)
        {
            Hybrid->init_PE_Profile();
        }
#endif
	
#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
	Hybrid->init_Neutral_Profile();
#endif	

#if defined(use_ion_production_as_field)
	Hybrid->init_IonProduction_Profile();
#endif	
	
#if defined(use_dust_species_as_field)
	Hybrid->init_Dust_Profile();
#endif		
	
// 	Hybrid->add_force(id_EField);

	Hybrid->negative_particles();


	Hybrid->output_all(TL);


	//! syncronize all process before time is measured,
	//! else conflict when writing state file may appear.
	Hybrid->synchronize_allProcesses_MPI();
	time(&Hybrid->run_start_phys);
	Hybrid->run_start_comp = clock();

	Hybrid->synchronize_allProcesses_MPI();
	Hybrid->global_MPI_lastState_time = MPI_Wtime();
	Hybrid->run_start_mpi = Hybrid->global_MPI_lastState_time;




	//! Time Loop defined as in TB 2004 PhD thesis pp. 27
	for(; TL<TL_MAX;)
	{


		TL++;

		//! measure time for statistics
		Hybrid->measure_time();

		
		//! --- MAIN HYBRID CYCLE ---
		//! 1.)
                if(run_calc_first_E) {
                        Hybrid->calc_first_E();
                }
                
		//! 2.)
		if(run_CAM) {
                        Hybrid->CAM();
                }

		//! 3.)
		if(run_calc_second_E) {
                        Hybrid->calc_second_E();
                }

// 		Hybrid->add_force(id_EField);

		//! 4.)
		if(run_accelerate_Particle) {
                        Hybrid->accelerate_Particle();
                }

		//! 5.) 
		if(run_collect_Ui_minus) {
                        Hybrid->collect_Ui_minus();
                }

		//! 6a.)
		if(run_move_Particle){
                        Hybrid->move_Particle();
                }

		//! 6b)
		if(run_Split_Merge_Particle){
                        Hybrid->Split_Particle();
                        Hybrid->Merge_Particle();
                }
                
		//! scheint zu Abstürzen auf SWARM zu führen:
		//! 20 < 20 ... exiting ... ?!?
// 		Hybrid->check_weight_sorting();

		//! 6b2)
                if(run_negative_Particles){
                        if(negative_particles_in_Simu)
                                Hybrid->negative_particles();		
                }
                
		//! 6b3)
                if(run_delete_too_light_Particles){
                                Hybrid->delete_too_lightweight_particles();		
                }
                
                if(run_collect_Rho_prepare_Recombination_Density){
                        if(!TestParticle_Simulation){
                                Hybrid->collect_RHOnp1_UIplus_LAM_GAM();
                                Hybrid->prepare_Recombination_Density();
                        }
                }
                
                //! 6d.)
                if(run_chemical_Reactions){
                        if(!TestParticle_Simulation) {
                                if(TL!=0) {                
                                        Hybrid->chemical_Reactions();
                                }
                        }
                }
                
		//! 6c.)
		if(run_inject_obstacle_ions){
                        Hybrid->inject_obstacle_ions();
                }	
                
		//! 6f.)
		if(run_resize_pArrays){
                        Hybrid->resize_pArrays();
                }
                
		//! 7.)
		if(run_collect_Rho_calc_Recombination_Density){
                        if(!TestParticle_Simulation){
                                Hybrid->collect_RHOnp1_UIplus_LAM_GAM();
                                Hybrid->calc_Recombination_Density();
                        }
                }
                

		//! 8.)
		if(run_average_Ui_rho_setREZrho){
                        Hybrid->average_Ui_rho_setREZrho();
                }

		//! subcycle plasma and obstacle B if required
		if(run_advanceB){
                        for(INT32 loop=0; loop<num_advance_B_loops; loop++)
                        {
                                //! 9.a)
                                Hybrid->advanceB_Plasma();

                                //! 9.b)
                                Hybrid->advanceB_Obstacle();
                        }
                }


		//! --- DATA ADMINISTRATION ---

		//! A.) Redistribute Blocks to processes
		Hybrid->redistribute_blocks(TL);

		//! B.) Refine Mesh
 		Hybrid->refine_Mesh(TL);

		//! C.) Output Data
		Hybrid->output_all(TL);

		//! D.) Reset Block timing
		Hybrid->reset_block_timing();

		//! E.) Build Statistic
		Hybrid->statistic();

		//! F.) Save State
		if(TL_SAVE_STATE && TL_SAVE_STATE<=TL_MAX) 
		{
			if( !(TL % TL_SAVE_STATE ==0) )
			log_file << " Next save state in "<< TL_SAVE_STATE -TL%TL_SAVE_STATE <<" TL." << endl;
			else
			   Hybrid->save_state();
		}				
		if(TL==TL_MAX)
                         Hybrid->save_state();
                
		//! G.) Resub job
		Hybrid->resubjob();

	}

	if(extended_upstream_boundary_conditions)
                Hybrid->delete_extended_upstream_boundary_conditions();
	
	log_file << " Simulation finished at TL "<<TL<<endl;

	//! finish MPI process
	Hybrid->finalize_MPI();

  #if defined useSPICE
        //unload SPICE Kernel at the end
	 unload_c(kernel_name);
#endif
	 
}

