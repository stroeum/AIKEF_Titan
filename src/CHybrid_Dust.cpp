


#include "CHybrid.h"

#include "CBlk.h"

#include "utils.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

#include <math.h>
#define ERROR_SMALL_STATEFILE  "Missmatch in blocksize of state file and active run \n"\
			       "Statefile to SMALL\n"\
			       "Exiting ...."

#define ERROR_LARGE_STATEFILE  "Missmatch in blocksize of state file and active run \n"\
			       "Statefile to LARGE\n"\
			       "Exiting ...."

			       

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "unistd.h"

#include <math.h>


//! time to wait in micro-seconds
#define TIME_TO_WAIT 200000

#if defined(use_dust_species_as_field)

// INT64 num_Blks_in_CrossSection[3];

// ofstream Particle_File;1
ifstream Particle_inFile_Dust;

void CHybrid::init_Dust_Profile(void)
{
	log_file << " INITIALIZING DUST PROFILE ... ";	
	
	INT32 type_to_read = 3;
	
	if(!analytical_dust_plume_only)
	for(INT32 dustSpec=0; dustSpec<num_Dust_Species; dustSpec++)
	read_extern_Profile_uniform_grid(dustSpec,id_density_dustSpecies1+dustSpec,id_velocity_dustSpecies1+dustSpec,type_to_read);
	
	if(analytical_dust_plume)
	set_analytical_dust_plume_extern(id_density_dustSpecies1,id_velocity_dustSpecies1);	
	
	log_file << "done." << endl<<endl;

}



//!------------------------------------------------------------------------
//!- set Plume for NanoDust to extern field 1:
//!------------------------------------------------------------------------
void CHybrid::set_analytical_dust_plume_extern(INT32 id_dustSpecDens, INT32 id_dustSpecVel)
{


	log_file << "  setting analytical dust plume to field " << Field_Name[id_dustSpecDens] << "...       ";
	
	
	//! in case density should not be added to existing density,
	//! field hast to be set to zero
	if(analytical_dust_plume_only)
	set_zero_field(id_dustSpecDens);
	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank==temp_Block->responsible_mpi_process)
			temp_Block->set_analytical_dust_field(id_dustSpecDens,id_dustSpecVel);			
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	
	}
	
	log_file << "done." << endl;

}



bool CHybrid::read_BlockFile_Dust(void)
{
	INT32 level;
	bool has_child[8];
	INT32 mpi_process[8];
	INT32 *root_mpi_process = new INT32[num_root_blocks];

	memset(has_child,        0,               8*sizeof(bool));
	memset(mpi_process,      0,               8*sizeof(INT32));
	memset(root_mpi_process, 0, num_root_blocks*sizeof(INT32));
	
	
	char filename[200];
// 	sprintf(filename,"State/State_Blocks_p%05d",mpi_myRank);

	
	sprintf(filename,"Dust/State_Blocks_p00000");
	
	
	ifstream Block_File;
	Block_File.open(filename, ifstream::binary);
	
	
	
	
	if(!Block_File)
	{

		log_file << "no Block File found." << endl;
		return false;

	}
	
	log_file << " found." << endl;
	//! first sizeof(INT32) BYTES are TL

	//! read time level
	Block_File.read(reinterpret_cast<char*> (&TL), sizeof(INT32));

	//! set TL 
	TL = 1;
	
	//! in order to build an average time of time level, time since prog start
	//! has to be devided by (TL-TL_at_last_restore_state)
	TL_at_last_restore_state = TL;
	
	log_file << " Continuing at TL " << TL << "." << endl << endl;
	
	log_file << " --------------------------------------------------" << endl;
	log_file << " >>>>>>>>>>>>>>>>Restoring Block Tree<<<<<<<<<<<<<<" << endl;
	

	log_file << "  -> restoring root_mpi_process" << endl << endl;

	//! restore mpi_processes of root array:
	Block_File.read(reinterpret_cast<char*> (root_mpi_process), num_root_blocks*sizeof(INT32));

	//! alloc root array using restored mpi processes
	init_RootBlocks_set_mpi_process(root_mpi_process);


	log_file << "  ->  restoring children " << endl;

	for(level=0; level<MAX_LEVEL; level++)
	{

	
		log_file << "  ->  L" << level << endl;

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{


			memset(has_child,   0, 8*sizeof(bool));
			memset(mpi_process, 0, 8*sizeof(INT32));

			Block_File.read(reinterpret_cast<char*>   (has_child), 8*sizeof(bool));
			Block_File.read(reinterpret_cast<char*> (mpi_process), 8*sizeof(INT32));
	
			for(INT32 oct=0; oct<8; oct++)
			 if(has_child[oct])
			 temp_Block->refine_Oct_set_mpi_process(oct, mpi_process[oct]);


			temp_Block = temp_Block->next_Blk_of_BlockList;
		}

	}
	
	
	Block_File.close();
	

	delete[] root_mpi_process;



	synchronize_allProcesses_MPI();





	//! Note:
	//! no MPI reduce required here since all blocks exist on
	//! every process
	
	//! Note:
	//! total_active_Blocks is incresed in refine_Oct(oct), so no additional
	//! summation is required here
	log_file << " Block Tree reconstructed ..." << endl;
	log_file << "  Info (global==local, Blocks exits on each proc)" << endl;
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	log_file << " -> L" << level << ":  " << total_Blocks_L[level] << endl;
	log_file << " total active Blocks: " << total_active_Blocks  << endl;
	log_file << " --------------------------------------------------" << endl << endl;
	log_file << endl;


	

	return true;


}

//!--------------------------------------------------------
//!- restore_state
//!--------------------------------------------------------
bool CHybrid::restore_dust(void)
{




	//! start timer for save state
	time(&last_save_state_time);
	
	particle_restored = false;
	
	log_file << endl << " Searching for state files ...     ";
	
	
	if(!read_BlockFile_Dust()) return false;


	//!-----------------------------------------------------------


	//! set default values for restore assuming number of processes
	//! this run equals number of of processes in statefile
	INT32 processes_to_restore_each_process = 1;

	INT32 restore_start_id = mpi_myRank +0;
	INT32 restore_end_id   = mpi_myRank +1;


	//! count number of involved mpi_processes
	INT32 number_of_processes_in_BlockFile = count_BlockTree_mpi_processes();



	log_file << " number of processes specified in Block File: " << number_of_processes_in_BlockFile << endl;
	log_file << " number of processes specified for this job: "  << mpi_num_processes << endl;
	log_file << endl;


	//! THREE CASES MAY OCCUR:

	//! 1) number_of_processes_in_BlockFile  < 'mpi_num_processes'
	//!	 -> DO RESTORE AS USUAL
	//!	 -> REASSIGNE BLOCKS IN REDISTRIBUTION STEP

	//! 2) number_of_processes_in_BlockFile == 'mpi_num_processes'
	//!	 -> DO RESTORE AS USUAL

	//! 3) number_of_processes_in_BlockFile  > 'mpi_num_processes'
	//!	 -> ASSIGN EACH PROCESS AN INTERVAL OF PROCESSES TO RESTORE
	//!	 -> REORGANIZE ALLOCATED FIELD MEMORY
	//!	 -> RESTORE FILES
	//!	 -> FINALLY SET 'mpi_myRank' TO EACH BLOCK OF INTERVAL

	//! check whether to few processes are available
	if(mpi_num_processes < number_of_processes_in_BlockFile)
	{

		log_file << " number of processes specified in Block File exeeds" << endl;
		log_file << " number of processes specified for this simulation" << endl;

		processes_to_restore_each_process = number_of_processes_in_BlockFile / mpi_num_processes;

		INT32 additional_processes = number_of_processes_in_BlockFile%mpi_num_processes;

		log_file << " Processes_to_restore_each_process: " <<   processes_to_restore_each_process << endl;
		log_file << " The first " << additional_processes  << " processes will restore one additional process." << endl;



		restore_start_id = (mpi_myRank +0) * processes_to_restore_each_process +additional_processes;
		restore_end_id   = (mpi_myRank +1) * processes_to_restore_each_process +additional_processes;

		if(mpi_myRank < additional_processes)
		{

			restore_start_id -= (additional_processes-mpi_myRank-0);
			restore_end_id   -= (additional_processes-mpi_myRank-1);

		}
		log_file << " process " << mpi_myRank << " will handle ["<<restore_start_id<<":"<<restore_end_id<< "[" << endl;



		reorganize_memory_for_assigned_blocks(restore_start_id, restore_end_id);


	}


	//!-----------------------------------------------------------

	//! loop across all files to restore
	for(INT32 restore_process_id = restore_start_id;  restore_process_id<restore_end_id; restore_process_id++)
	{
	
		log_file << " **********************************" << endl;
		log_file << " * Restoring data of process " << restore_process_id <<  endl;
		log_file << " * This Process: " << mpi_myRank <<  endl;
		log_file << " **********************************" << endl;

		//! - carefull, id_BTotal and id_BEven may have the same name
		//! - restore of id_BTotal should not be required
// 		read_Field3D(  id_BEven, restore_process_id);
// 		read_Field3D( id_EField, restore_process_id);
// 
// 		if(read_extern_field_from_state)
// 		{
// 			read_Field3D( id_externRho1, restore_process_id);
// 			
// 			//! only normalize at beginning
// 			if(TL<=TL_read_extern)
// 			multiply_field(id_externRho1,norm_externRho[0]);
// 			
// 			read_Field3D( id_extern_Ui1, restore_process_id);
// 			
// 			//! only normalize at beginning
// 			if(TL<=TL_read_extern)
// 			multiply_field(id_extern_Ui1,norm_externUi[0]);
// 		}
// 		
// 		read_Field3D(  id_PhiDC, restore_process_id);
// 		
// 		for(INT32 average_field=0; average_field<num_average_fields; average_field++)
// 		read_Field3D(id_average_Field1 +average_field, restore_process_id);
			
		read_DustParticle(restore_process_id);
	
	}

	if(mpi_num_processes < number_of_processes_in_BlockFile)
	reset_responsible_mpi_process(number_of_processes_in_BlockFile);

	set_fields();
	
	//! build BTotal rather than restoring it
	//! (done in post mesh refinement)
	build_B_total();
	
	log_file << " done." <<endl;
	synchronize_allProcesses_MPI();
	
	return true;

}


//!------------------------------------------------------------
//!- read_Particle:
//!------------------------------------------------------------
void CHybrid::read_DustParticle(INT32 restorer_mpi_process)
{

	INT32 level;
	CBlock *temp_Block;
	
	char buffer[200];
	num_pArrayReallocs = 0;

	bool missmatch_in_file_size = false;


	for(INT32 species=0; species<num_Particle_Species; species++)
	{
		//! Particle_inFile is global ifstream
		//! (defined at top of this file)
		sprintf(buffer,"Dust/State_PartSpecies%d_p%05d", species, restorer_mpi_process);
		Particle_inFile_Dust.open(buffer, ofstream::binary);
		if(!Particle_inFile_Dust)
		{
			log_file << "no Particle Species "<< species << " File found." << endl;
// 			return;
			//! Do not return, otherwise some processes will
			//! reach collective mpi messages below but
			//! others won't -> deadlock
			//! (may only happen when number of procs was increased after
			//! (restore -> some procs dont have Fields/Particles till then)
		}
		else
		{
		
		
			memset(num_Blocks_transfered_L, 0,(MAX_LEVEL+2)*sizeof(INT64));
			memset(num_particle_processed_L, 0,(MAX_LEVEL+2)*sizeof(INT64));
			
			
		// 	log_file << " restoring Particle of species " << species << " ..." << endl;
			for(level=0; level<=MAX_LEVEL; level++)
			{
				temp_Block = BlockList_of_Lev[level];
				while(temp_Block)
				{
			
					if(restorer_mpi_process == temp_Block->responsible_mpi_process)
					{
			
						//! TODO:
						//! In general it should be sufficient 
						//! only to use top level Blks, however 
						//! not all particles are written down 
						//! when doing so !
						//! -> some refined Blks contain particle ?!?
		// 				if(!temp_Block->child_array)
						{
						
							if(!read_DustParticle_of_Block(species, temp_Block))
							missmatch_in_file_size = true;

							num_Blocks_transfered_L[level]++;
						}
					}
			
					temp_Block = temp_Block->next_Blk_of_BlockList;
				}
			}
			
			//! check whether end of file is reached
			double final;
			Particle_inFile_Dust.read(reinterpret_cast<char*> (&final), sizeof(D_REAL));
		
			if(!Particle_inFile_Dust.eof())
			{
		
				log_file     << ERROR_LARGE_STATEFILE << endl;
				log_file << ERROR_LARGE_STATEFILE << endl;
				missmatch_in_file_size = true;
// 				exit(1);
			}
		
			
			Particle_inFile_Dust.close();

		}


		//! set restored to true, otherwise box will be filled
// 		particle_restored = true;


	
	// 	log_file << "num_pArrayReallocs: "<< num_pArrayReallocs << endl;
		
	
	// 	for(INT32 level=0; level<=MAX_LEVEL; level++)
	// 	log_file << " Top Level Blks in L" << level <<":  "
	// 	         << num_Blocks_transfered_L[level] << endl;
	// 	log_file << endl;
		
	
		log_file << " -----------------------------" << endl;
		log_file << " >>>>>>>>>>Species "<< species << "<<<<<<<<<<<" << endl;
		log_file << "  particle read: " << endl;
	
		//! show information
		INT64 local_info_values[(MAX_LEVEL+1)];
		stringstream info_names[INFO_ARRAY_SIZE];
	
		//! blocks in respective level:
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		{
	
			info_names[level] << " -> L" << level << ":  ";
			local_info_values[level] = num_particle_processed_L[level];
		}
	
		show_information(local_info_values,
				info_names,
				(MAX_LEVEL+1), BUILD_SUM);
		
	
		INT64 total_particle_read = 0;
	
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		total_particle_read +=  local_info_values[level];
		log_file << " Particle in species "<<species<< ": " << total_particle_read  << endl;
		log_file << " -----------------------------" << endl << endl;



		log_file << " checking for missmatch in file size  ..." << endl;
		log_file << " -> (will be marked as 'NaN')" << endl;
		check_for_NaN_MPI(missmatch_in_file_size);
		log_file << " done." << endl;

    	}


	log_file << " *****************************" << endl;
	log_file << " Particle read (all species): " << endl;
	//! show information
	INT64 local_info_values[(MAX_LEVEL+1)];
	stringstream info_names[INFO_ARRAY_SIZE];

	//! particles in respective level:
   	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		info_names[level] << " -> Particle in L" << level << ":  ";
		local_info_values[level] = num_total_particles_in_L[level];
	}

	show_information(local_info_values,
			 info_names,
			 (MAX_LEVEL+1), BUILD_SUM);

	INT64 total_particle_read = 0;

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	total_particle_read +=  local_info_values[level];

	log_file << " -> total Particles: " << total_particle_read << endl;
	log_file << " *****************************" << endl << endl;


	
}


//!------------------------------------------------------------
//!- read_Particle_of_Block:
//!------------------------------------------------------------
bool CHybrid::read_DustParticle_of_Block(INT32 species, CBlock *active_block)
{


	for(INT32 cell=0; cell < num_nodes_in_block; cell++)
	{

		//! check whether end of file is already reached 
		if(Particle_inFile_Dust.eof())
		{
			log_file     << ERROR_SMALL_STATEFILE << endl;
			log_file << ERROR_SMALL_STATEFILE << endl;
			return false;
// 			exit(1);
		}


		Particle_inFile_Dust.read( reinterpret_cast<char*>
				      (active_block->num_MPiC[species]+cell),
				      sizeof(INT32)
				     );
				    
/*		log_file << " num_MPiC[species]+cell " << active_block->num_MPiC[species]+cell <<endl;*/
			
			
		//! increase size/set new size if required
		if(active_block->size_of_pArray[species][cell] < active_block->num_MPiC[species][cell])
		active_block->resize_pArray(species, cell, int( min_pArray_Size_factor *active_block->num_MPiC[species][cell]));


		Particle_inFile_Dust.read(reinterpret_cast<char*>
				    (active_block->pArray[species][cell]),
				     active_block->num_MPiC[species][cell] *(sizeof(particle))
				   );

				   
		//! temprary pointer to particle in list
		particle *active_particle;		   
				    
		for(INT32 part_index=0; part_index<active_block->num_MPiC[species][cell]; part_index++)
		{

			active_particle = active_block->pArray[species][cell] +part_index;
			
			active_particle->v[0] *= 1.e-3;
			active_particle->v[1] *= 1.e-3;
			active_particle->v[2] *= 1.e-3;
		}	
				   
		num_total_particles += 				  active_block->num_MPiC[species][cell];
		num_total_particles_in_L[active_block->RLevel] += active_block->num_MPiC[species][cell];
		num_particle_processed_L[active_block->RLevel] += active_block->num_MPiC[species][cell];
	

		
		
		
	}


	return true;
}


//!-----------------------------------------------------------
//!-----------------------------------------------------------
void CHybrid::add_dust_to_field(void)
{
	
	
	D_REAL *dust_rho_extern, *dust_rho; 
	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{	
	
				for(short dust_field=0; dust_field<num_Dust_Species; dust_field++)
				{	
				
					//! total ion density
					dust_rho_extern = temp_Block->Field_Type[id_density_dustSpecies1+dust_field];					
			
					//! set pointer to respective species
					dust_rho = temp_Block->Field_Type[id_rhoSpecies1 +num_Particle_Species+dust_field];					

				
					for(INT32 node=0; node<num_nodes_in_block; node++)
					{
							if(TL<TL_DUST_MAX)
								dust_rho[node] += dust_rho_extern[node]*TL/TL_DUST_MAX;
							else
								dust_rho[node] += dust_rho_extern[node];
					}
				
				}
					
			}
			
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	
	
	
	
}	
//!-----------------------------------------------------------



//!-----------------------------------------------------------
//!-----------------------------------------------------------
void CHybrid::add_dust_current_to_field(INT32 id_speciesJi)
{
	
	
	D_REAL *dust_rho_extern, *dust_rho, *dust_u_extern_x, *dust_Ji_x; 
	
	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{	
				for(short dust_field=0; dust_field<num_Dust_Species; dust_field++)
				{
					//! total ion density
					dust_rho_extern = temp_Block->Field_Type[id_density_dustSpecies1+dust_field];					
			
					//! set pointers to ion velocity (in scratch field)
					dust_u_extern_x = temp_Block->Field_Type[id_velocity_dustSpecies1+dust_field];
					D_REAL* dust_u_extern_y = dust_u_extern_x + 1*num_nodes_in_block;
					D_REAL* dust_u_extern_z = dust_u_extern_x + 2*num_nodes_in_block;

					//! set pointers to ion velocity (in scratch field)
					dust_Ji_x = temp_Block->Field_Type[id_speciesJi];
					D_REAL* dust_Ji_y = dust_Ji_x + 1*num_nodes_in_block;
					D_REAL* dust_Ji_z = dust_Ji_x + 2*num_nodes_in_block;
					
						
					D_REAL charge = Ion_Charges[num_Particle_Species+dust_field];

					if(id_speciesJi == id_Gam)
					charge = Ion_Charges[num_Particle_Species+dust_field]*Ion_Charges[num_Particle_Species+dust_field]/Ion_Masses[num_Particle_Species+dust_field];


					for(INT32 node=0; node<num_nodes_in_block; node++)
					{
						
						if(TL<TL_DUST_MAX)
						{	
							dust_Ji_x[node] += charge*dust_rho_extern[node]*dust_u_extern_x[node]*TL/TL_DUST_MAX;
							dust_Ji_y[node] += charge*dust_rho_extern[node]*dust_u_extern_y[node]*TL/TL_DUST_MAX;
							dust_Ji_z[node] += charge*dust_rho_extern[node]*dust_u_extern_z[node]*TL/TL_DUST_MAX;
						}	
						else
						{	
							dust_Ji_x[node] += charge*dust_rho_extern[node]*dust_u_extern_x[node];
							dust_Ji_y[node] += charge*dust_rho_extern[node]*dust_u_extern_y[node];
							dust_Ji_z[node] += charge*dust_rho_extern[node]*dust_u_extern_z[node];
						}	
					}
				}	
			}
			
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	
	
	
	
}	


#endif
