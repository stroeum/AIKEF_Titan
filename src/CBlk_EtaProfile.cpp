

#include "CBlk.h"
#include "parameters.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <math.h>


using namespace std;




//!-------------------------------------------------------------//
//! init_Block: 									//
//!-------------------------------------------------------------//
void CBlock::init_Eta_Profile(void)
{

	//! Obstacle Resistivity
	//! (if Diffusion shall compensate lack of convection within
	//!  the obstacle, ETA has to be larger than V_sw*2*R_Obstacle)
	//!  NOTE: either 
	//!        1) switch on "advance_obstacle_B"
	//!           when using high eta values (>>1.)
	//!           -> Crank Nicholson discretisation for BField
	//!     or 2) "define ETA_TERM_BField"
	//!     	 defines.h for low eta values (~1.)
	//!           -> leap Frog discretisation for BField


	const D_REAL Eta_Obstacle =  25;//10.*vec_len(V_sw)*2*R_Obstacle;

	//! define radius in which Eta_Obstacle shall be used
	const D_REAL R_Eta_Obstacle = 1*R_Moon;
	
	
	//! Either use fermi profile or smoothing for transition from
	//! Eta_Obstacle to background eta (Eta_sw)
	const bool use_eta_fermi = true;//true;
	//! fermi profile for resistivity given by
	//! Eta= Eta_Obstacle *1./(exp((r - R_Eta_Obstacle)*fermi_slope)+1.);
	const D_REAL fermi_slope = 40;
	
	
	//! Decide if a inverse fermi profile should be used
	const bool use_eta_comet = false;
	//! Resistivity in the solar wind at infinity 
	//! the value of Eta_comet_innerComa will be added
	const D_REAL Eta_comet_inf = 3.e-1;
	//! Resistivity in the inner Coma
	const D_REAL Eta_comet_innerComa = 1.e-3;
	//! Determine the radius of low resistivity area
	const D_REAL R_Eta_innerComa = 4.;
	//! slope of change
	const D_REAL fermi_slope_comet = 0.025;

        //! Use smooth framebox function
        const bool use_eta_tanh = false;
        const D_REAL fermi_slope_core_obstacle = 20;
        const D_REAL fermi_slope_obstacle_sw   = 20;

        const D_REAL R1 = obstacle_core_fraction * R_Obstacle;
        const D_REAL R2 = R_Eta_Obstacle;
        const D_REAL eta1 = Eta_Obstacle;
        const D_REAL eta2 = Eta_sw;
        const D_REAL alpha1 = fermi_slope_core_obstacle;
        const D_REAL alpha2 = fermi_slope_obstacle_sw;

	D_REAL r_vec[3] = {0.,0.,0.};
	D_REAL null_vec[3] = {0.,0.,0.};

	
	//! Definition for Star-Planet system
	const D_REAL Eta_core_Obstacle = 0.;
	const D_REAL Eta_core_Star     = 0.;

	const D_REAL Eta_Star 	  = 0.5*Eta_Obstacle;//10.*vec_len(V_sw)*2*R_Star;
  
	//! define radius in which Eta_Obstacle shall be used
	const D_REAL R_Eta_Star     = 0.9*R_Star;

	D_REAL star_core_fraction = obstacle_core_fraction;
	
	//! Star Planet Radii
	D_REAL R[2][2] = {
	  {obstacle_core_fraction*R_Obstacle, R_Eta_Obstacle},
	  {star_core_fraction    *R_Star    , R_Eta_Star    }
	};
	
	D_REAL O[2][3]={
	  {0.,0.,0.},
	  {Position_Star[0],Position_Star[1],Position_Star[2]}
	};
	
        D_REAL OM[2][3]={
	  {0.,0.,0.},
	  {0.,0.,0.}
	};
	
	D_REAL dist_to_O[2] = {0.,0.};

	D_REAL eta[2][3] = {
	  {Eta_core_Obstacle, Eta_Obstacle	, Eta_sw},
	  {Eta_core_Star    , Eta_Star		, Eta_sw}
	};
	for (int i=0; i<2; i++) {
	  for (int j=0; j<3; j++) {
	    eta[i][j] -= Eta_sw/2;
	  }
	}
	D_REAL alpha[2][2] = {
	  {fermi_slope_core_obstacle,fermi_slope_obstacle_sw},
	  {fermi_slope_core_obstacle,fermi_slope_obstacle_sw},
	};
		
	
	
	

	D_REAL* Eta = Field_Type[id_Eta];

	for(INT32 u=0; u< BlkNds_X; u++)
	  for(INT32 v=0; v< BlkNds_Y; v++)
	    for(INT32 w=0; w< BlkNds_Z; w++)
	    {

		INT32 u_v_w  = u*BlkNds_Y*BlkNds_Z 
			      +v*BlkNds_Z 
			      +w;
		

		INT32 cell_indices[3] = {u,v,w};

		if(mesh_type == STAGGERED)
		intern2normedCoords_HIMesh(r_vec, null_vec, cell_indices);
		else
		intern2normedCoords(r_vec, null_vec, cell_indices);


		D_REAL dist_to_orig = vec_len(r_vec);


		//! set Eta_sw as background-eta value
                
		Eta[u_v_w] = Eta_sw;

		//! better use fermi function than define L0 & smooth:
		//! In case only defined on L0, Eta has to be interpolated
		//! to higher levels which lead to "rectangular shape".
	
		//! NOTE:
		//! In case "smooth ETA" version is used, do not forget to
		//! include field_from_parent(Eta) function in refine_Block !
		if(dist_to_orig < R_Eta_Obstacle && !use_eta_fermi)
                    Eta[u_v_w] = Eta_Obstacle;
	
		if(use_eta_fermi)
                    Eta[u_v_w]= Eta_Obstacle *1./(exp((dist_to_orig -R_Eta_Obstacle)*fermi_slope)+1.) + Eta_sw;
	
                if(use_eta_comet)
                    Eta[u_v_w]= Eta_comet_inf *(1.-1/(exp((dist_to_orig -R_Eta_innerComa)*fermi_slope_comet)+1.)) + Eta_comet_innerComa;

		//! eta should be zero inside obstacle core
		if(dist_to_orig < obstacle_core_fraction* R_Obstacle)
                    Eta[u_v_w] = 0.;

		//! Use smooth framebox function
                //!                                   0 < r < obstacle_core_fraction * R_Obstacke => eta = 0;
                //! obstacle_core_fraction * R_Obstacke < r < R_Eta_Obstacle                      => eta = Eta_Obstacle;
                //!         R_Eta_Obstacle < R_Obstacle < r                                       => eta = Eta_sw;
                if(use_eta_tanh){
			if(!use_Star) {
				Eta[u_v_w] = (1+tanh((dist_to_orig-R1)*alpha1))/2 * ( eta2 + (eta1-eta2) * (1-tanh((dist_to_orig-R2)*alpha2))/2  );
			} else {
				Eta[u_v_w] = 0.;
				for(int i=0; i<2; i++){
					for (int n=0; n<3; n++) { 
						OM[i][n] = r_vec[n] - O[i][n];
					}
					dist_to_O[i] = vec_len(OM[i]);
					Eta[u_v_w] += Eta_tanh(dist_to_O[i],R[i],eta[i],alpha[i]);
				}    
			}
		}
	}

	//! set Eta boundaries
	if(eta_Alfven_Wing_boundary)
	for(INT32 u=0; u< BlkNds_X; u++)
	  for(INT32 v=0; v< BlkNds_Y; v++)
	    for(INT32 w=0; w< BlkNds_Z; w++)
	    {

		INT32 u_v_w  = u*BlkNds_Y*BlkNds_Z 
			     +v*BlkNds_Z 
			     +w;
		


		INT32 cell_indices[3] = {u,v,w};
		intern2normedCoords(r_vec, null_vec, cell_indices);
		
		D_REAL eta ;
		
// 		//! upper boundary
// 		eta = resistive_bound_eta[4]/resistive_bound_dist[4] *( r_vec[2] - ( LZ - Box_Origin[2] - resistive_bound_dist[4]));
// 		
// 		if(eta>Eta_sw)
// 		Eta[u_v_w] = eta;
// 		
// 		//! lower boundary
// 		eta = resistive_bound_eta[5]/resistive_bound_dist[5] *( -r_vec[2] + ( -Box_Origin[2] + resistive_bound_dist[5]) );
// 		

		D_REAL xhalfaxis = 6*R_Moon;  
		D_REAL zhalfaxis = 7*R_Moon;  
		D_REAL superellipseexponent = 6;     
		
		D_REAL asym_of_box = Box_Origin[2]-0.5*LZ;
		
		eta = resistive_bound_eta[4]*(pow((r_vec[0]+ Box_Origin[0])/xhalfaxis,superellipseexponent) + pow((r_vec[2]+ asym_of_box)/zhalfaxis,superellipseexponent) - 1)
			/(pow((LX-Box_Origin[0]+ Box_Origin[0])/xhalfaxis,superellipseexponent) + pow((LZ-Box_Origin[2]+ asym_of_box)/zhalfaxis,superellipseexponent) - 1);
	
		if(eta>Eta_sw && vec_len(r_vec)>5.*R_Moon)
		Eta[u_v_w] = eta;
	
	    }


	//! set Eta boundaries
	//!----------- -X boundary ------------------------
	if(is_box_boundary[_Im1_] && resistive_bound_cells[_Im1_])
	 for(INT32 u=0; u<resistive_bound_cells[_Im1_]; u++)
	  for(INT32 v=0; v<BlkNds_Y; v++)
	   for(INT32 w=0; w<BlkNds_Z; w++)
	   {
		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Im1_];

	    }



	//!----------- +X boundary ------------------------
	if(is_box_boundary[_Ip1_] && resistive_bound_cells[_Ip1_])
	 for(INT32 u=BlkNds_X-resistive_bound_cells[_Ip1_]; u<BlkNds_X; u++)
	  for(INT32 v=0;					     v<BlkNds_Y; v++)
	   for(INT32 w=0;				      w<BlkNds_Z; w++)
	   {
		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Ip1_];

	   }


	//!----------- -Y boundary ------------------------
	if(is_box_boundary[_Jm1_]  && resistive_bound_cells[_Jm1_])
	 for(INT32 u=0; u<BlkNds_X; u++)
	  for(INT32 v=0; v<resistive_bound_cells[_Jm1_]; v++)
	   for(INT32 w=0; w<BlkNds_Z; w++)
	   {
		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Jm1_];

	    }

	//!----------- +Y boundary ------------------------
	if(is_box_boundary[_Jp1_]  && resistive_bound_cells[_Jp1_])
	 for(INT32 u=0;					   u<BlkNds_X; u++)
	  for(INT32 v=BlkNds_Y-resistive_bound_cells[_Jp1_]; v<BlkNds_Y; v++)
	   for(INT32 w=0;				      w<BlkNds_Z; w++)
	   {

		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Jp1_];

	   }

	//!----------- -Z boundary ------------------------
	if(is_box_boundary[_Km1_]  && resistive_bound_cells[_Km1_])
	 for(INT32 u=0; u<BlkNds_X; u++)
	  for(INT32 v=0; v<BlkNds_Y; v++)
	   for(INT32 w=0; w<resistive_bound_cells[_Km1_]; w++)
	   {

		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Km1_];

	  }


	//!----------- +Z boundary ------------------------
	if(is_box_boundary[_Kp1_]  && resistive_bound_cells[_Kp1_])
	 for(INT32 u=0;					    u<BlkNds_X; u++)
	  for(INT32 v=0;					     v<BlkNds_Y; v++)
	   for(INT32 w=BlkNds_Z -resistive_bound_cells[_Kp1_]; w<BlkNds_Z; w++)
	   {

		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Kp1_];

	   }





}

D_REAL CBlock::Eta_tanh(D_REAL r, D_REAL R[], D_REAL eta[], D_REAL alpha[])
{
//! Use smooth framebox function
//!    0 < r < R[0] => eta = eta[0];
//! R[0] < r < R[1] => eta = eta[1];
//! R[1] < r        => eta = eta[2];
	//return (eta[2] + (1+tanh((r-R[0])*alpha[0]))/2 * (  eta[2]- eta[0]  +     (eta[1]-eta[2]) * (1-tanh((r-R[1])*alpha[1]))/2));
	return (  eta[2]- eta[0]  +     (eta[1]-eta[2]) * (1-tanh((r-R[1])*alpha[1]))/2);
 
}
