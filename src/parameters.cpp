/*! @file
 * Main AIKEF Parameter file.
 * 
 * This file is used to define most of the AIKEF settings.
 * All variables used here must also be defined in parameters.h.
 *
 *
 * @par GENERAL REMARKS:
 * - Since virtually all parameters are given in normalized units corresponding to the common Hybrid Code Normalization (See eg. Thorsten Bagdonats (TB) PhD Thesis pp.9), it is not mentioned in the variable description separately.
 * - In case parameters use different units, it is explained in the respective comment
 * - There are two further files, which can be modified by the user:
 *  -# "CHybrid_FieldNames.cpp"  to name Fields as they shall appear in visualisation
 *  -# "CHybrid_IonProfiles.cpp" to define obstacle ion profiles
 
 * @par SOME GENERAL TERMS:
* - Node:                 Point on the numerical Mesh on which physical fields are defined.
*                         (abbriviated by Nd).
* - Ghost Node:    Node at the boundary of the numarical mesh (abbriviated by GN).
*                         The values on GN are not calculated, but copied or fixed.
* - Resolution:    Nodes per Length.
* - Level:        Always meaning the level of refinement. Incrementing
*                         the Level by one results in doubling the resouloution.
* - Cell:                 The cuboid Volume surrounded by 8 Nodes
* - Zone:                 Identical to cell.
* - Block:        Numerical Mesh, consisting of a defined number of Nodes
*                         for each dimension. Each Block is surrounded by a layer
*                         of Ghost Nodes.
*                         (Block is abbriviated by Blk, its Nodes by BlkNds).
* - Root Block:     Block at Level 0 (abbriviated by RB).
*                          Number of Root Blocks remains constant within the entire
*                          simulation run. Root Blocks CANNOT be added or removed after
*                          run has started.
* - Oct:                  One eighth of a block the is build by halfing the block
*                         length in each direction
* - Simulation-Box: Collectivity of all Blocks, defining the entire discretised scenario.
*                          The initial Simulation Box is exclusively build up by Root Blocks.
*                          Even though the RootBlocks are refined during the simulation run,
*                          the length (size/Volume) of the Simulation Box does NOT change, only
*                          resolution is spacially increased.
*                          Initially the number of Nodes in one dimension (eg. X) of the
*                          Simulation Box is given by by RB_X *(BlkNds_X-2).
* - Macro Particle: Physical Ions are represented by Macro Particles (MP). Each MP
*                         represents a certain number of Ions depending on the MP's weight
*                         Do not confuse Weight and (Ion)Mass!
* - Cross Section:  Cut through simulation box data of (CS)
* - Child Block:    Each Block can be refined in ...
*
* @par Parameter are sorted in terms of:
* - 0) Normalization and Background Parameter
* - 1) Geometry Parameter
* - 2) Output Parameter
* - 3) Numerical Parameter
* - 4) Plasma Parameter
* - 5) Solar Wind Background Parameter (removed in current version)
* - 6) Obstacle Parameter
 */
 
#define X         0
#define Y         1
#define Z         2
#define MAGNITUDE 3

#include "defines.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <math.h>

#include "version.h"

/*! Major version number of aikef parameter file.
 *
 * This version number is incremented if changes in the source occur
 * which might break backward compatibility and should match the value in version.h.
 */
#define AIKEF_PARAM_VERSION_MAJOR 6



/*! Check if Major version number of Code and Param. file match.
 * 
 * This ist only to prevent compatibility issues with different versions. 
 * If the versions do not match and no issues are to be expected, the version number can be changeg to the vlaue in version.h to fix this error.
 *
 */
#if AIKEF_PARAM_VERSION_MAJOR!=AIKEF_VERSION_MAJOR
#error "Major AIKEF Version of Parameter File and Source-Code do not match! Please check for compatibility issues and adjust Parameter File version."
#endif

//! To enable NASA NAIF SPICE support uncomment appropriate line in defines.h
//! The kenel_name ist in most cases the name of the SPICE Meta Kernel, that list all other used Kernels
//! for more information please refer to the official SPICE documentation
#if defined useSPICE
        //load SPICE Kernel
	const char* kernel_name="test.txt";
#endif


//!-------------------------------------------------------------//
//!-- 0) Normalization and Solar Wind (sw) Background Parameter //
//!-------------------------------------------------------------//


//! Natural and fundamental constants
//! only defines in parameters.cpp
const D_REAL au = 149597900000.; // astronomical unit

const D_REAL R_Sun = 6.958e+8; // Sun radius, m
const D_REAL M_Sun = 1.98855e30; // Sun mass, kg
const D_REAL T_Sun = 24.47*86400; // Period of rotation of the Sun, s

const D_REAL R_Jupiter = 69911e+3; // Jupiter radius, m
const D_REAL M_Jupiter = 1.898e27; // Jupiter mass, kg
const D_REAL T_Jupiter = 9.92*86400; // Day of Jupiter, s

const D_REAL R_HD179949 = 1.19*R_Sun; // HD179949 radius, m
const D_REAL M_HD179949 = 1.28*M_Sun; // HD179949 mass, kg
const D_REAL T_HD179949 = 9.*86400; // HD179949 period, s

const D_REAL a_HD179949b = 0.045*au; // HD179949b semi-major axis, m
const D_REAL R_HD179949b = 1.05*R_Jupiter; // HD179949b radius, m
const D_REAL M_HD179949b = 0.92*M_Jupiter; // HD179949b mass, kg
const D_REAL T_HD179949b = 3.0925*86400; // Period of rotation of HD179949b, s

const D_REAL pi = M_PI;
const D_REAL mu_0 = 4 * pi * 1.e-7;
const D_REAL m_p = 1.672631e-27;
const D_REAL e = 1.602e-19;
const D_REAL kB = 1.38065e-23;
const D_REAL amu = 1.660538782e-27;
const D_REAL G = 6.67384e-11; // Gravitational constant m3 kg-1 s-2;
#define c_SI 299792458.
#define eps_0 8.854187871e-12

//! parameters of plasma

//! Normalisation quantities
//! Magnetic Field in Telsa
const D_REAL SI_B0 = 5.93e-9;
//! Number Density in Particle per Cubicmeter
const D_REAL SI_n0 = 0.15e+6;
//! Mass of normalization species in amu
const D_REAL SI_m0 = 16*m_p;
//! Alfvén Velocity in Meters per Seconds
const D_REAL SI_v0 = SI_B0/sqrt(mu_0*SI_n0*SI_m0);
//! Inverse Ion Gyration Frequency in Seconds
const D_REAL SI_t0 = SI_m0/(e*SI_B0);
//! Ion Inertia Length in Meters
const D_REAL SI_x0 = SI_v0*SI_t0;
//! Normalizazion acceleration in Meters per Seconds^2
const D_REAL SI_a0 = SI_v0/SI_t0;
//! Normalizazion pressure in pa
const D_REAL SI_p0 = SI_B0*SI_B0/mu_0;
//! Normalizazion magnetic moment in A.m2
const D_REAL SI_MM0 = 4*pi/mu_0*SI_B0*pow(SI_x0,3.0);

//! IMF (Magnetic Background BField) in normalized units
//! Angle of Magnetic Field to Solar Wind Velocity in deg
const D_REAL B_angle = 90;
// const D_REAL B_sw[3] = {cos(B_angle*M_PI/180.), sin(B_angle*M_PI/180.), 0};
const D_REAL B_sw[3] = {0., 0., -1.};

//! Velocity of Solar Wind / upstream plasma
//! in Metre Per Second;
const D_REAL SW_v = 120.0e+3;

//! Alfvénic Mach Number
const D_REAL MA = SW_v/SI_v0;
//! upstream plasma velocity in normalized units

//! Calculation for Plasma Betas
//! Beta = thermalPressure / magnetic Pressure
//! Beta = calcBeta*Temperatur (in Kelvin)
//! calcBeta = SI_n0 * kB /(SI_B0*SI_B0/(2*mu_0))
const D_REAL calcBeta = SI_n0*kB*2*mu_0/(SI_B0*SI_B0);
//! normalized electron mass for electron pressure calculation
const D_REAL mass_electron_norm = 9.10938215e-31/SI_m0;

//! ---------------------------------------------------------------------------
//! hendrik's Enceladus parameters
//! ---------------------------------------------------------------------------
//! use realistic Enceladus plume model or simple test geometry
const bool use_test_scenario=false;//true;
//! which Enceladus flyby
const short flyby = 13; 
//! electron impact ionization frequency in 1/s
const D_REAL Ence_nu_e_total = 1.5e-09; 
//! upstream ion density of all flybys (cf. table 2.2 of phd thesis H.Kriegel)
const D_REAL var_ion_dens[] = {90e+06,60e+06,65e+06,90e+06,55e+06,90e+06,45e+06,60e+06,45e+06,
					75.e+06,65e+06,60e+06,55e+06,65e+06,85e+06,90e+06,105e+06,110e+06,100e+06,70e+06};
  //! SI normalization values
const D_REAL frac_of_species[5]={0.225 + 0.775*use_test_scenario,0.225,0.225,0.225,0.1}; //! O, OH, H2O, H3O, H  NOTE: sum has to be equal to one
const D_REAL average_ion_mass = (frac_of_species[0]*16.+frac_of_species[1]*17.+frac_of_species[2]*18.+frac_of_species[3]*19.+frac_of_species[4]*1)
				*!use_test_scenario + 18.*use_test_scenario;
//! Parameter analytical Plume (if used)
const D_REAL base_density = 1e+15;	//! neutral peak density in SI (m^3)
const D_REAL xi_nD_0 = 1.e+09; 		//! dust peak density in SI (m^3)
const D_REAL H_d = 0.25*252.e+03/SI_x0;
const D_REAL H_d_dust = 948.e+03/SI_x0;
const D_REAL opening_angle_gas = 10;	//! in degrees
const D_REAL opening_angle_dust = 20.;
const D_REAL theta0=0;			//! tilt of plume
const D_REAL phi0= 0;
//! ---------------------------------------------------------------------------


//! ---------------------------------------------------------------------------
//! christoph's comet parameters
//! ---------------------------------------------------------------------------
//! outgassing velocity in Meters per seconds
const D_REAL COM_v0 = 1.e+3;
//! Gas Production Rate of the Comet in 1 per Second
const D_REAL COM_Q[] = {5.e27,0.,0.};
const D_REAL COM_Q_Total = 5.e27;
//! Coll. freq. neutrals/ions in m^3/s
// const D_REAL SI_nu = 1.1e-15;
//! Radius of Cometary Ion Injection
const D_REAL COM_a_SI = 2000.;
//! Raw Density for Neutral Density
//! n_N = COM_RawDensity*1/r^2
const D_REAL COM_RawDensity = 1/(4*M_PI*COM_v0*SI_x0*SI_x0*SI_n0);
//! Calc Com_a
const D_REAL COM_a = COM_a_SI/SI_x0;
//! ---------------------------------------------------------------------------


//!-------------------------------------------------------------//
//!------------- 1) Geometry Parameter -------------------------//
//!-------------------------------------------------------------//


//!-------------------------------------------------------------//
//!-------------- Frame Parameters -----------------------------//
//!-------------------------------------------------------------//
#ifdef use_rotating_frame
D_REAL O_Frame[3] = {0., 0., 0.}; 
//D_REAL w_Frame[3] = {0., 0., 1.e-2}; 
D_REAL w_Frame[3] = {0., 0., 2*pi*SI_t0/T_HD179949b};
#endif
//! Number of Nodes in each Block:
//! - only even values allowed (otw. false interpolation)
//! - Ghost Nodes are included
//!   -> only BlkNds_-2 Nodes are used for physical Field !
//!   NOTE:
//! ------------------------------------------------------------------
//! - USE EQUAL NUMBER IN EACH DIRECITON IN CASE REFINEMENT IS USED, -
//! - ELSE INTERPOLATION AT BOUNDARIES FAILS FOR SOME REASON !!!     -
//!   (maybe child to parent field method ?!?)
//! ------------------------------------------------------------------
const INT32 BlkNds_X = 8;
const INT32 BlkNds_Y = 8;
const INT32 BlkNds_Z = 8;

//! Number of Root Blocks in Szenario
const INT32 RB_X = 16; 
const INT32 RB_Y = 16;
const INT32 RB_Z = 16;

//! Radius of Enceladus
//! (is used instead of R_Obstacle for geometry of neutral cloud since some
//!  test simulations had R_Obstable=0 but boxsize and axis should still
//!  be given in Enceladus radii)
const D_REAL R_Moon = 2575.e+3/SI_x0;

//! Size of simulation box
const D_REAL LX = 10.*R_Moon;
const D_REAL LY = 10.*R_Moon;
const D_REAL LZ = 10.*R_Moon;

//! use periodic Boundaries;
//! every particle leaving the simulation box at any side is inserted
//! at the opposite side of the simulation box. Field derivatives at
//! boundaties are calculated using values from the oposite side of the
//! simulation box.
const bool use_periodic_bounds[3] = {false, false, false};

//! specify number of cells respective boundary set are
//! set to eta value defined below
const INT32 resistive_bound_cells[]  = {3,3,
                                        3,3,
                                        3,3};

//! specify value of eta at boundary
const F_REAL resistive_bound_eta[]  = { 0., 0.,
                                        0., 0.,
                                        0., 0.};

//! alfven wings are strongly reflected at the boundaries
//! therefore the wings have to be damped by using eta
//! over long distances (even longer than one block)
//! also optimal_MPiC is reduced in that block(s)
//! NOTE: works so far only for z-direction
const bool eta_Alfven_Wing_boundary = false;

//! specify distance of respective boundary for eta to go
//! linearly from eta value above to zero
const F_REAL resistive_bound_dist[]  = {0.,0.,
                                        0.,0.,
                                        0.,0.}; 

//! set cometary boundary conditions
//! In case of using the extended upstream boundary conditions 
//! the values of massloading needs to be calculate before starting the time loop
const bool extended_upstream_boundary_conditions = false;

//! General remark on Boundary conditions:
//! 1. entry: x_min-side        2. entry: x_max-side
//! 3. entry: y_min-side        4. entry: y_max-side
//! 5. entry: z_min-side        6. entry: z_max-side
//! in case "use_periodic_bounds" is set true, this parameter
//! will be ignored.


//!-------------------------------------------------------------//
//!-------------- Main Obstacle Parameter ----------------------//
//!-------------------------------------------------------------//
//! specify angular velocity for obstacle
const D_REAL omega_rotating_obs[3] = {0., 0., 0.};

//! Radius of Obstacle
const D_REAL R_Obstacle = R_Moon+1000000/SI_x0;

//! Mass of Obstacle 
//const D_REAL M_Obstacle = (R_Obstacle*SI_x0*SI_v0*SI_v0)/(2.*G*SI_m0);
const D_REAL M_Obstacle = M_HD179949b/SI_m0; 

//! Origin of simulationbox (= Position of Obstacle)
const D_REAL Box_Origin[3] = { 0.501 * LX, 0.501 * LY, 0.501 * LZ};
//const D_REAL Box_Origin[3] = { 5./6. * LX, 0.5 * LY, 0.5 * LZ};

//! Angular velocity of the obstacle. 
const D_REAL w_Obstacle[3] = {0., 0., 0.}; //sqrt((G*M_Obstacle*SI_m0)/(pow(R_Obstacle*SI_x0,3.0)))*SI_t0 

//! Obstacle's intrinsic Dipole Field
//! In order to fix it inside obstacle use
//! a obstacle_core_fraction>0.

//! The variable is_intern_Magnetic_Moment must be set to:
//! (1) TRUE: if Magnetic Moment is from obstacle (then the corresponding field
//! MUST NOT be set at ghostnodes!) 
//! (2) FALSE: if Magnetic Field is extern (such as Saturn's
//! dipole field for Enceladus) (then ghost nodes have to be initialized with 
//! the dipole field)
const bool is_intern_Magnetic_Moment = false;

//! MM is defined such that at the origin MM and BField are parallel.
const D_REAL Magnetic_Moment[3] = { 0., 0., 0*1.e7*pow(R_Obstacle,3.0)/(100*pow(800-R_Obstacle,2.0)) };// see calculation (notebook JAR#14pr639)
//                                       Factor*MMstar*R_pl^3/R_star*(a-R_pl)^2

//! Tilt angle in degree
const D_REAL Magnetic_Moment_Angle[2] = {0., 0.};

//! Position Offset in normalized units
const D_REAL Magnetic_Moment_Offset[3] = {0., 0., 0.};

//!-------------------------------------------------------------//
//!-------------- Secondary Obstacle (Star) Parameter ----------//
//!-------------------------------------------------------------//

//! for two body simulations: activate second obstacle
const bool use_Star = false;

//! Radius of Obstacle
const D_REAL R_Star = 0.51*R_Moon;

//! Position of Second Obstacle
const D_REAL Position_Star[3] = {16.5*R_Moon, 0.,0.}; 
//const D_REAL Position_Star[3] = { .75, 0, 0};

//! Mass of Second Obstacle
//const D_REAL M_Star = 2e30/SI_m0; 
//const D_REAL M_Star = pow(0.25,2.0)*(R_Star*SI_x0*SI_v0*SI_v0)/(2.*G*SI_m0);
const D_REAL M_Star = ( pow(2*M_PI/T_HD179949b,2.0)*pow(a_HD179949b,3.0)/G - M_HD179949b ) /SI_m0; 

//! Angular velocity of the star. 
//! use_synced_rot == true: w_Star=w_Frame, Orbital period of obstacle is synchronized with stellar rotation (the value of w_Star is discarded);
//! use_synced_rot == false: w_Star is given by the value below.
const bool use_synced_rot = false;
const bool T_Star = T_HD179949;
D_REAL w_Star[3] =  {0., 0., 0*2*pi*SI_t0/T_HD179949}; //  sqrt((G*M_Star*SI_m0)/(pow(R_Star*SI_x0,3.0)))*SI_t0

//! Escape velocity: 
//! By default M_Obstacle yields V_esc = .25. 
//! V_sw < V_esc => particles can't escape stellar gravitational pull
//! V_sw < Vo ≈ v_A => particles are subalfvenic at R_Star
//! V_sw > Vo ≈ v_A => particles are superalfvenic at R_Star
//! If all values are normalized 
//! V_sw<1 ⇒ subalfvenic at R_Star
//! V_sw>1 ⇒ superalfvenic at R_Star
const D_REAL V_esc = sqrt((2.*G*M_Star*SI_m0)/(R_Star*SI_x0))/SI_v0; 

// Stellar wind velocity
const D_REAL V_sw[3] = {MA, 0., 0.};
//const D_REAL V_sw[3] = {V_esc*.5, 0., 0.};

//! MM is defined such that at the origin MM and BField are parallel.
const D_REAL beta_sw = 1.e-2;
const D_REAL    m_sw = 1.;
const D_REAL    n_sw = 1.e-6;
const D_REAL   Bh_sw = sqrt(m_sw*n_sw*(V_sw[0]*V_sw[0]+V_sw[1]*V_sw[1]+V_sw[2]*V_sw[2])/beta_sw);
const D_REAL MM_Star = Bh_sw*pow(Position_Star[0]*Position_Star[0] + Position_Star[1]*Position_Star[1] + Position_Star[2]*Position_Star[2],3./2.);
//const D_REAL Magnetic_Moment_Star[3] = { 0., 0., MM_Star};// see calculation (notebook JAR#14p39)
//const D_REAL Magnetic_Moment_Star[3] = { 0., 0., 1.e6};// see calculation (notebook JAR#14p39)
//const D_REAL Magnetic_Moment_Star[3] = { 0., 0., a_HD179949b*a_HD179949b*R_HD179949b/pow(SI_x0,3.0)};// see calculation (notebook JAR#14p39)

const D_REAL Magnetic_Moment_Star[3] = { 0., 0., 0.*1e7};

//! Tilt angle in degree
const D_REAL Magnetic_Moment_Angle_Star[2] = {0., 0.};

//! Position Offset in normalized units
const D_REAL Magnetic_Moment_Offset_Star[3] = {0., 0., 0.};

//!---------------------------------------------------------------------
//!-------------- Particle Boundary Conditions -------------------------
//!---------------------------------------------------------------------


//! 1) First decide, whether particle boundaries have to homogenous
//!      or not
//! 2) In case
//!      a) true  -> further specify whether in- or outflow shall be used
//!      b) false -> inflow & outflow specification will be ignored and
//!         inhomogenous boundaries will be used that are specified
//!         in "utils_boundaries.cpp"
//!NOTE: 7th entry is used for simulation initialization
const bool use_hom_particle_bounds[7] = {
    true,  true,
    true,  true,
    true,  true,
    true
  };
  
  
//! 3) In case homogenous boundaries are used, set use_particle_inflow_bounds to
//!    a) true  ->  inject particle at respective side (inflow)
//!    b) false -> do not inject particle at respective side (outflow)
//! NOTE: Up to now applying outflow boundaries for particle at min-sides
//!       does not work well (Due to assymetry in gather method).
const bool use_particle_inflow_bounds[6] = {
    true, false,
    true, true,
    true, true
  };


//! 3) In case outflow boundaries are used, set use_particle_copy_bounds to
//!    a) true  ->  copy particles from penultimate cell to ultimate
//!    b) false -> to erase all particle in boundary cell
//!        In case respective value in "use_particle_inflow_bounds" is
//!        set true as well, values set here will be ignored.
const bool use_particle_copy_bounds[6] = {
    false, false,
    false, false,
    false, false
  };

//!---------------------------------------------------------------------
//!----------------- Field Boundary Conditions -------------------------
//!---------------------------------------------------------------------
//! 1) First decide, whether BField boundaries have to homogenous
//!      or not
//! 2) In case
//!      a) true  -> further specify whether in or outflow shall be used
//!      b) false -> inflow & outflow specification will be ignored and
//!         inhomogenous boundaries will be used that are specified
//!         in "utils_boundaries.cpp"
//!NOTE: 7th entry is used for simulation initialization
const bool use_hom_B_bounds[7] = {
    true, true,
    true, true,
    true, true,
    true
  };
  

//! set new IMF after TL_new_Bsw (e.g. sector border, magnetopause crossing)
const D_REAL new_B_sw[3] = {0.26, -0.48, -0.84};

//! TL when new BField boundaries shall be set
//! (set zero to switch off)
const INT32 TL_new_Bsw = 100001;

//! in order to avois sharp transition rotate old to
//! new field within finite time
const INT32 num_TL_adapt =  -1;
  
  

//! - in case inhomogenous BField boundaries are used, they can be set
//!   depending on whether plasma is streaming into or out of the boundary
//!   at the respective cell
//! -> it is tested wether the boundary-normal is parallel to the velocity
//!    that is specified in "utils_boundaries.cpp"
//! NOTE:
//! This only makes sense in case inhomogenous particle conditions are
//! specified at the respective boundary
//! NOTE:
//! works for X-MIN/MAX boundaries only
const bool set_BFieldInflow_at_VelocityInflow = false;



//! In case homogenous boundaries are used:
//! For each side of the simulation box has to be specified, which
//! type of boundary conditions shall be used for the EM-Fields:
//! true:  set fields constant to background value (inflow)
//! false: set fields derivative to zero (outflow)
//! in case "use_periodic_bounds" is set true, this parameter
//! will be ignored:

//! specify boundaries for Bfield
const bool use_B_inflow_bounds[6] = {
    true, false,
    false, false,
    false, false
  };

//! specify boundaries for EField
  const bool use_E_inflow_bounds[6] = {
    true, false,
    true, true,
    true, true
  };
  
//! specify boundaries for UField
  const bool use_U_inflow_bounds[6] = {
    true, false,
    true, true,
    true, true
  };
  
#if defined nonadiabatic_gradPE_TERM
  const bool use_PE_inflow_bounds[6] = {
    true, false,
    false, false,
    false, false
  };
#endif
				      
//! In case inhomogenous boundaries are used they must be specied in:
//! "utils_boundaries.cpp"
//! Either
//! 1) An anylytic function can be specified.
//! 2) A one dimensional file can be applied.

//! In case 2) following values have to be specified
//! It is implicitely assumed that extern boundary conditions are 1 dimensional
//! which means that they are either a function of X or Y or Z, but not a
//! function of eg. X and Y at the same time.

//! ----- EBF = Extern Boundary Field -----------------------------------

//! specify how many files shall be read
const INT32 num_EBF = 0;

//! specify how many coordinate values have to be read for each data point
//! (this indicates whether the dataset boundary values are a 1D, 2D or 3D function)
const INT32 num_coords_EBF[] = {3,3,3};

//! specify how many components the fields inside the files have
const INT32 num_comps_EBF[] = {1,1,1};

//! specify name of the respective file
const char file_names_EBF[][RUN_NAME_SIZE] = {"nSW_bounds.txt" , "nCI_bounds", "u_Bounds.txt"};



//!-------------------------------------------------------------//
//!------------- 2) Output Parameter ---------------------------//
//!-------------------------------------------------------------//

//! name of simulation run (used as a prefix for each file)
const char Run_Name[RUN_NAME_SIZE] = "T36";

//! set output path relative to directory of executable
const char data_output_path[RUN_NAME_SIZE] = "../output";


//! End Run when TL_MAX is reached
const INT32 TL_MAX = 200000;


//! SET 0 TO SWITCH OFF RESPECTIVE OUTPUT 

//! vtk    is fileformat for Paraview
//! silo   is fileformat for VisIt
//! native is fileformat for own Visualization

//! zonal Data: Each Zone/Cell is rendered using the color value of the
//!                 lower left Node, resulting in a mosaic-like plot.
//!                 Even though VisIt offers th possibility do visualize
//!                 nodal data in a zonal fashion, some strange interpolatoion is used.
//!                 Hence it is impossible to identify which value is defined on a
//!                 certain node (-> not usefull for e.g. debugging).
//! nodal Data: Each Node is asociated with one color value,
//!                 color of areas inbetween nodes are interpolated.
//! set true for zonal style.
//! set false for nodal style.
bool SILO_STYLE_ZONAL = false;


const INT32  TL_OUTPUT_2D_SILO         =   00000;
const INT32  TL_OUTPUT_2D_NATIVE       =   0;
const INT32  TL_OUTPUT_2D_GNUPLOT      =   0;

const INT32  TL_OUTPUT_3D_SILO         = 100;
const INT32  TL_OUTPUT_3D_NATIVE       =   0;
const INT32  TL_OUTPUT_3D_uniform_grid =  100;
const INT32  TL_OUTPUT_3D_ASCII        =   0;

//! Enable HDF5 GZIP Compresion for SILO-Output
//! Approx. 30% smaller output files but needs slightly more computation time

bool COMPRESS_SILO = true;

//! in case Courant-Criterion is violated (particle faster than one cell per TL)
//! write silo_3D output to enable better debugging
//! NOTE: works only if ESTIMATE_FASTEST_PARTICLE is activated
const bool write_output_if_tofast       = false;

//! write file with total enegy, momentum and mass of ions.
const INT32 TL_OUTPUT_ErgMomMass        = 0;

//! write Line Output

const INT32 TL_OUTPUT_LINEOUT           = 0;
//!  Particle Detctor Output
const INT32 TL_OUTPUT_PARTICLE_DETECTOR = 00000;
//! track path of trajectory and write 1dim field data
const INT32 TL_OUTPUT_TRAJECTORY        = 00000;

const INT32 TL_OUTPUT_PARTICLE_TRACKS   = 0;

const INT32 TL_OUTPUT_GETMESH           = 0;

//!-----------------------------------------------------------------
//!---------- statefile related ------------------------------------
//!-----------------------------------------------------------------

//! save state every xx time level (set to 0 to switch off)
//! NOTE:
//! For each block the responsible MPI process will be stored, hence simulations that
//! were operated on N process must not be re-started with a number of processes
//! smaller than N, yet they can be re-started with a number larger than N.
const INT32 TL_SAVE_STATE  = 2000;

//! OPTOIONS:
//! - LEAVE_UNCHANGED
//! - CONVERT_TO_TRACK
//! - CONVERT_TO_ORDINARY
const INT32 convert_TRACK_PARTICLE_State = LEAVE_UNCHANGED;

//! in general all processes should write their state file at the same time, however
//! some filesystem cannot handle that much data at once. In such case this parameter 
//! should be set true in order to force the processes to write their state file 
//! one after the other. In case HDF5 is used for output this parameter should be set to false,
//! because all filesystem related tasks are handled by the hdf5 api.
const bool serialize_writing_of_StateFile = false;


//! These parameters can be used to switch between the new HDF5 Statefile Format and the old binary format.
//! It is also possible to read old binary files and output new HDF5 Statefiles
const bool write_hdf_state=false;
const bool read_hdf_state=false;

//! Secure State File
//! after a successfull restore the code mv the folder State to StateTLXXXXXX
const bool secure_state_file = false;

//! resubmit of jobfile
//!resubing job if TL_Max is not reached at the end of the Walltime
//!Walltime must be in second
bool resubmit_enable = false;
const INT32 resubmit_walltime =.5*3600;
const char resubmit_jobscript[100] = "msub job.sh";
const INT32 resubmit_security_factor = 60*12;

//!-----------------------------------------------------------------
//!---------- silo output related -------------------------------
//!-----------------------------------------------------------------

//! use factor eg. to scale coordinated from normalized to SI units
//! in visit visualisation.
const F_REAL factor_scale_mesh = 1./R_Moon;

//! Only for 2D output
//! Do not cut through GN Planes in depth direction
//! GN surroundig the Blks are always printed !!!
//! in case set to false all GN have to be
//! update before output.
//! In case set to true update GN are not required,
//! however the cross section might be slightly shifted
//! versus th correct position.
const INT32  avoid_GhostCXS   =  true;

//! Cross Section (CS) to cut (in normed units with respect to origin=(0,0,0))
const F_REAL CROSS_SECTION[3] =   {0., 0., 0.};

//! specify label name for VisIt visualization
const char visit_xlabel[20] = " x";
const char visit_ylabel[20] = " y";
const char visit_zlabel[20] = " z";

//! specify length unit of coordinate axis
const char visit_length_unit[20] = "R_T";


//!-----------------------------------------------------------------
//!---------- trajectory output related -------------------------------
//!-----------------------------------------------------------------

//! Name of trajectory file
//! - Specify trajectory in normalized units relative to Box-Origin
//! - the following example plots the sub-solar line
//!   ending at the box origin
//! - fieldxyz means the data measured by a spacecraft 
//! (time_string)               -30.    0.      0.      (fieldx) (fieldy) (fieldz)
//! (time_string)               -20.    0.      0.      (fieldx) (fieldy) (fieldz)
//! (time_string)               -10.    0.      0.      (fieldx) (fieldy) (fieldz)
//! (time_string)                 0.    0.      0.      (fieldx) (fieldy) (fieldz)


//! specify whether additional time string should be read
const bool do_read_timestring = true;

//! reading the measured field (column 5,6,7) and write into new file
const bool do_read_SpaceCraftField = true;

//! specify how many trajectories shall be traced
const INT32 num_trajectories = 0;

//! specify names of files in which positions of trajectories are strored
const char Trajectory_FileName[][RUN_NAME_SIZE] = { "T36_tiis_1s",
                                                    "Orbit_m6_invY_Z0",
                                                    "CA_MP_BS_m6",
                                                    "CA_MP_BS_m6_invY"};


//! if trajctory file is not given in normalized units, the positions
//! have to be normalized (by ion inertia length)
const D_REAL trajectory_pos_conversion_fact=2575./R_Moon;

//! decide whether to trace trajectory or only to display in 2D out
//! (latter one applyies especially to marks such as MP or BS)
const bool trajectory_2D_out_only[] = {false, false, false, false};

//! compose single file from all files of parallel trajectory output
//! -> will make the code somewhat slower but result in significantly
//!    less files
const bool pack_parallel_output_to_single_file = true;

//! How many fields shall be traced each trajectory
const INT32 num_fields_to_trace = 1;


//! NOTE:
//! In order to trace average fields, use
//! for averaged field 1: id_average_Field1 +0
//! for averaged field 2: id_average_Field1 +1
//! for averaged field 3: id_average_Field1 +2
//! etc ...

//! specify which fields to trace via field id defined in defines.h
const INT32 ID_of_Fields_to_trace[] = {id_BTotal, id_rho_np1};
                                      /*
                                        id_BTotal,
                                        id_average_Field1,
                                        id_BEven,
                                        id_UI_plus,
                                        id_rho_np1,
                                        id_PISpecies1,
                                        id_rotB
                                      */

//! decide whether new file shall be created each time level
//! or output shall be appended to existing file
//! (e.g. for time series)
const bool create_new_file_each_TL = true;

//!-----------------------------------------------------------------
//!---------- Lineout Output related ----------------------------
//!-----------------------------------------------------------------


//! Number of lines for output
const INT32 num_lines = 0;

//! Conversation Factor to normalised Units
//! if lines is not given in normalized units, the positions
//! have to be normalized (by ion inertia length)
const D_REAL lineout_pos_conversion_fact=1.;

//! Names of Lines
//! these will be the names of the output files
//! since these names are also used as variable names for
//! silo ouput, only alphanumeric characters are allowed
//! (especially no "." )
const char LineOut_FileName[][RUN_NAME_SIZE] = {"New_Horizons",
    "New_Horizons2",
    "New_Horizons3",
    "New_Horizons4",
    "New_Horizons5",
    "New_Horizons6"
  };


//! Number of fields in Lineout Output                                              
const INT32 num_fields_to_trace_lineout = 5;

//! specify which fields to trace via field id defined in defines.h
const INT32 ID_of_Fields_to_trace_lineout[] = {id_rhoSpecies1, id_rhoSpecies1+1,id_rhoSpecies1+2,id_UI_Species1,id_UI_plus };
  /*
    ENERGY_SPECTRUM,
    id_PEtotal
    , id_rho_np1
    , id_UI_plus
    , id_divU
    , id_rotB
    , id_average_Field1 + 0 // id_PEtotal
    , id_average_Field1 + 1 // id_rho_np1 
    , id_average_Field1 + 2 // id_UI_plus+0
   */

//! compose single file from all files of parallel trajectory output
//! -> will make the code somewhat slower but result in significantly
//!    less files
const bool pack_parallel_output_to_single_file_lineout = true;

//! decide whether new file shall be created each time level
//! or output shall be appended to existing file
//! (e.g. for time series)
const bool create_new_file_each_TL_lineout = true;


//!-----------------------------------------------------------------
//!---------- Average fields related -------------------------------
//!-----------------------------------------------------------------

//! specify number of fields to average
const INT32 num_average_fields = 0;


//! NOTE --- HOW TO TRACE AVERAGE FIELDS ---
//! When average fields shall be traced specify:
//! id_average_Field1 +0,
//! id_average_Field1 +1
//! id_average_Field1 +2
//! etc.
//! where id_average_Field1 is the first field specified
//! below in "IDs_of_Fields_to_average"


//! specify how many time levels averaging shall
//! start before "TL_FINISH_AVERAGE_FIELDS"
const INT32  num_TL_average_Fields = 75;

//! specify prefix that is added to field name in order to
//! distinguish original vs avereged field
const char average_field_name_prefix[] = "avg";

//! specify when averaging shall be finished
//! (most likely equal to output)
//! set zero to turn off averaging
const INT32  TL_FINISH_AVERAGE_FIELDS = TL_OUTPUT_3D_SILO;


//! eg. if TL_FINISH_AVERAGE_FIELDS = 20 and
//!        num_TL_average_Fields = 10
//!        fields of time levels [10:19], [30:39], ...
//! will be averaged.
//! For this example average-fields will be normalized
//! (ie. devided bynum_TL_average_Fields) at TL=19,39, ...


//! specify fields that shall be averaged via their ID
//! as specified in defines.h
const INT32 IDs_of_Fields_to_average[] = {
  id_rho_np1,
  id_UI_plus, 
  id_BTotal,
  id_EField,
  id_rotB
};
//id_PEtotal


//!-----------------------------------------------------------------
//!---------- uniform_grid Output related -------------------------------
//!-----------------------------------------------------------------

//! some specifications for uniform_grid output
const F_REAL uniform_grid_B0SI = SI_B0;
const F_REAL uniform_grid_n0SI = SI_n0;
const F_REAL uniform_grid_x0SI = SI_x0;
const F_REAL uniform_grid_v0SI = SI_v0;

//! Since the uniform grid cannot handle block adaptive meshes,
//! values of the simulation code's adaptive mesh are interpolated
//! to a catesian mesh for the uniform_grid Tool. The number of mesh nodes
//! of the cartesian mesh has to be specified in "numNds_uniform_gridMesh".
//! The resolution should be choosen in such way that it matches up
//! with the resolution of the adaptive mesh at the highest level
//! of refinement (RB_X*(BlkNds_X-2)*MAX_LEVEL), but one has to be 
//! careful about the size of the resulting file!
const INT32 numNds_uniform_gridMesh[3] = {
  RB_X*(BlkNds_X-2)*(MAX_LEVEL+1) ,  
  RB_Y*(BlkNds_Y-2)*(MAX_LEVEL+1) ,
  RB_Z*(BlkNds_Z-2)*(MAX_LEVEL+1)}; 

//! If the uniform_grid Output should be just a part of the simulation boxsize
//! (e.g. for smaller filesize), output box length could be scaled
const F_REAL uniform_grid_shrink_domain[3] = {1., 1., 1.};
					
//! number of fields which should appear in uniform_grid visualization
const INT32 uniform_grid_num_fields_to_trace = 2;

//! id's of fields which should appear in uniform_grid visualization
const INT32 uniform_grid_ID_of_Fields_to_trace[] = {id_BTotal, id_EField,
  id_rho_np1,
  id_UI_plus,
  id_BTotal,
  id_EField,
  id_average_Field1,
  id_average_Field1+1,
  id_average_Field1+2,
  id_average_Field1+3
};

//! type of field for id's selected above, neccessary to assign correct
//! units to fields in visualization
//! possible values are uniform_grid_TYPE_BFIELD for nT, uniform_grid_TYPE_EFIELD
//! for V/km,
//! uniform_grid_TYPE_VELOCITY_FIELD for km/s and uniform_grid_TYPE_DENSITY_FIELD
//! for 1/cm3
const INT32 uniform_grid_Field_type[] = {uniform_grid_TYPE_BFIELD,uniform_grid_TYPE_EFIELD,
  uniform_grid_TYPE_DENSITY_FIELD,
  uniform_grid_TYPE_VELOCITY_FIELD,
  uniform_grid_TYPE_BFIELD,
  uniform_grid_TYPE_EFIELD,
  uniform_grid_TYPE_DENSITY_FIELD,
  uniform_grid_TYPE_VELOCITY_FIELD,
  uniform_grid_TYPE_BFIELD,
  uniform_grid_TYPE_EFIELD};


//!--------------------------------------------------------------
//!------------ Track Particle (Test Particle Simulation) -------
//!--------------------------------------------------------------

//!--------------------------------------------
//! NOTE:                                     -
//! TRACK_PARTICLE has to be set in defines.h -
//!--------------------------------------------

//! In test particle simulations any function is deactivated except
//! for particle movement and acceleration
const bool TestParticle_Simulation = false;


//! Trace Particles:
const INT32  num_mark_particle_in_species[] = {0, 0, 0, 0,250,0,0,0};
//! NOTE: this number is only exact, if INJECT_PARTICLE_TO_TRACK
//!	  is set to true, else it is rather a proxy

//! volume to mark particles (xmin,xmax,ymin,ymax,zmin,zmax)
const D_REAL mark_particle_volume[6] ={-5.*R_Moon,-3.*R_Moon ,
					-1.*R_Moon, 1.*R_Moon,
					-3.*R_Moon,-1.*R_Moon };
					
//! Trace Particles:
// const INT32  num_mark_particle_in_species[] = {100, 0, 0, 0};

//! if particles shall be tracked in full selfconsistent hybrid
//! simulations, some particles can be "marked"
//! -> choose in function "mark_particle" which to mark
//! set TL_MARK_PARTICLE_TRACKS to -1 in case
//! particles shall not be marked
const INT32  TL_MARK_PARTICLE_TRACKS   = -1;

//! in case a Test Particle Simulation is carried ou in which
//! particle are not initially available, inject particle via
//! function "inject_particle_to_track" right after fields are
//! restored
const bool INJECT_PARTICLE_TO_TRACK = true; 

//! in case a Test Particle Simulation is carried ou in which
//! particle are not initially available, inject particle via
//! function "inject_particle_to_track" right after fields are
//! restored
const INT32 TL_INJECT_PARTICLE_TO_TRACK = -1;

const INT32 TL_CONVERT_PARTICLE_TRACKS_TO_SILO_MESH = 0;

const bool binary_particle_tracks = true;

const INT32 num_tracks_each_group = 100;

//!--------------------------------------------------------------

//!-----------------------------------------------------------------
//!---------- Particle Detector output related ---------------------
//!-----------------------------------------------------------------

//! Number of detector boxes
const INT32 num_particle_detector = 2;
const char Detector_FileName[][RUN_NAME_SIZE] = { "upstream",
                                                    "pickup",
                                                    "detector3"};

//! Detector Box Geometry

const D_REAL detector_box_xmin[] = {-4.5*R_Moon,4.5*R_Moon,0};
const D_REAL detector_box_xmax[] = {-4*R_Moon,5*R_Moon,0};
const D_REAL detector_box_ymin[] = {-0.25*R_Moon,-0.25*R_Moon,0};
const D_REAL detector_box_ymax[] = {0.25*R_Moon,0.25*R_Moon,0};
const D_REAL detector_box_zmin[] = {-0.25*R_Moon,4*R_Moon,0};
const D_REAL detector_box_zmax[] = {0.25*R_Moon,4.5*R_Moon,0};

//!-------------------------------------------------------------//
//!-------------- 3a) Numerical Parameter ----------------------//
//!-------------------------------------------------------------//


//! Numerical time step at Level 0.
//! NOTE:
//! dt should be at least 5 times smaller than
//! Courant Criteria suggest. Otherwise density 
//! lacks at Blk Levl borders occur in flow direction.
//! TODO MAX_LEVEL is set BELOW this parameter!
// const D_REAL  dt =  0.1* RB_X*(BlkNds_X-2)/LX / (SW_v/SI_v0);
const D_REAL  dt =  .0002;//0.00125 * 5. * 4.;

//! specify mesh type
//! 1) UNIFORM
//! 2) STAGGERED
const INT32 mesh_type = UNIFORM;


//! one step:
//! 1) dtB = rot(uxB-...)
//! two step
//! 1)   E = -uxB+...
//! 2) dtB = -rotE
const bool LF_one_step = false;

//! Choose advance_B_Algorithm for Plasma Equations:
//! 0: Leap Frog   (fast)
//! 1: Runge Kutta (a bit more stable than LF)
const INT32 advance_B_Algorithm = 0;

//! advance_B_loop includes:
//! - advance B Plasma 
//!   -> num_sub_cycle times, see defines.h)
//! - advance B Obstacle
//! - divergence cleaning
const INT32 num_advance_B_loops = 3;


//! plasma resistivity
//! could be used instead of smoothing E and B
//! similar numerical diffusion is achieved if
//! smooth_EB = Eta_sw * dt / (cell size)^2 (all in normalized units)
//! (see phd thesis kriegel, pages 71 - 73)
//! rule of thumb: Eta_sw = 0.1 * cell size 
const D_REAL Eta_sw = 0.0;

//! strength of EM-Field smoothing
//! (set value for each refinement level)
const bool several_smooth_in_levels = true; //! if true, smoothing is carried out 2^Level times
const D_REAL smooth_Factor = 0.001;
const D_REAL smooth_E[] = { smooth_Factor * 1., smooth_Factor * 2., smooth_Factor *4. , smooth_Factor * 8.};
const D_REAL smooth_B[] = { smooth_Factor * 1., smooth_Factor * 2., smooth_Factor *4. , smooth_Factor * 8.};

//! number & strength of smoothing the Resistivity profile of the obstacle
//! values are used only if use_resistive_obstacle (section 6) is set to true 
const INT32  num_smooth_eta = 10;
const D_REAL smooth_eta[] = {1./8. *0.5,  1./4. *0.5, 1./2. *0.5, 1. * 0.5,  1.* 0.5};

#if defined(nonadiabatic_gradPE_TERM_smooth_PE)
const D_REAL smooth_PE[] = { smooth_Factor * 1, smooth_Factor * 2., smooth_Factor *4. , smooth_Factor * 8.};
#endif



//! Minimal Charge Density (MCD_J2U) is used in "convert_Ji2Ui" function
//! as an lower limit for rho. Applying MCD_BField results in strongly
//! underestimated velocities.
//! (rho meaning the total density of all ion species)
//! NOTE: this has to be larger than the smallest density different from zero.
//! NOTE: Has to be MCD_BField or comparable for negative species
//!       since a low density cell can have large currents 
const D_REAL MCD_J2U = 1e-4;

//! Minimal Charge Density (MCD_BField) is used 
//! in "magnetice Field" function as an lower limit for rho.
//! (rho meaning the total density of all ion species)
const D_REAL MCD_BField = 0.10;



//!-------------------------------------------------------------//
//!-------------- 3b) Refinement Parameter ---------------------//
//!-------------------------------------------------------------//

//! Number of maximal refinement-levels used
const INT32 MAX_LEVEL = 1;   


//! Apply particle splitting and merging every x TL
//! set 0 to switch off
//! NOTE:
//! if merged particle enter split area "stripes" occur !?!
//! (at least for cold plasma)
const INT32 TL_SPLIT = 1; 
const INT32 TL_MERGE = 1;


//! NOTE:
//! In case to many particle are inside one cell (>8000), the computational
//! speed decreases exponentially hence the code gets AWFULLY slow!
//! Additionally using several thousands particle inside a single cell does
//! not make sense anyway, so "num_merge_each_cell" must be increased.

//! define how many merge/split processes should be carried out in each cell
//! -> the higher oPiC the higer these parameters should be set
const int num_merge_each_cell = 50;
const int num_split_each_cell = 40;

//! ---- DEFINE HOW TO SPLIT/MERGE -------------------------------
//! Decide wheter to split heaviest or most centred particle
//! -> in case most centred is chosen, CoM is conserved by
//!    default an parameters below do not affect splitting
const bool split_heaviest_particle = true;


//! In order to speedup the merging process (which is O(N²)), the
//! best triple of particle is estimated beyond "num_particle_in_MergeList"
//! particle. Since all particle are weight sorted it is ensured, that only
//! the lightest particle on one cell are inside the mergelist
const INT32 num_particle_in_MergeList = 100;//0.5*(optimal_MPiC[0]+optimal_MPiC[1]+optimal_MPiC[2]);//+optimal_MPiC[2]+optimal_MPiC[3]);

//! In order to avoid very high macroparticle numbers in the tail
//! if cell position > merge_tail_distance => merge
const D_REAL merge_tail_distance = -4.;

//! To slightly violate centre of mass conservation drastically
//! reduces noise and avoids numerical artifacts
//! -> in summary this seems to reflect the physics much better
//! -> THIS SHOULD ALWAYS BE SET TO FALSE EXCEPT FOR TESTING !!!!
const bool conserve_CoM_at_merging = false;
const bool conserve_CoM_at_splitting = false;



//! new randomly generated positions may be out of cell
//! -> define how often to retry generating random position
//! -> this parameter will be ignored in case 'conserve_CoM'
//!    is set false
const INT32 num_randomize_position_merge = 200;
const INT32 num_randomize_position_split = 100;

//! level individual activation of splitting/merging
//! in general it should not be required to merge
//! or split in the highest level when "fac_oMPiC_at_LevBoundBlks"
//! is set to 8 or higher
const bool do_merge_in_L[] = {true, true, false, true, true, true};
const bool do_split_in_L[] = {true, true, false, true, true, true};

//! decide whether to split exclusively particle of certain
//! species, eg. do not split planetary Ion species
const bool do_split_in_species[] = {true, true, false, false, false};
const bool do_merge_in_species[] = {true, true, false, false, false};
//! --------------------------------------------------------------





//! ---- GENERAL REFINEMENT OPTIONS -------------------------------


//! increase oMPiC at boundary blocks neighbours in order to
//! OBTAIN smoother transition
const F_REAL fac_oMPiC_at_LevBoundBlks = 4.;


//! decide wheter to smooth or to simply inject moments to parent:
//! Assume L0 is refined to L1 hence particle are in L1.

//! Using the smooth version gives exactly the same density profile for L0,
//! as if particle would have been collected in L0.
//! This only cannot be fullfilled at level boundaries:
//! At minus boundaries for first node inside refined block (SURE ???)
//! At plus boundaries for first node outside refined block
//! (can be seen in visu, follow reason using sketch)
const bool smooth_moments_to_parent = false;

// 
//! This in-accuricy described above appears with and without using gather blocks,
//! so there seems to be no real advantage in using gather blocks.
const bool use_gather_blocks = false;
//! ----------------------------------------------------------------





//! ---- INITIAL HIERARCHICAL REFINEMENT -------------------------------

//! A) CUBOID
//! Specify minimal and maximal Coordinate relative to Box Origin
const INT32 TL_REFINE_STATIC_CUBOID[] = {-1, -1, -1, -1, -1};

//! "left, lower corner"
const D_REAL minX_refBox_of_L[] = {-2.* R_Moon, -6.0* R_Obstacle, -3.5* R_Obstacle};
const D_REAL minY_refBox_of_L[] = {-2.* R_Moon, -2.0* R_Obstacle, -2.0* R_Obstacle};
const D_REAL minZ_refBox_of_L[] = {-2.* R_Moon, -2.0* R_Obstacle, -0.6666* R_Obstacle};
//! "right, upper corner"
const D_REAL maxX_refBox_of_L[] = { 2.* R_Moon,  0.0* R_Obstacle,  0.0* R_Obstacle};
const D_REAL maxY_refBox_of_L[] = { 2.* R_Moon,  2.0* R_Obstacle,  2.0* R_Obstacle};
const D_REAL maxZ_refBox_of_L[] = { 2.* R_Moon,  2.0* R_Obstacle,  0.6666* R_Obstacle};


//! B) SPHERE
//! Specify radius of sphere in which Blks shall be refined

const INT32 TL_REFINE_STATIC_SPHERE[] = {0, 1, 1, -1, -1};


const D_REAL radius_refSphere_of_L[] = { 2.5*R_Moon, 75.0, 50.0, 25.0, 10.0};


//! C) Zylinder along x
const INT32 TL_REFINE_STATIC_ZYLINDER_X[] = {-1, -1, -1, -1, -1};

const D_REAL radiusX_refZylinder_of_L[] = {5.* R_Obstacle, 2.5* R_Obstacle};

const D_REAL minX_refZylinder_of_L[] = { -LX, -LX};
const D_REAL maxX_refZylinder_of_L[] = {  LX,  0.};

//! D) Zylinder along y
const INT32 TL_REFINE_STATIC_ZYLINDER_Y[] = {-1, -1, -1, -1, -1};

const D_REAL radiusY_refZylinder_of_L[] = { 0., 0., 0., 0.};

const D_REAL minY_refZylinder_of_L[] = { 0., 0.};
const D_REAL maxY_refZylinder_of_L[] = { 0., 0.};



//! E) Zylinder along z
const INT32 TL_REFINE_STATIC_ZYLINDER_Z[] = {-1, -1, -1, -1, -1};

const D_REAL radiusZ_refZylinder_of_L[] = { 0., 0.,  3.*R_Obstacle};

const D_REAL minZ_refZylinder_of_L[] = {0., 0.,  -1.3333*R_Obstacle};
const D_REAL maxZ_refZylinder_of_L[] = {0., 0.,  +1.3333*R_Obstacle};

//! forbid removal of initial refined blks
const bool never_remove_static_refinement = true;

//! ----------------------------------------------------------------

//! ---- TIME ADAPTIVE REFINEMENT ----------------------------------
//! general remarks:
//! - refcrit = refinement criteria
//! - any field that has an id (see defines.h) can be used
//!   as refinement criteria

//! REFINEMENT PROCEDURE
//! 1)
//! -  for every oct within one block it is tested, whether a
//!    certain threshhold is exceeded.
//! -> if this is true for any refcrit, the oct becomes "falgged for refinement"

//! 2)
//! - if "force_refine_environment" or "force_refine_entire_block" is set, additional
//!   octs will be falgged for refinement.

//! 3)
//! -  afterwards for any existing oct it is tested, whether a
//!    certain threshhold is undercut.
//! -> if this is true for any refcrit, the oct becomes "falgged for removal"

//! 4)
//! a) octs that are flagged for refinement will be refined
//! b) octs that are flagged for removal will be removed
//! c) octs that are not flagged at all will remain unaffected
//! -  In case refcrit A marks an oct for fefinement while refcrit B marks the
//!    oct for removal, the oct will become refinent


//! Apply Mesh Refinement Algorithm every x TL
//! set 0 to switch off
const INT32 TL_REFINE_MESH =   0;


//! specify how many criterias should be used
const INT32 num_refcrit = 4;


//! specify field IDs for all refcrits
const INT32 refcrit_field_IDs[] = {id_rotB, id_UI_plus, id_UI_plus, id_UI_plus, id_UI_plus};

//! component of respective field that should be used as refcrit:
//! X, Y, Z or MAGNITUDE
//const INT32 refcrit_field_comp[] = {MAGNITUDE,  2, 2, 1, 1};
const INT32 refcrit_field_comp[] = {X,Y,Z,MAGNITUDE};


//! Choose for each refcrit whether it is maximum or minimum based.

//! A)
//! In case maximum based is choosen, for each oct the following condition will be checked:
//! if(oct_average/global_oct_average_maximum > refine_threshold) flag for refinemnt
//! -> the refine_threshold should be a value between 0 and 1 since the value of
//!    "oct_average/global_oct_average_maximum" will be a value between 0 and 1

//! B)
//! In case maximum based is set false, for each oct the following condition will be checked:
//! if(oct_average/global_minimum < refine_threshold) flag for refinemnt.
const bool refcrit_maximum_based[] = {true, false, true, false, true};


//! decide wether oct_average of respective field shall be normalized
//! by its simulation global extremum of wheter the absolute value
//! shall be used
const bool refcrit_normalize_to_global_extremum[] = {false, false, false, false, false};


//! define refine_threshold for each every refcrit and every level:
//! -> each column lists the refine_thresholds for all levels of one refcrit
//! -> each   row  lists the refine_thresholds for all refcrits of one level
const F_REAL refine_threshold[][NUM_CRITERIA] ={{      0.04,    -0.3,    0.3,   -0.6,    0.6},
                                                {      0.12,    -88.0,    88.0,   -113.0,    113.0},
                                                {      0.5,  -100.0,  100.0, -100.0,  100.0},
                                                {     10.00,  -100.0,  100.0, -100.0,  100.0},
                                                {     10.00,  -100.0,  100.0, -100.0,  100.0}};


//! define remove_threshold for each every refcrit and every level:
//! -> each column lists the remove_thresholds for all levels of one refcrit
//! -> each   row  lists the remove_thresholds for all refcrits of one level
const F_REAL remove_threshold[][NUM_CRITERIA] ={{  0.8 *0.04, -0.8* 0.3, 0.8 *0.3, -0.8*0.6, 0.8*0.6},
                                                {  0.8 *0.12, -0.8* 88.0, 0.8 *88.0, -111.8* 3.0, 111.8 *3.0},
                                                {  0.8 *0.5,     100.0,   -100.0,     100.0,   -100.0},
                                                {  0.8 *1.00,     100.0,   -100.0,     100.0,   -100.0},
                                                {  0.8 *1.00,     100.0,   -100.0,     100.0,   -100.0}};




//! In order to avoid discontonuities on level boundaries environment
//! of flagged block s refined
//! -> In case set false lots of noise is introduced
//! NOTE:
//! -> THIS SHOULD ALWAYS BE SET TRUE, OTHERWISE REFINEMENT BOUNDARIES WILL
//!    BE TO CLOSE TO FEATURES OF INTEREST !!!!
//!    (only set false for testing or quick previews)
const bool force_refine_environment =  true;


//! in case one Block is flagged, every Block of
//! the respective child_array will be flagged
//! resulting in standard block AMR
//! -> should always be set FALSE, except for comparison of
//!    "standard block AMR" vs. "hybrid block AMR"
//! NOTE:
//! WILL NOT BE APPLIED FOR INITIAL REFINEMENT
//! BUT ONLY REFINEMENT DURING RUNTIME
const bool force_refine_entire_block  = false;
//! ----------------------------------------------------------------










//!-------------------------------------------------------------//
//!-------------- 3c) Optimization Parameter -------------------//
//!-------------------------------------------------------------//

//! decide whether to order blocks along
//! Space Filling Curve (SFC) or to use the
//! common linear schema
bool use_SFC = true;


//! in order to distribute blocks to intervals of SFC,
//! processing time of each block is recorded.
const INT32 TL_reset_block_timing = 0;

//! in order to obtrain imroved massload distribute
//! blocks among mpi processes
//! Attention: Only some fields will be redistrubuted!!! 
//! Problems with Neutral gas velocity and electron pressure!!!
const INT32 TL_REDISTRIBUTE_BLOCKS = 5000;


//! distribute blocks with some time level offset
//! -> this should make the code redist blocks
//!    some time levels after refinement took
//!    place
//! (for some reason crashes are more likely
//! (when values different from zero are used ???)
//! (-> better set to zero)
const INT32 OFFSET_REDISTRIBUTE_BLOCKS = 0;

//! decide to assign children to respective parent or not
const bool distribute_RB_based = false;

//! force redistribution aof blocks after restore
const bool redistribute_after_restore = true;

//! choose redistibrution criteria:
//! available are:
//! TOTAL_TIME
//! FIELD_TIME
//! BLOCK_NUMBER
//! PARTICLE_TIME
//! PARTICLE_NUMBER

//! below should be used for root block based distribution
//! TOTAL_TIME_INC_CHILDREN
//! FIELD_TIME_INC_CHILDREN
//! BLOCK_NUMBER_INC_CHILDREN
//! PARTICLE_TIME_INC_CHILDREN
//! PARTICLE_NUMBER_INC_CHILDREN
const INT32 distribution_criteria = BLOCK_NUMBER;

//! Synchronisation of MPI Sends and Recieves
//! true = send -> recieve -> wait
//! false = send -> recieve -> wait_all
const bool sync_mpi_send_rec = true;

//! Some processes will operate many blocks with few particles each, other
//! processes will operate few blocks with many particles each.
//! 
//! It may happen that inside block XY are as many particles as in
//! 100 other blocks. Hence one process would deal with block XY
//! where another process would operate all 100 other blocks which leads
//! to a bad load balance at the magnetic field routine.

//! This can be avoided by limiting the maximal number of blocks/each_process
//! to a multiple of the average block number each process (num_blocks/num_processes)

//! (this parameter is ignored in case BLOCK_NUMBER is
//! specified as distribution_criteria)
const F_REAL multiple_of_average = 3.;


//! it is iteratively estimate how many block a given process has to operate,
//! depending on the number of particles inside this block
//! The higher the ireration the better the distribution should be, of course
//! it will sature for very high numbers (>100000)
INT32 max_iteration_redistribute = 400000;


//! show global information (particles moved, collected, moved ...)
//! since all information have to be gathered collecively
//! together which is expensive, this SHOUld not be done to frequently
const INT32  TL_LOGFILE_globalINFO = 1;

//! decide whether each process should write its own log file
//! (in general only useful for debugging)
const bool  show_proc0_logfile_only = true;

//! one particle_array is allocated each cell.
//! -> In case particles enter the cell but the array is to small, a new array is
//!    allocated, particles are copied into it and the old array is deleted
//! -> in order to avoid to frequent re-allocation, min_pArray_Size_factor*num_particle
//!    array elemts are allocated when ever the particle array is to small
const D_REAL min_pArray_Size_factor = 1.1;

//! - In order to avoid to high memory consumption, a small array is allocated in
//!   case num_elements*max_pArray_Size_factor > num_particle
//! - particles are copied and the large array deleted
const D_REAL max_pArray_Size_factor = 1.3;

//! specify how often the it should be checked whether the particle arrays are too large
const INT32 TL_RESIZE_pARRAYS = 1;


//!-------------------------------------------------------------//
//!-------------- 3d) Simulations Configuration ----------------//
//!-------------------------------------------------------------//

//! These variables allows a fast configuration of A.I.K.E.F. for simple test, e.g. test of the ionisation routine
//! This can be used for an automatic test routine. ...
const bool run_calc_first_E = true;
const bool run_CAM = true;
const bool run_calc_second_E = true;
const bool run_accelerate_Particle = true;
const bool run_collect_Ui_minus = true;
const bool run_move_Particle = true;
const bool run_Split_Merge_Particle = true;
const bool run_negative_Particles = false;
const bool run_delete_too_light_Particles = false;
const bool run_inject_obstacle_ions = true;
const bool run_collect_Rho_prepare_Recombination_Density = true;
const bool run_chemical_Reactions = true;
const bool run_resize_pArrays = true;
const bool run_collect_Rho_calc_Recombination_Density = true;
const bool run_average_Ui_rho_setREZrho = true;
const bool run_advanceB = true;


//!-------------------------------------------------------------//
//!-------------- 4) Ion / Ionospheric Parameter ---------------//
//!-------------------------------------------------------------//


//! Total Number of ion and dust species (NOT electrons), has to be at least one.
//! ("Total" meaning Inflow-Ions plus Obstacle-Ions)
const INT32 num_Charged_Species = NUM_PARTICLE_SPECIES;

//! Total number of species that are treated as PARTICLES
//! may be different from num_Charged_Species in case
//! a dust (or ion) species is added to the density from
//! some field
//! NOTE: has to be set in defines.h
const INT32 num_Particle_Species = NUM_PARTICLE_SPECIES;  

//! Number of ion species emmitted from each Box-Boundary.
//! In case set to zero, no plasma will be Initialzed.
//! Still obstacle ions can be inserted, which is a 
//! good (and especially fast) opportunity to find out,
//! how an obstacle ion profile develops.
const INT32 num_Inflow_Species = 2;

//! Indices of inflow species.
//! Exactly "num_Inflow_Species" integers have to be provided.
//! Be carefull not to use the same species for inflow- and obstacle ion
//! at the same time when editing "CHybrid_IonProfiles.cpp".
//! Indices starting at 0 !!!
const INT32 index_Inflow_Species[] = {0,1};


//! define when start to inject respective species
//! - will be ignored in case is no inflow species
//! - set to zero for continues inflow
const D_REAL start_inject_species_each_xx_t0[] = {0 ,0, 0, 0, 0};

//! define duration of injection
//! - will be ignored in case start is zero
const D_REAL duration_inject_species_each_xx_t0[] = {0 ,0, 0, 0, 0};


//! Special Velocity Distribution
//! true: ions which will injected at the boundaries or into empty cell will have a special distribution in velocity space. (c.f. CBlk_Part-Init.cpp -> fill_empty_Cell  for more details)
//! false: ions will have a normal maxwellian velocity distribution
const bool special_Velocity_Distribution[] = {false, false, false, false,false};


//! CONCERNING BRACKETS BELOW:
//! Only the first entry up to num_Particle_Species will be used,
//! all other values will be ignored.

//! optimal_MPiC is the aspired number of Macro Particle for a certain species in each cell.
//! In case of inflow_species, each cell is initialized with this number of
//! particles for the respective species.
//! During simulation the number of particles in a cell will change, so splitting
//! and merging procedures try to adjust the active number to optimal_MPiC.
// const INT32 optimal_MPiC[] = {50, 0, 0, 0,  0, 0};
const INT32 optimal_MPiC[] = {10,10,10, 10,  10, 0};


//! in case particles are split several times, their weight may decrease to values 
//! that give less than 0.01% of the density at this cell. thus delete these particles
//! that have a weight smaller than min_particle_weight
//! switch on/off with run_delete_too_light_Particles
const D_REAL min_particle_weight = 1e-10;

//! Specify Mass, Charge, density, beta/thermal velocity and Obstacle rho for all species:
//!   Mass of 1. denotes one ion Mass of the normalization species (eg. 1 proton mass)
//! Charge of 1. denotes one ion Charge of the normalization species (eg. 1 proton Charge)

//! SPECIES 0 is reference species (!!!)


const D_REAL  Ion_Masses[] = {
    16./16.,
    1./16.,
    2./16.,
    16./16.,
    28./16.,
    0.,
    0,
  };

const D_REAL Ion_Charges[] = { 
    1.,
    1.,
    1.,
    1.,
    1.
  };
			     
//! Plasma Densities for each inflow species.
//! (values at indices of obstacle species will be ignored)
const D_REAL rho_sw[] = {0.1/0.1, 0.01/0.15,  0.,  0.,  0., 0.};



//! Ion plasma betas for each Species:
//! THIS VARIABLE MAY ONLY BE USED, IN CASE "rho_sw" AND "Ion_Masses" ARE
//! SPECIFIED FOR THE RESPECTIVE SPECIES (as mass-density is needed to derive
//! thermal velocity from plasma beta).
//! In case "rho_sw" or "Ion_Masses" are not set for respective
//! species (e.g. in case of obstacle ions), directly specify
//! temperature by using "Ti" variables.
const D_REAL Ion_Betas[]  =  { 
    calcBeta*2000*e/kB,
    calcBeta*200.*e/kB,
    0.,
    0.,
    0
  };

//! allow for different temperature parallel and perpendicular to B
//! TODO currently it is assumed that B0 points in z direction...
//! "temperatures" have to be inserted in Joule
const D_REAL Ti_para[] = {kB*0e6,kB*0*e,0,0,0,0};
const D_REAL Ti_perp[] = {kB*0e6,kB*0*e,0,0,0,0};


//!electron temperature in K
const D_REAL Te = e/kB*100;//1.35*e/kB;      
//! electron plasma betas for each Species.
// const D_REAL Electron_Betas[]  = {calcBeta*Te, 0., 0., 0.};
const D_REAL Electron_Betas[]  = {calcBeta*100*e/kB, calcBeta*100*e/kB, calcBeta*0.01*e/kB,calcBeta*0.01*e/kB,calcBeta*0.01*e/kB};

//! Adiabatic exponent for equation of state
const D_REAL kappa_electron =  2.0;



#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)

	//! Number of neutral species:
	//! NOTE: has to be set in defines.h
	const INT32 num_Neutral_Species = NUM_NEUTRAL_SPECIES;

	//! use analytical expression for neutral density
	const bool analytical_neutral_profile = true; 	
	//! use neutral density profile from file
	const bool neutral_profile_from_file = false;
	//! NOTE parameters above usually exclude each other. However, there may be
	//! rare cases in which both, analytical profile and profile from file are used.
	//! In that case, you should be careful which is set first...

	//! Particle-masses of neutral species (total number of species is defined above)
	// const D_REAL Neutral_Masses[] = { 18*amu/SI_m0, 1., 1., 1. };
	const D_REAL Neutral_Masses[] = { 2.*amu/SI_m0, 16.*amu/SI_m0, 28.*amu/SI_m0, 1., 1. };

	//! Beta of newly ionised electrons
        const D_REAL NewElectron_Betas[]  = {calcBeta*0.01*e/kB, 0.01*calcBeta*e/kB, 0.01*calcBeta*e/kB, calcBeta*180.};
// 	const D_REAL NewElectron_Betas[]  = {0.*calcBeta*100., 0*calcBeta*180., calcBeta*180., calcBeta*180.};

	//! temperature of neutral gas in Kelvin!
	const D_REAL Neutral_Temperature[] = {150,150,150,0};
	
	//! degrees of freedom of neutral species
	//! for example 3 for a pointmass
	//! water = 12 see Page 120
	const D_REAL Neutral_DOF[] = { 3., 3., 3., 3. };

	//! below necessary if neutral field is read from file
	const D_REAL norm_externNeutralDensity[] = {1.,1,1};

	const D_REAL norm_externNeutralVelocity[] = {1.e-03,1,1};

	//! file names for reading extern fields
	const char extern_NeutralField_name[][80] = { "../../../1236.dat",
							"lineY",
							"lineZ"};
							
	
#endif

		
//! -------------- Reactions, collisions, charge exchange and recombination --------------------------------------------------------
							
//! Reaction Rates for Reactions of the type							
//! spec X + NeutralSpecies Y1 -> DestSpec Z + NeutralSpecies Y2
//! ReactionRate[X][][] := Source Ion Species
//! ReactionRate[][Z][] := Destination Ion Species
//! ReactionRate[][][Y] := Neutral Species
//! note that reactions include elastic collisions and charge exchange
//! array entries have to be in normalized units
//! since rates are usually given in cm^3/s, normalize by multiplying with	
D_REAL RateNorm = SI_t0*1e-6*SI_n0;
const D_REAL ReactionRate[][num_Particle_Species][NUM_NEUTRAL_SPECIES] =
  { 
    {         
       {1.94e-10*RateNorm, 7.33e-10*RateNorm,6.82e-10*RateNorm },{ 0.,0,0 },{ 0.,0,0 },{ 0.,0,0 },{ 0.,0,0 }
    },
    {       
       { 0.,0,0 },{1.91e-9*RateNorm, 4.02e-9*RateNorm,3.36e-9*RateNorm },{ 0.,0,0 },{ 0.,0,0 },{ 0.,0,0 }
     },
    {
	{ 0.,0,0 },{ 0.,0,0 },{ RateNorm*1.16e-10,RateNorm*2.76e-9,RateNorm*2.33e-9 },{ 0.,0,0 },{ 0.,0,0 }
    },
    {
	{ 0.,0,0 },{ 0.,0,0 },{ 0.,0,0 },{ RateNorm*1.94e-10,RateNorm*7.33e-10,RateNorm*6.82e-10 },{0, 0,0}
    },
    {
	{ 0.,0,0 },{ 0.,0,0 },{ 0.,0,0 },{ 0.,0,0 },{ RateNorm*1.13e-10,RateNorm*4.72e-10,RateNorm*4.57e-10 }
    }
  };

//! apart from the constant rates in the array above, it is also possible to include a velocity-dependence of the reaction rate
//! the expressions for the velocity dependence have to be provided in 
const bool use_velocity_dependent_reaction_rates = false;
								
//! activate recombination for respective ion species
//! recombination rate is calculated in CBlk_Chemical_Reaction.cpp
//! in function void CBlock::calc_RecombinationAlphaField(void)
const bool Recombination_for_Species[] = {0, 0, 0, 0, 0, 0};

//! check max value of reaction rates
const bool check_max_reaction_rates = true;
//! check max value of p = dt/tau
const bool check_max_reaction_probability = true;

//! --------------------------------------------------------------------------------------------------------------------------------


//!----------------- Photoionisation, eletron impact ionization and stellar wind production ----------------------------------------

//! two ways to describe ionization processes:
//! a) specify ionisation rate if neutral density is set
//	advantage: no global value for production rate has to be known
//! b) specify total ion production rate

//! Characteristic Frequency of Photoionisation in 1/s
//! normalization with SI_t0
// const D_REAL PhotoionisationRate[][NUM_NEUTRAL_SPECIES] =	{ 	 
	// //	   /*Neutralspec 0*/, /*Neutralspec 1*/ ----
	// /* Species 0*/ /*{*/    6.4e-11*SI_t0   ,  /* 0    },*/
	// /* Species 1*/ /*{*/    6.1e-10*SI_t0   ,  /* 0    },*/
	// /* Species 2*/ /*{*/    3.7e-09*SI_t0   ,  /* 0    },*/
	// /* Species 3*/ /*{*/     0		,  /* 0    },*/
	// /* Species 4*/ /*{*/    1.4e-10*SI_t0  /* ,   0    }*/
	// };
const D_REAL PhotoionisationRate[][NUM_NEUTRAL_SPECIES] =	{ 	 
	//	   /*Neutralspec 0*/
	/* Species 0*/ /*{*/ {   0.,0.,0.    }, /* 0    },*/
			{ 0.,0.,0.},
			{ 0.,0.,0.},
			{ 0.,0.,0.},
			{ 0.,0.,5.e-9}
	};
	
//! Characteristic Frequency of Photoionisation in 1/s
const D_REAL COM_nu = 5.88e-7;
	
//! values Enceladus solar maximum	
// const D_REAL PhotoionisationRate[][NUM_NEUTRAL_SPECIES] =	{ 	 
// 	//	   /*Neutralspec 0*/, /*Neutralspec 1*/ ----
// 	/* Species 0*/ {    2.4e-10*SI_t0   ,   0    },
// 	/* Species 1*/ {    1.6e-09*SI_t0   ,   0    },
// 	/* Species 2*/ {    9.1e-09*SI_t0   ,   0    },
// 	/* Species 3*/ {     0		    ,   0    },
// 	/* Species 4*/ {    4.5e-10*SI_t0   ,   0    }
// 	};	
				
	
//! decide if electron ionisation rate depends on electron temperature
//! NOTE not yet implemented for electron pressure equation!
const bool ElectronionisationRate_Te_dependent = false;
//! if false, normalization with SI_t0
//! if true, ElectronionisationRate[][] is multiplied with value of
//! field id_ElectronTemperature (Te in Kelvin)!
//! if expression for ElectronionisationRate is more complex (e.g. sqrt(Te)*const)
//! this has to be changed in CBlk_IonProfile.cpp in insert_ions_from_neutral_profile
// const D_REAL ElectronionisationRate[][NUM_NEUTRAL_SPECIES] =	{ 	 
	// //	                           /*Neutralspec 0*/, /*Neutralspec 1*/ ----
	// /* Species 0*/ /*{*/ Ence_nu_e_total*0.007/1.26*SI_t0   ,  /* 0    },*/
	// /* Species 1*/ /*{*/ Ence_nu_e_total*0.222/1.26*SI_t0   ,   /*0    },*/
	// /* Species 2*/ /*{*/ Ence_nu_e_total*0.958/1.26*SI_t0   ,   /*0    },*/
	// /* Species 3*/ /*{*/ 		0 		    ,   /*0    },*/
	// /* Species 4*/ /*{*/ Ence_nu_e_total*0.0759/1.26*SI_t0  ,   /*0    }*/
	// };	
const D_REAL ElectronionisationRate[][NUM_NEUTRAL_SPECIES] =	{ 	 
	{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}
	};	
	
	
//! this replaces the old and ugly parameter obstacle_MP_weight_each_t0[]
//! entries in ions/s
const D_REAL Global_IonProduction_Rate[][NUM_NEUTRAL_SPECIES] =	{ 	 
	//	     /*Neutralspec 0*/, /*Neutralspec 1*/ ----
	/*Species 0*/ /*{*/    0 ,    /*0.0     },*/
	/*Species 1*/ /*{*/    0 ,    /*0.0     }*/
	};	


//! num macro Particle each TL
//! representing the heavy ions
//! NOTE: if ions are inserted by ionization of neutral profile,
//! the number of inserted ions is adjusted to this value
const INT32 obstacle_MP_num_each_TL[] = {     
    200,
    200,
    200,
    200,
    200,
    0
  };
//! Next tracer currently only work with inject_sphere_source() function. It uses 
const bool use_1st_MP_as_tracer = false;
		

#if defined nonadiabatic_gradPE_COLLISION_TERM

//! collision parameter for electron-neutral-interactions
//! the product of the electron density (ion charge density in normalized units)
//! and the neutral density
//! and the en_coll_param is the factor
//! n_e*m_e /(m_e+m_n) * nu_{e,n}
//! -> see Bachelorthesis (Hendrik Ranocha)
//! since rates are usually given in cm^3/s, normalize by multiplying with RateNorm               
//! implementation requires addition factor m_e /(m_e+m_n)                                          
const D_REAL en_coll_param[] = {                                                                                                                          
    4.427e-11*RateNorm*mass_electron_norm/(mass_electron_norm+Neutral_Masses[0]), 
    0.,
    0.,
    0.
  };


#endif  
  
  //! Inelastic Collisions between electrons and Water (H20)
const bool en_use_inelastic_collisions = false;

//! Input in Kelvin!
const D_REAL en_in_coll_delta_t = 50;

//! Since Inelastic Collisions only with water, specify which neutral species is water.
const INT32 en_water_neutral_species = 0;
						    
#if defined(use_dust_species_as_field)

	//! Number of neutral species:
	//! used as parameters in nonadiabatic_gradPE_TERM computation
	const INT32 num_Dust_Species = 0;

	//! NOTE parameters should be renamed
	//! use analytical expression for dust density, but add this density to density from file
	const bool analytical_dust_plume = false; 	
	//! use analytical expression and no file
	const bool analytical_dust_plume_only = false; 

	//! below necessary if dust field is read from file

	const D_REAL norm_externDustDensity[] = {1.,1,1};

	const D_REAL norm_externDustVelocity[] = {1.e-03,1,1};

	//! file names for reading extern fields
	const char extern_DustField_name[][80] = { "../../../StaubFeld10B.dat",
							"lineY",
							"lineZ"};
							
	//! Full dust density at TL
	const INT32 TL_DUST_MAX = -1;//40000*(1.-0.5*LOWRES)*(1.- 49./50.*VOYAGER);
							
	//! if density corrections according to dusty plasma shall be considered, potential of 
	//! dust cloud may be included in extern rho file
	const bool phi_cloud_from_extern = false;						    
							
	//! smooth dust field
	const INT32 num_smooth_extern = 0;
	const D_REAL smooth_extern[] = {0.0,0.,0.,0.4,0.8};

#endif

	
#if defined(use_ion_production_as_field)
	
	//! ion production field is photoionization frequency * neutral density
	//! due to shadow or absorption of photons (chapman profile), this may be
	//! different from the neutral density profile
	//! the neutral density itself is, however, necessary for the electron pressure
	//! and therefore two different fields are necessary
	
	const INT32 num_ion_prod_fields = 3;

	//! use ion_prod density profile from file
	const bool ion_prod_profile_from_file = true; 
	
	//! file names for reading extern fields
	const char extern_IonProdField_name[][80]= { "../../profiles/Titan/H2_Titan_TA_new",
							"../../profiles/Titan/CH4_Titan_TA_new",
							"../../profiles/Titan/N2_Titan_TA_new"};

	//! calculate the ion production profile of the neutral field
	//! that results from photoabsorption (Chapman profile)	
	const bool calc_ionProd_from_neutral_field = false;

	//! set analytical ion production profile to field
	//! e.g. analytical Chapman profile or modify neutral profile
	//! to consider shadow of moon or planet
	const bool set_analytical_ionProd_profile = false;
							
#endif		


//!--------------------------------------------------------
//!----- extern fields related  --------------------------
//!---------------------------------------------------------

//! This parameter applies to both, neutral and dust fields
const bool serialize_reading_of_extern_field = true;

//! Extern Fields for future usages
//! NOTE extern fields are currently included as neutral fields, ion production fields 
//! and dust fields. However there may be future usage as fields for multi-scale simulations.
//! TODO Horia...

//! specify number of extern density fields for Ion injection
const INT32 num_externRhoVelocityFields = 0;

//! Indices of inflow species.
//! Exactly "num_externRhoVelocityFields" integers have to be provided.
//! Be carefull not to use the same species for inflow,obstacle and extern ions
//! at the same time.
const INT32 index_externRho_Species[5] = {};

//! normalization of extern Density:
//! e.g. Mercury:
//! n0 = 32/cm^3
const D_REAL norm_externRho[] = {32.,1,1};

//! normalization of extern Velocity:
//! e.g. Mercury:
//! B0 = 21nT
//! n0 = 32.e6 1/m^3
//! -> v0 = 80.9728 km/s = 8.09728*10^6 cm/s
const D_REAL norm_externUi[] = {8.09728 *1.e6,1,1};

//! file names for reading extern fields
const char extern_Field_name[][RUN_NAME_SIZE] = { "/velocity_ionization_field.txt",
                                                    "lineY",
                                                    "lineZ"};




//!-------------------------------------------------------------//
//!-------------- 6) Main Obstacle Parameter -------------------//
//!-------------------------------------------------------------//

//! decide if use resistivity inside obstacle
//! details are specified in CBlk_EtaProfiles.cpp
const bool use_resistive_obstacle = true;


//! OBSTACLE MAGNETIC Field calculation.
//! OM1) Switch on/off calculation inside the Obstacle.
//! TODO: does not work with AMR up to now
const bool advance_obstacle_B = true;

//! OM2) choose relaxations parameter for SOR
const D_REAL B_SOR_omega_L[] = {1.15, 1.2, 1.3, 1.};
//{1.2, 1.3, 1.35, 1.4};
//const D_REAL B_SOR_omega_L[] = {1.3, 1.4, 1., 1.};
//! OM3) choose break criteria for SOR
const D_REAL B_SOR_max_error = 1.e-3;

//! OM4) B_SOR_num_cycle_base^L cycles will be performed
//!      in level L each iterations
const INT32 B_SOR_num_cycle_base = 2;

//! OM5) choose maximal number of iterations for SOR
const INT32 B_SOR_max_Iteration = 2000;

//! OM6) calc error every xx TL
const INT32 B_SOR_calc_error_step = 3;


//! DIVERGENCE CLEANER (DC)
//! DC1) Switch on/off BField-divergence cleaning inside the Obstacle.
const bool div_cleaner = true;

//! DC2) choose relaxations parameter for SOR
//const D_REAL DC_omega_L[] = {1.72, 1.78, 1.88, 1.92};
const D_REAL DC_omega_L[] = {1.72, 1.78, 1.75, 1.92};
//! DC3) choose break criteria for SOR
const D_REAL DC_max_error = 1.e-3;

//! DC4) B_SOR_num_cycle_base^L cycles will be performed
//!      in level L each iterations
const INT32 DC_num_cycle_base = 2;

//! DC5) choose maximal number of iterations for SOR
const INT32 DC_max_Iteration = 1000;

//! DC6) calc error every xx TL
const INT32 DC_calc_error_step = 3;


//! The relaxation parameters B_SOR_omega_L[] and DC_omega_L[]
//! crucially determine the convergence speed of the solvers.
//! But they can be determined only empirically for each simulation setup.
//! Thus check which value yields fastest convergence (smallest error)
//! testing values from B_SOR_omega_L-steps/2 in steps of 0.01
//! set to zero to switch optimization off
const INT32 num_optimize_SOR_steps = 0;
//! TL starting the optimization 
//! At TLstart_optimize_B_SOR_omega + num_optimize_SOR_steps, the
//! optimal value will be written in the logfile and this value
//! should be used for next simulation.
const INT32 TLstart_optimize_B_SOR_omega = 10;




//! Switch on/off electric field calculation inside the Obstacle.
//! should usually be switched on. In case of inert moons end high
//! resistivities its more stable to switch it off.
const bool advance_obstacle_E = true;


//! within the core no magnetic field is advanced nor smoothed
//! (-> it will remain the initial field forever)
//! in % of obstacle radius
const D_REAL obstacle_core_fraction = 0.;


//! densities for each ion species to eliminate strong electron pressure 
//! Gradients at Obstacle surface.
//! - inner densities are set INTERNALLY for gradPE calculation
//!   and NOT visible in visualization!
const D_REAL obstacle_rho[] = { MCD_BField,  MCD_BField,  MCD_BField,  MCD_BField,  MCD_BField};


//!-------------------------------------------------------------//
//!-------------- 9) Miscellaneous -----------------------------//
//!-------------------------------------------------------------//

//! Abreviations below were originally defined in absolute_Globals.h.
//! In rare cases this leads to incorrect values, since the order
//! of value assinement is not well defined.
//! Defining these parameters at the end of this file always yields
//! correct assinement.
const INT32 num_root_blocks    = RB_X*RB_Y*RB_Z;
const INT32 num_nodes_in_block = BlkNds_X *BlkNds_Y *BlkNds_Z;


 

