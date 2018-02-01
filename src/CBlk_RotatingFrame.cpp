
#include "CBlk.h"
#include "parameters.h"
#include "utils.h"


#include <math.h>
#include <ctime>

using namespace std;


//!-----------------------------------------------------------
//! Centrifugal and Coriolis accelerations
//!-----------------------------------------------------------

void CBlock::apply_RotationEffects(PARTICLE_REAL *r, PARTICLE_REAL *w, PARTICLE_REAL *v, D_REAL dt_particle)
{
#if defined use_rotating_frame
  D_REAL a_cor[3] = {0.,0.,0.}; // Coriollis acceleration
  D_REAL v_ang[3] = {0.,0.,0.}; // Angular velocity of the reference frame
  D_REAL a_cfg[3] = {0.,0.,0.}; // Centrifugal acceleration
  
  vec_cross(a_cor, w, v    );
  vec_cross(v_ang, w, r    );
  vec_cross(a_cfg, w, v_ang);

  for (int m=0; m<3; m++) {
    a_cor[m] *= 2.;
    v[m] -= (a_cor[m]+a_cfg[m])*dt_particle;
  }
  
  //! The frame rotates around the rotation axis of the star, not that of the obstacle.
  //! The centrifugal acceleration of the center of the obstacle must be accounted for.
  if (use_Star==true){
    D_REAL OframeO[3]   = {
      0.-O_Frame[0],
      0.-O_Frame[1],
      0.-O_Frame[2]
    }; // Position of the star
    D_REAL v_org[3] = {0.,0.,0.}; // Velocity of obstacle center in the rest frame of the star-obstacle.
    D_REAL a_org[3] = {0.,0.,0.}; // Acceleration of obstacle center in the rest frame of the star-obstacle.
    vec_cross(v_org,w,OframeO);
    vec_cross(a_org,w,v_org);
    for (int m=0; m<3; m++)
      v[m] -= a_org[m]*dt_particle;
  }
#endif
}
