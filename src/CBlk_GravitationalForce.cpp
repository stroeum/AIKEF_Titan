 
#include "CBlk.h"
#include "parameters.h"
#include "utils.h"


#include <math.h>
#include <ctime>

using namespace std;


//!-----------------------------------------------------------
//! Gravitational acceleration forces
//!-----------------------------------------------------------

//void CBlock::apply_Gravitation(PARTICLE_REAL *r, PARTICLE_REAL *v, D_REAL dt_particle)
void CBlock::apply_Gravitation(PARTICLE_REAL *r1, PARTICLE_REAL *r2, D_REAL gM, PARTICLE_REAL *v, D_REAL dt_particle)
{
  D_REAL r[3]  = {0.,0.,0.};
  D_REAL r_abs = 0.;
  D_REAL fa[3] = {0.,0.,0.};

  for (int m=0; m<3; m++)
    r[m] = r2[m] -r1[m];
  r_abs = vec_len(r);

  if (r_abs>=0) {
    for(int m=0;m<3;m++) {
      fa[m] = -gM*r[m]/pow(r_abs,3.0);
      v[m] += fa[m]*dt_particle;
    }
  } else {
    printf("Error in CBlock::apply_Gravitation.\n");
    exit(0);
  }
}
