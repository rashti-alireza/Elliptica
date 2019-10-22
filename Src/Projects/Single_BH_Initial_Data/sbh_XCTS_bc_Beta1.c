#include "sbh_headers.h"
#include "maths_equation_solvings_lib.h"
#include "sbh_XCTS_equations_lib.h"


void *sbh_bc_Beta_U1(void *vp1,void *vp2)
{
  DDM_SCHUR_BC_DECLARE
  unsigned ijk;/* node index */

  /* declaring: */
  GET_FIELD(B0_U1)
  GET_FIELD(_gammaI_U0U1);
  GET_FIELD(_gammaI_U1U1);
  GET_FIELD(_gammaI_U1U2);
  const double M_BH = GetParameterD_E("BH_mass");
  const double a    = GetParameterD_E("BH_X_U2")*M_BH;
  const double a2   = SQR(a);
  

  if (patch->outerB)/* at outer boundary */
  {
  DDM_SCHUR_BC_OPEN
  double x   = patch->node[ijk]->x[0];
  double y   = patch->node[ijk]->x[1];
  double z   = patch->node[ijk]->x[2];
  double r2 = SQR(x)+SQR(y)+SQR(z);
  double rbar2  = 0.5*(r2-a2+sqrt(SQR(r2-a2)+4*a2*SQR(z)));
  double rbar   = sqrt(rbar2);
  double k0 = (rbar*x+a*y)/(rbar2+a2);
  double k1 = (rbar*y-a*x)/(rbar2+a2);
  double k2 = z/rbar;
  double H  = M_BH*rbar/(rbar2+a2*SQR(k2));

  double outerB_F = 
B0_U1[ijk]
-2*H*(_gammaI_U0U1[ijk]*k0+_gammaI_U1U1[ijk]*k1+_gammaI_U1U2[ijk]*k2);

  F[map[ijk]] = outerB_F;

  DDM_SCHUR_BC_CLOSE
  }/* end of if (patch->outerB) */
  else if (patch->innerB)/* at inner boundary */
  {
  DDM_SCHUR_BC_OPEN
  double x   = patch->node[ijk]->x[0];
  double y   = patch->node[ijk]->x[1];
  double z   = patch->node[ijk]->x[2];
  double r2 = SQR(x)+SQR(y)+SQR(z);
  double rbar2  = 0.5*(r2-a2+sqrt(SQR(r2-a2)+4*a2*SQR(z)));
  double rbar   = sqrt(rbar2);
  double k0 = (rbar*x+a*y)/(rbar2+a2);
  double k1 = (rbar*y-a*x)/(rbar2+a2);
  double k2 = z/rbar;
  double H  = M_BH*rbar/(rbar2+a2*SQR(k2));

  double innerB_F = 
B0_U1[ijk]
-2*H*(_gammaI_U0U1[ijk]*k0+_gammaI_U1U1[ijk]*k1+_gammaI_U1U2[ijk]*k2);

  F[map[ijk]] = innerB_F;

  DDM_SCHUR_BC_CLOSE
  }/* end of else if (patch->innerB) */
  return 0;
}
