/*
// Alireza Rashti
// July 2019
*/

#include "free_data.h"

/* populating the free data part of initial data that we chose ourself */
void populate_free_data(Grid_T *const grid)
{
  /* populate conformal metric and its inverse */
  _gammas(grid);
  
  /* Christoffer symbols made up of conformal metric */
  //_Gamma(grid);
  
  /* Ricci scalar made up of conformal metric _gamma */
  //_R(grid);
  
  /* extrinsic curvature */
  //K(grid);
}

/* populate conformal metric and its inverse */
static void _gammas(Grid_T *const grid)
{
  /* roll off rate at exp(-(r/r0)^4)  */
  const double r0   = 0.5*GetParameterD_E("BH_NS_separation");
  const double M_BH = GetParameterD_E("BH_mass");
  const double a    = GetParameterD_E("BH_dimensionless_spin")*M_BH;
  const double a2   = SQR(a);
  double H,k1,k2,k3;/* in ds^2 = (eta_ij+2*H*ki*kj)dx^i*dx^j */
  /* center of BH */
  const double C_BH = 0.5*GetParameterD_E("BH_NS_separation");
  unsigned p,ijk,nn;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    nn = patch->nn;
    
    /* _gamma */
    double *_gamma_D1D1 = patch->pool[Ind("_gamma_D1D1")]->v;
    double *_gamma_D1D2 = patch->pool[Ind("_gamma_D1D2")]->v;
    double *_gamma_D1D3 = patch->pool[Ind("_gamma_D1D3")]->v;
    double *_gamma_D2D2 = patch->pool[Ind("_gamma_D2D2")]->v;
    double *_gamma_D2D3 = patch->pool[Ind("_gamma_D2D3")]->v;
    double *_gamma_D3D3 = patch->pool[Ind("_gamma_D3D3")]->v;
    
    /* _gamma inverse */
    double *_gammaI_U1U1 = patch->pool[Ind("_gammaI_U1U1")]->v;
    double *_gammaI_U1U2 = patch->pool[Ind("_gammaI_U1U2")]->v;
    double *_gammaI_U1U3 = patch->pool[Ind("_gammaI_U1U3")]->v;
    double *_gammaI_U2U2 = patch->pool[Ind("_gammaI_U2U2")]->v;
    double *_gammaI_U2U3 = patch->pool[Ind("_gammaI_U2U3")]->v;
    double *_gammaI_U3U3 = patch->pool[Ind("_gammaI_U3U3")]->v;
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x   = patch->node[ijk]->x[0];
      double y   = patch->node[ijk]->x[1]-C_BH;
      double z   = patch->node[ijk]->x[2];
      double r2 = SQR(x)+SQR(y)+SQR(z);
      double r  = sqrt(r2);
      double rbar2  = 0.5*(r2-a2+sqrt(SQR(r2-a2)+4*a2*SQR(z)));
      double rbar   = sqrt(rbar2);
      double e   = exp(-pow(r/r0,4));
      
      k1 = (rbar*x+a*y)/(rbar2+a2);
      k2 = (rbar*y-a*x)/(rbar2+a2);
      k3 = z/rbar;
      H  = M_BH*rbar/(rbar2+a2*SQR(k3));
      double C = 2.*H*e;
      double A = 1./(1+C*(SQR(k1)+SQR(k2)+SQR(k3)));
      
      _gamma_D1D1[ijk] = 1.+C*k1*k1;
      _gamma_D1D2[ijk] = C*k1*k2;
      _gamma_D1D3[ijk] = C*k1*k3;
      _gamma_D2D2[ijk] = 1+C*k2*k2;
      _gamma_D2D3[ijk] = C*k2*k3;
      _gamma_D3D3[ijk] = 1+C*k3*k3;
      
      _gammaI_U1U1[ijk] = A*(1+C*(SQR(k2)+SQR(k3)));
      _gammaI_U1U2[ijk] = -A*C*k1*k2;
      _gammaI_U1U3[ijk] = -A*C*k1*k3;
      _gammaI_U2U2[ijk] = A*(1+C*(SQR(k1)+SQR(k3)));
      _gammaI_U2U3[ijk] = -A*C*k2*k3;
      _gammaI_U3U3[ijk] = A*(1+C*(SQR(k1)+SQR(k2)));
      
    }
  }
}

