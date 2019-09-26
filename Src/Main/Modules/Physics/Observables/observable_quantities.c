/*
// Alireza Rashti
// September 2019
*/

/* synopsis:
// =========
//
// * initialize observable *
// Observable_T *obs = init_observable(grid);
//
// * specifiy which obeservable *
// obs->quantity = "ADM_momentums" # means ADM momentums
//
// * plan observable *
// plan_observable(obs);# it finds out the related patches, physical metric etc.
//
// * calculate the observable *
// double Px_ADM = obs->ADM_momentum_x(obs);# x component
// double Py_ADM = obs->ADM_momentum_y(obs);# y component
// double Pz_ADM = obs->ADM_momentum_z(obs);# z component
//
// *free*
// free_observable(obs);
*/

#include "observable_quantities.h"

/* initialzing stuct Observable_T */
Observable_T *init_observable(void *grid)
{
  Observable_T *const obs = calloc(1,sizeof(*obs));
  pointerEr(obs);

  obs->grid = grid;
  
  return obs;
}

/* planning observable according to each grid and project */
void plan_observable(Observable_T *const obs)
{
  Grid_T *const grid = obs->grid;
  
  if (strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
  {
    populate_observable_BBN_CS(obs);
  }
  else
    abortEr(NO_OPTION);
  
}

/* populating Observable_T for BBN in cubed spherical case */
static void populate_observable_BBN_CS(Observable_T *const obs)
{
  Grid_T *const grid = obs->grid;
  const unsigned N_outermost = (unsigned) GetParameterI_E("Number_of_Outermost_Split");
  Patch_T **patches = 0,*patch = 0;
  char stem[1000];
  unsigned n,N,ijk,nn;
  
  if (N_outermost == 0)
    abortEr("No outermost patch for integration.\n");
    
  if (strcmp_i(obs->quantity,"ADM_momentums") || strcmp_i(obs->quantity,"ADM_momentum"))
  {
    /* first collect all of the patches required */
    for (n = 0; n < N_outermost; ++n)
    {
      patches = realloc(patches,6*(n+1)*sizeof(*patches));
      pointerEr(patches);
      
      sprintf(stem,"outermost%u_up",n);
      patches[6*n]   = GetPatch(stem,grid);
      
      sprintf(stem,"outermost%u_down",n);
      patches[6*n+1] = GetPatch(stem,grid);
      
      sprintf(stem,"outermost%u_back",n);
      patches[6*n+2] = GetPatch(stem,grid);
      
      sprintf(stem,"outermost%u_front",n);
      patches[6*n+3] = GetPatch(stem,grid);
      
      sprintf(stem,"outermost%u_left",n);
      patches[6*n+4] = GetPatch(stem,grid);
      
      sprintf(stem,"outermost%u_right",n);
      patches[6*n+5] = GetPatch(stem,grid);
    }
    N = 6*n;
    /* alloc memory for each patch */
    obs->ADM = calloc(N,sizeof(*obs->ADM));
    pointerEr(obs->ADM);
    obs->N_ADM = N;
    
    //test
    for (n = 0; n < N; ++n)
      printf("%s added\n",patches[n]->name);
    //end
    
    /* fill ADM struct for each patch */
    for (n = 0; n < N; ++n)
    {
      obs->ADM[n] = calloc(1,sizeof(*obs->ADM[n]));
      pointerEr(obs->ADM[n]);
      patch = patches[n];
      nn    = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);
      
      GET_FIELD(_gamma_D2D2)
      GET_FIELD(_gamma_D0D2)
      GET_FIELD(_gamma_D0D0)
      GET_FIELD(_gamma_D0D1)
      GET_FIELD(_gamma_D1D2)
      GET_FIELD(_gamma_D1D1)
      GET_FIELD(psi);
      
      obs->ADM[n]->patch = patch;
      /* populate metric components */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = SQR(psi[ijk])*SQR(psi[ijk]);
        g00[ijk] = psi4*_gamma_D0D0[ijk];
        g01[ijk] = psi4*_gamma_D0D1[ijk];
        g02[ijk] = psi4*_gamma_D0D2[ijk];
        g11[ijk] = psi4*_gamma_D1D1[ijk];
        g12[ijk] = psi4*_gamma_D1D2[ijk];
        g22[ijk] = psi4*_gamma_D2D2[ijk];
      }
      obs->ADM[n]->g00 = g00;
      obs->ADM[n]->g01 = g01;
      obs->ADM[n]->g02 = g02;
      obs->ADM[n]->g11 = g11;
      obs->ADM[n]->g12 = g12;
      obs->ADM[n]->g22 = g22;
      
      /* outward normal vector is only needed for outermost0 patches
      // where the surface integral carried out */
      if (strstr_i(patch->name,"_outermost0_"))
      {
        obs->ADM[n]->surface_integration_flg = 1;
        obs->ADM[n]->Z_surface = 1;
        obs->ADM[n]->K = 0;

        double *n_U0 = alloc_double(nn);
        double *n_U1 = alloc_double(nn);
        double *n_U2 = alloc_double(nn);
        
        for (ijk = 0; ijk < nn; ++ijk)
        {
          n_U0[ijk] = dq2_dq1(patch,_c_,_x_,ijk);
          n_U1[ijk] = dq2_dq1(patch,_c_,_y_,ijk);
          n_U2[ijk] = dq2_dq1(patch,_c_,_z_,ijk);
          
          /* normalization */
          double psi4 = SQR(psi[ijk])*SQR(psi[ijk]);
          double Norm2 = 
psi4*(_gamma_D0D0[ijk]*pow(n_U0[ijk], 2) + 2.0*_gamma_D0D1[ijk]*
n_U0[ijk]*n_U1[ijk] + 2.0*_gamma_D0D2[ijk]*n_U0[ijk]*n_U2[ijk] +
_gamma_D1D1[ijk]*pow(n_U1[ijk], 2) + 2.0*_gamma_D1D2[ijk]*n_U1[ijk]*
n_U2[ijk] + _gamma_D2D2[ijk]*pow(n_U2[ijk], 2));

          double Norm = sqrt(Norm2);
          
          n_U0[ijk] /= Norm;
          n_U1[ijk] /= Norm;
          n_U2[ijk] /= Norm;
          
        }
        obs->ADM[n]->n_U0 = n_U0;
        obs->ADM[n]->n_U1 = n_U1;
        obs->ADM[n]->n_U2 = n_U2;
      }/* end of if (strstr_i(patch->name,"_outermost0_")) */
      
    }
    
    free(patches);
  }/* end of if (strcmp_i(obs->quantity,"ADM_momentums") || strcmp_i(obs->quantity,"ADM_momentum")) */
  else
    abortEr(NO_OPTION);
  
}

/* free stuct Observable_T */
void free_observable(Observable_T *obs)
{
 if (!obs)
    return;
  
  unsigned i;
    
  for (i = 0; i < obs->N_ADM; ++i)
  {
    _free(obs->ADM[i]->g00);
    _free(obs->ADM[i]->g01);
    _free(obs->ADM[i]->g02);
    _free(obs->ADM[i]->g11);
    _free(obs->ADM[i]->g12);
    _free(obs->ADM[i]->g22);
    _free(obs->ADM[i]->n_U0);
    _free(obs->ADM[i]->n_U1);
    _free(obs->ADM[i]->n_U2);
    
    free(obs->ADM[i]);
    
  }
  _free(obs->ADM);
  free(obs);
}
