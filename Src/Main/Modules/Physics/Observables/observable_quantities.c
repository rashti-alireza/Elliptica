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
// double Px_ADM = obs->Px_ADM(obs);# x component
// double Py_ADM = obs->Py_ADM(obs);# y component
// double Pz_ADM = obs->Pz_ADM(obs);# z component
// double Jz_ADM = obs->Jz_ADM(obs);# z component of angular momentum
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
  else if (strcmp_i(grid->kind,"SBH_CubedSpherical_grid"))
  {
    populate_observable_BBN_CS(obs);/* works as well for SBH */
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
    obs->Px_ADM = ADM_momentum_x_BBN_CS;
    obs->Py_ADM = ADM_momentum_y_BBN_CS;
    obs->Pz_ADM = ADM_momentum_z_BBN_CS;
    obs->Jz_ADM = ADM_angular_momentum_z_BBN_CS;
    
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
    Patch_T *patch = obs->ADM[i]->patch;
    if (patch)
    {
      if (_Ind("ADM_integrand_P_U0") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_P_U0);
        REMOVE_FIELD(ADM_integrand_P_U0);
      }
      if (_Ind("ADM_integrand_P_U1") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_P_U1);
        REMOVE_FIELD(ADM_integrand_P_U1);
      }
      if (_Ind("ADM_integrand_P_U2") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_P_U2);
        REMOVE_FIELD(ADM_integrand_P_U2);
      }
      if (_Ind("ADM_integrand_G_U0") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_G_U0);
        REMOVE_FIELD(ADM_integrand_G_U0);
      }
      if (_Ind("ADM_integrand_G_U1") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_G_U1);
        REMOVE_FIELD(ADM_integrand_G_U1);
      }
      if (_Ind("ADM_integrand_G_U2") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_G_U2);
        REMOVE_FIELD(ADM_integrand_G_U2);
      }
    }
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

/* calculating ADM momentum in x component */
static double ADM_momentum_x_BBN_CS(Observable_T *const obs)
{
  double Px = 0;
  const unsigned N = obs->N_ADM;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  populate_ADM_momentums_integrand_PdS_GdV(obs);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = obs->ADM[p]->patch;
    
    if (!obs->ADM[p]->surface_integration_flg)
      continue;
    
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_P_U0")];
    I->g00 = obs->ADM[p]->g00;
    I->g01 = obs->ADM[p]->g01;
    I->g02 = obs->ADM[p]->g02;
    I->g11 = obs->ADM[p]->g11;
    I->g12 = obs->ADM[p]->g12;
    I->g22 = obs->ADM[p]->g22;
    
    I->Spectral->Z_surface = obs->ADM[p]->Z_surface;
    I->Spectral->K         = obs->ADM[p]->K;
    
    plan_integration(I);
    Px += execute_integration(I);
    free_integration(I);
  }
  
  /* volume integration */
  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = obs->ADM[p]->patch;
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_G_U0")];
    I->g00 = obs->ADM[p]->g00;
    I->g01 = obs->ADM[p]->g01;
    I->g02 = obs->ADM[p]->g02;
    I->g11 = obs->ADM[p]->g11;
    I->g12 = obs->ADM[p]->g12;
    I->g22 = obs->ADM[p]->g22;
    
    plan_integration(I);
    Px -= execute_integration(I);
    free_integration(I);
  }
  
  Px /= (8*M_PI);
  return Px;
}

/* calculating ADM momentum in x component */
static double ADM_momentum_y_BBN_CS(Observable_T *const obs)
{
  double Py = 0;
  const unsigned N = obs->N_ADM;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  populate_ADM_momentums_integrand_PdS_GdV(obs);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = obs->ADM[p]->patch;
    
    if (!obs->ADM[p]->surface_integration_flg)
      continue;
    
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_P_U1")];
    I->g00 = obs->ADM[p]->g00;
    I->g01 = obs->ADM[p]->g01;
    I->g02 = obs->ADM[p]->g02;
    I->g11 = obs->ADM[p]->g11;
    I->g12 = obs->ADM[p]->g12;
    I->g22 = obs->ADM[p]->g22;
    
    I->Spectral->Z_surface = obs->ADM[p]->Z_surface;
    I->Spectral->K         = obs->ADM[p]->K;
    
    plan_integration(I);
    Py += execute_integration(I);
    free_integration(I);
  }
  
  /* volume integration */
  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = obs->ADM[p]->patch;
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_G_U1")];
    I->g00 = obs->ADM[p]->g00;
    I->g01 = obs->ADM[p]->g01;
    I->g02 = obs->ADM[p]->g02;
    I->g11 = obs->ADM[p]->g11;
    I->g12 = obs->ADM[p]->g12;
    I->g22 = obs->ADM[p]->g22;
    
    plan_integration(I);
    Py -= execute_integration(I);
    free_integration(I);
  }
  
  Py /= (8*M_PI);
  return Py;
}

/* calculating ADM momentum in x component */
static double ADM_momentum_z_BBN_CS(Observable_T *const obs)
{
  double Pz = 0;
  const unsigned N = obs->N_ADM;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  populate_ADM_momentums_integrand_PdS_GdV(obs);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = obs->ADM[p]->patch;
    
    if (!obs->ADM[p]->surface_integration_flg)
      continue;
    
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_P_U2")];
    I->g00 = obs->ADM[p]->g00;
    I->g01 = obs->ADM[p]->g01;
    I->g02 = obs->ADM[p]->g02;
    I->g11 = obs->ADM[p]->g11;
    I->g12 = obs->ADM[p]->g12;
    I->g22 = obs->ADM[p]->g22;
    
    I->Spectral->Z_surface = obs->ADM[p]->Z_surface;
    I->Spectral->K         = obs->ADM[p]->K;
    
    plan_integration(I);
    Pz += execute_integration(I);
    free_integration(I);
  }
  
  /* volume integration */
  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = obs->ADM[p]->patch;
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_G_U2")];
    I->g00 = obs->ADM[p]->g00;
    I->g01 = obs->ADM[p]->g01;
    I->g02 = obs->ADM[p]->g02;
    I->g11 = obs->ADM[p]->g11;
    I->g12 = obs->ADM[p]->g12;
    I->g22 = obs->ADM[p]->g22;
    
    plan_integration(I);
    Pz -= execute_integration(I);
    free_integration(I);
  }
  
  Pz /= (8*M_PI);
  return Pz;
}

/* calculating ADM angular momentum in z component */
static double ADM_angular_momentum_z_BBN_CS(Observable_T *const obs)
{
  double Jz = 0;
  const unsigned N = obs->N_ADM;
  double y_CM = 0;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  populate_ADM_momentums_integrand_PdS_GdV(obs);
  
  if (get_parameter("y_CM"))
    y_CM = GetParameterD_E("y_CM");
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = obs->ADM[p]->patch;
    unsigned nn = patch->nn;
    unsigned ijk;
    
    if (!obs->ADM[p]->surface_integration_flg)
      continue;
    
    Field_T *Py = patch->pool[Ind("ADM_integrand_P_U1")];
    Field_T *Px = patch->pool[Ind("ADM_integrand_P_U0")];
    ADD_FIELD(xPy_adm_integrand);
    ADD_FIELD(yPx_adm_integrand);
    DECLARE_FIELD(xPy_adm_integrand);
    DECLARE_FIELD(yPx_adm_integrand);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x = patch->node[ijk]->x[0];
      double y = patch->node[ijk]->x[1]-y_CM;
      xPy_adm_integrand->v[ijk] = x*Py->v[ijk];
      yPx_adm_integrand->v[ijk] = y*Px->v[ijk];
    }
    
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->g00 = obs->ADM[p]->g00;
    I->g01 = obs->ADM[p]->g01;
    I->g02 = obs->ADM[p]->g02;
    I->g11 = obs->ADM[p]->g11;
    I->g12 = obs->ADM[p]->g12;
    I->g22 = obs->ADM[p]->g22;
    I->Spectral->Z_surface = obs->ADM[p]->Z_surface;
    I->Spectral->K         = obs->ADM[p]->K;
    
    I->Spectral->f = xPy_adm_integrand;
    plan_integration(I);
    Jz += execute_integration(I);
    
    I->Spectral->f = yPx_adm_integrand;
    plan_integration(I);
    Jz -= execute_integration(I);
    
    free_integration(I);
    REMOVE_FIELD(xPy_adm_integrand);
    REMOVE_FIELD(yPx_adm_integrand);
  }
  
  /* volume integration */
  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = obs->ADM[p]->patch;
    unsigned nn = patch->nn;
    unsigned ijk;
    
    Field_T *Gy = patch->pool[Ind("ADM_integrand_G_U1")];
    Field_T *Gx = patch->pool[Ind("ADM_integrand_G_U0")];
    ADD_FIELD(xGy_adm_integrand);
    ADD_FIELD(yGx_adm_integrand);
    DECLARE_FIELD(xGy_adm_integrand);
    DECLARE_FIELD(yGx_adm_integrand);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x = patch->node[ijk]->x[0];
      double y = patch->node[ijk]->x[1]-y_CM;
      xGy_adm_integrand->v[ijk] = x*Gy->v[ijk];
      yGx_adm_integrand->v[ijk] = y*Gx->v[ijk];
    }
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->g00 = obs->ADM[p]->g00;
    I->g01 = obs->ADM[p]->g01;
    I->g02 = obs->ADM[p]->g02;
    I->g11 = obs->ADM[p]->g11;
    I->g12 = obs->ADM[p]->g12;
    I->g22 = obs->ADM[p]->g22;
    
    I->Spectral->f = xGy_adm_integrand;
    plan_integration(I);
    Jz -= execute_integration(I);
    
    I->Spectral->f = yGx_adm_integrand;
    plan_integration(I);
    Jz += execute_integration(I);
    
    free_integration(I);
    REMOVE_FIELD(xGy_adm_integrand);
    REMOVE_FIELD(yGx_adm_integrand);
  }
  
  Jz /= (8*M_PI);
  return Jz;
}
