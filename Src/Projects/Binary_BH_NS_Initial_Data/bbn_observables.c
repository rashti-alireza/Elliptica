/*
// Alireza Rashti
// September 2019
*/

/* synopsis:
// =========
//
// * initialize observable *
// Observable_T *obs = init_observable(grid,plan_items_func,free_items_func);
//
// * specifiy which obeservable *
// obs->quantity = "ADMs" # means calculate P and J ADM for the system 
// obs->quantity = "NS_ADMs" # means calculate P and J ADM for single NS 
// obs->quantity = "BH_ADMs" # means calculate P and J ADM for single BH
//
// * plan observable *
// plan_observable(obs);# it finds out the related patches, physical metric etc.
//
// * calculate the observable *
// double Px_ADM = obs->Px_ADM(obs);# x component
// double Py_ADM = obs->Py_ADM(obs);# y component
// double Pz_ADM = obs->Pz_ADM(obs);# z component
// double Jx_ADM = obs->Jx_ADM(obs);# x component of angular momentum
// double Jy_ADM = obs->Jy_ADM(obs);# y component of angular momentum
// double Jz_ADM = obs->Jz_ADM(obs);# z component of angular momentum
//
// *free*
// free_observable(obs);
*/

#include "bbn_observables.h"
#define VOLUME_INTEGRAL 1 /* put it to 1 if you want \int{Gdv} */

/* plan and populate PsJs_ADM_S sturct and obs struct
// for binary and single objects */
void bbn_plan_ADMs_CS(Observable_T *obs)
{
  Grid_T *const grid = obs->grid;
  
  if (!strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
    Error0(NO_OPTION);
  
      
  if (strcmp_i(obs->quantity,"ADMs"))
  {  
    const unsigned N_outermost = (unsigned) Pgeti("Number_of_Outermost_Split");
    Patch_T **patches = 0,*patch = 0;
    char stem[1000];
    struct PsJs_ADM_S **adm = 0;
    unsigned p = 0;
    unsigned n,N,ijk,nn;
    
    if (N_outermost == 0)
      Error0("No outermost patch for integration.\n");
    N = 6*N_outermost/* outermosts */ +
        4/* 4 filling boxes */        +
        10/* 10 sides for surroundings */;
    patches = calloc(N,sizeof(*patches));
    IsNull(patches);  
    
    /* alloc memory for all patches */
    adm = calloc(N,sizeof(*adm));
    IsNull(adm);
    /* this is where we link to obs struct */
    obs->items = adm;
    obs->Nitems = N;
    
    /* first collect all of the patches required */
    p = 0;
    for (n = 0; n < N_outermost; ++n)
    {
      sprintf(stem,"outermost%u_up",n);
      patches[p++]   = GetPatch(stem,grid);
      
      sprintf(stem,"outermost%u_down",n);
      patches[p++] = GetPatch(stem,grid);
      
      sprintf(stem,"outermost%u_back",n);
      patches[p++] = GetPatch(stem,grid);
      
      sprintf(stem,"outermost%u_front",n);
      patches[p++] = GetPatch(stem,grid);
      
      sprintf(stem,"outermost%u_left",n);
      patches[p++] = GetPatch(stem,grid);
      
      sprintf(stem,"outermost%u_right",n);
      patches[p++] = GetPatch(stem,grid);
    }
    /* filling box for vol integrals */
    patches[p++] = GetPatch("filling_box_up",grid);
    patches[p++] = GetPatch("filling_box_down",grid);
    patches[p++] = GetPatch("filling_box_back",grid);
    patches[p++] = GetPatch("filling_box_front",grid);
    
    /* surroundings for surface integrals */
    patches[p++] = GetPatch("left_NS_surrounding_up",grid);
    patches[p++] = GetPatch("left_NS_surrounding_down",grid);
    patches[p++] = GetPatch("left_NS_surrounding_left",grid);
    patches[p++] = GetPatch("left_NS_surrounding_back",grid);
    patches[p++] = GetPatch("left_NS_surrounding_front",grid);
    
    patches[p++] = GetPatch("right_BH_surrounding_up",grid);
    patches[p++] = GetPatch("right_BH_surrounding_down",grid);
    patches[p++] = GetPatch("right_BH_surrounding_right",grid);
    patches[p++] = GetPatch("right_BH_surrounding_back",grid);
    patches[p++] = GetPatch("right_BH_surrounding_front",grid);
    
    assert(p==N);
    
    /* fill ADM struct for each patch */
    for (n = 0; n < N; ++n)
    {
      adm[n] = calloc(1,sizeof(*adm[n]));
      IsNull(adm[n]);
      patch = patches[n];
      nn    = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);
      
      READ_v(_gamma_D2D2)
      READ_v(_gamma_D0D2)
      READ_v(_gamma_D0D0)
      READ_v(_gamma_D0D1)
      READ_v(_gamma_D1D2)
      READ_v(_gamma_D1D1)
      READ_v(psi);
      
      adm[n]->patch = patch;
      /* populate metric components */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*_gamma_D0D0[ijk];
        g01[ijk] = psi4*_gamma_D0D1[ijk];
        g02[ijk] = psi4*_gamma_D0D2[ijk];
        g11[ijk] = psi4*_gamma_D1D1[ijk];
        g12[ijk] = psi4*_gamma_D1D2[ijk];
        g22[ijk] = psi4*_gamma_D2D2[ijk];
      }
      adm[n]->g00 = g00;
      adm[n]->g01 = g01;
      adm[n]->g02 = g02;
      adm[n]->g11 = g11;
      adm[n]->g12 = g12;
      adm[n]->g22 = g22;
      
      /* surface integrals params */
      if (regex_search(".+(left|right)_(NS|BH)_surrounding.+",patch->name))
      {
        adm[n]->surface_integration_flg = 1;
        adm[n]->Z_surface = 1;
        adm[n]->K = patch->n[2]-1;
        populate_normal_surrounding(adm[n],_c_);
      }
    }
    obs->Px_ADM = ADM_momentum_x_BBN_CS;
    obs->Py_ADM = ADM_momentum_y_BBN_CS;
    obs->Pz_ADM = ADM_momentum_z_BBN_CS;
    obs->Jx_ADM = ADM_angular_momentum_x_BBN_CS;
    obs->Jy_ADM = ADM_angular_momentum_y_BBN_CS;
    obs->Jz_ADM = ADM_angular_momentum_z_BBN_CS;
    bbn_populate_ADM_integrand_PdS_GdV(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"NS_ADMs"))
  {  
    Patch_T **patches = 0,*patch = 0;
    struct PsJs_ADM_S **adm = 0;
    unsigned p = 0;
    unsigned n,N,ijk,nn;
    
    N = 6/* 6 sides for surroundings */;
        
    patches = calloc(N,sizeof(*patches));
    IsNull(patches);  
    
    /* alloc memory for all patches */
    adm = calloc(N,sizeof(*adm));
    IsNull(adm);
    /* this is where we link to obs struct */
    obs->items = adm;
    obs->Nitems = N;
    
    /* first collect all of the patches required */
    p = 0;
    /* surroundings for surface integrals */
    patches[p++] = GetPatch("left_NS_surrounding_up",grid);
    patches[p++] = GetPatch("left_NS_surrounding_down",grid);
    patches[p++] = GetPatch("left_NS_surrounding_left",grid);
    patches[p++] = GetPatch("left_NS_surrounding_right",grid);
    patches[p++] = GetPatch("left_NS_surrounding_back",grid);
    patches[p++] = GetPatch("left_NS_surrounding_front",grid);
    
    assert(p==N);
    
    /* fill ADM struct for each patch */
    for (n = 0; n < N; ++n)
    {
      adm[n] = calloc(1,sizeof(*adm[n]));
      IsNull(adm[n]);
      patch = patches[n];
      nn    = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);
      
      READ_v(_gamma_D2D2)
      READ_v(_gamma_D0D2)
      READ_v(_gamma_D0D0)
      READ_v(_gamma_D0D1)
      READ_v(_gamma_D1D2)
      READ_v(_gamma_D1D1)
      READ_v(psi);
      
      adm[n]->patch = patch;
      /* populate metric components */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*_gamma_D0D0[ijk];
        g01[ijk] = psi4*_gamma_D0D1[ijk];
        g02[ijk] = psi4*_gamma_D0D2[ijk];
        g11[ijk] = psi4*_gamma_D1D1[ijk];
        g12[ijk] = psi4*_gamma_D1D2[ijk];
        g22[ijk] = psi4*_gamma_D2D2[ijk];
      }
      adm[n]->g00 = g00;
      adm[n]->g01 = g01;
      adm[n]->g02 = g02;
      adm[n]->g11 = g11;
      adm[n]->g12 = g12;
      adm[n]->g22 = g22;
      
      /* surface integral */
      adm[n]->surface_integration_flg = 1;
      adm[n]->Z_surface = 1;
      adm[n]->K = 0;
      populate_normal_surrounding(adm[n],_c_);
    }
    obs->Px_ADM = ADM_momentum_x_BBN_CS;
    obs->Py_ADM = ADM_momentum_y_BBN_CS;
    obs->Pz_ADM = ADM_momentum_z_BBN_CS;
    obs->Jx_ADM = ADM_angular_momentum_x_BBN_CS;
    obs->Jy_ADM = ADM_angular_momentum_y_BBN_CS;
    obs->Jz_ADM = ADM_angular_momentum_z_BBN_CS;
    bbn_populate_ADM_integrand_PdS_GdV(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"BH_ADMs"))
  {  
    Patch_T **patches = 0,*patch = 0;
    struct PsJs_ADM_S **adm = 0;
    unsigned p = 0;
    unsigned n,N,ijk,nn;
    
    N = 6/* 6 sides for surroundings */;
        
    patches = calloc(N,sizeof(*patches));
    IsNull(patches);  
    
    /* alloc memory for all patches */
    adm = calloc(N,sizeof(*adm));
    IsNull(adm);
    /* this is where we link to obs struct */
    obs->items = adm;
    obs->Nitems = N;
    
    /* first collect all of the patches required */
    p = 0;
    /* surroundings for surface integrals */
    patches[p++] = GetPatch("right_BH_surrounding_up",grid);
    patches[p++] = GetPatch("right_BH_surrounding_down",grid);
    patches[p++] = GetPatch("right_BH_surrounding_left",grid);
    patches[p++] = GetPatch("right_BH_surrounding_right",grid);
    patches[p++] = GetPatch("right_BH_surrounding_back",grid);
    patches[p++] = GetPatch("right_BH_surrounding_front",grid);
    
    assert(p==N);
    
    /* fill ADM struct for each patch */
    for (n = 0; n < N; ++n)
    {
      adm[n] = calloc(1,sizeof(*adm[n]));
      IsNull(adm[n]);
      patch = patches[n];
      nn    = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);
      
      READ_v(_gamma_D2D2)
      READ_v(_gamma_D0D2)
      READ_v(_gamma_D0D0)
      READ_v(_gamma_D0D1)
      READ_v(_gamma_D1D2)
      READ_v(_gamma_D1D1)
      READ_v(psi);
      
      adm[n]->patch = patch;
      /* populate metric components */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*_gamma_D0D0[ijk];
        g01[ijk] = psi4*_gamma_D0D1[ijk];
        g02[ijk] = psi4*_gamma_D0D2[ijk];
        g11[ijk] = psi4*_gamma_D1D1[ijk];
        g12[ijk] = psi4*_gamma_D1D2[ijk];
        g22[ijk] = psi4*_gamma_D2D2[ijk];
      }
      adm[n]->g00 = g00;
      adm[n]->g01 = g01;
      adm[n]->g02 = g02;
      adm[n]->g11 = g11;
      adm[n]->g12 = g12;
      adm[n]->g22 = g22;
      
      /* surface integral */
      adm[n]->surface_integration_flg = 1;
      adm[n]->Z_surface = 1;
      adm[n]->K = 0;
      populate_normal_surrounding(adm[n],_c_);
    }
    obs->Px_ADM = ADM_momentum_x_BBN_CS;
    obs->Py_ADM = ADM_momentum_y_BBN_CS;
    obs->Pz_ADM = ADM_momentum_z_BBN_CS;
    obs->Jx_ADM = ADM_angular_momentum_x_BBN_CS;
    obs->Jy_ADM = ADM_angular_momentum_y_BBN_CS;
    obs->Jz_ADM = ADM_angular_momentum_z_BBN_CS;
    bbn_populate_ADM_integrand_PdS_GdV(obs);
    free(patches);
  }
  else
    Error1("This is no plan for %s.\n",obs->quantity);
  
}

/* populating normal outward vector for surrounding according to the given dir */
static void populate_normal_surrounding(struct PsJs_ADM_S *const adm,const Dd_T dir)
{
  Patch_T *const patch = adm->patch;
  const unsigned nn = patch->nn;
  double *n_U0 = alloc_double(nn);
  double *n_U1 = alloc_double(nn);
  double *n_U2 = alloc_double(nn);
  unsigned ijk;
  
  READ_v(_gamma_D2D2)
  READ_v(_gamma_D0D2)
  READ_v(_gamma_D0D0)
  READ_v(_gamma_D0D1)
  READ_v(_gamma_D1D2)
  READ_v(_gamma_D1D1)
  READ_v(psi);
    
  for (ijk = 0; ijk < nn; ++ijk)
  {
    n_U0[ijk] = dq2_dq1(patch,dir,_x_,ijk);
    n_U1[ijk] = dq2_dq1(patch,dir,_y_,ijk);
    n_U2[ijk] = dq2_dq1(patch,dir,_z_,ijk);
    
    /* normalization */
    double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
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
  adm->n_U0 = n_U0;
  adm->n_U1 = n_U1;
  adm->n_U2 = n_U2;
}

/* free stuct Observable_T */
void bbn_free_ADMs_CS(Observable_T *obs)
{
  if (!obs)
    return;
    
  if (!strcmp_i(obs->grid->kind,"BBN_CubedSpherical_grid"))
    Error0(NO_OPTION);
  
  struct PsJs_ADM_S **adm = obs->items;
  unsigned i;
    
  for (i = 0; i < obs->Nitems; ++i)
  {
    Patch_T *patch = adm[i]->patch;
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
    _free(adm[i]->g00);
    _free(adm[i]->g01);
    _free(adm[i]->g02);
    _free(adm[i]->g11);
    _free(adm[i]->g12);
    _free(adm[i]->g22);
    _free(adm[i]->n_U0);
    _free(adm[i]->n_U1);
    _free(adm[i]->n_U2);
    
    free(adm[i]);
  }
  _free(adm);
}

/* calculating ADM momentum in x component */
static double ADM_momentum_x_BBN_CS(Observable_T *const obs)
{
  double Px = 0;
  struct PsJs_ADM_S **const adm = obs->items;
  const unsigned N = obs->Nitems;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  //populate_ADM_momentums_integrand_PdS_GdV(obs);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    if (!adm[p]->surface_integration_flg)
      continue;
      
    Patch_T *patch = adm[p]->patch;
    
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_P_U0")];
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
    I->Spectral->X_surface = adm[p]->X_surface;
    I->Spectral->Y_surface = adm[p]->Y_surface;
    I->Spectral->Z_surface = adm[p]->Z_surface;
    I->Spectral->I         = adm[p]->I;
    I->Spectral->J         = adm[p]->J;
    I->Spectral->K         = adm[p]->K;
    
    plan_integration(I);
    Px += execute_integration(I);
    free_integration(I);
  }
  
  /* volume integration */
  if (VOLUME_INTEGRAL)
  for(p = 0; p < N; ++p)
  {
    if (adm[p]->surface_integration_flg)
      continue;
    
    Patch_T *patch = adm[p]->patch;
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_G_U0")];
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
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
  struct PsJs_ADM_S **const adm = obs->items;
  const unsigned N = obs->Nitems;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  //populate_ADM_momentums_integrand_PdS_GdV(obs);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    if (!adm[p]->surface_integration_flg)
      continue;
      
    Patch_T *patch = adm[p]->patch;
    
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_P_U1")];
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
    I->Spectral->X_surface = adm[p]->X_surface;
    I->Spectral->Y_surface = adm[p]->Y_surface;
    I->Spectral->Z_surface = adm[p]->Z_surface;
    I->Spectral->I         = adm[p]->I;
    I->Spectral->J         = adm[p]->J;
    I->Spectral->K         = adm[p]->K;
    
    plan_integration(I);
    Py += execute_integration(I);
    free_integration(I);
  }
  
  /* volume integration */
  if (VOLUME_INTEGRAL)
  for(p = 0; p < N; ++p)
  {
    if (adm[p]->surface_integration_flg)
      continue;
      
    Patch_T *patch = adm[p]->patch;
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_G_U1")];
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
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
  struct PsJs_ADM_S **const adm = obs->items;
  const unsigned N = obs->Nitems;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  //populate_ADM_momentums_integrand_PdS_GdV(obs);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    if (!adm[p]->surface_integration_flg)
      continue;
      
    Patch_T *patch = adm[p]->patch;
    
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_P_U2")];
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
    I->Spectral->X_surface = adm[p]->X_surface;
    I->Spectral->Y_surface = adm[p]->Y_surface;
    I->Spectral->Z_surface = adm[p]->Z_surface;
    I->Spectral->I         = adm[p]->I;
    I->Spectral->J         = adm[p]->J;
    I->Spectral->K         = adm[p]->K;
    
    plan_integration(I);
    Pz += execute_integration(I);
    free_integration(I);
  }
  
  /* volume integration */
  if (VOLUME_INTEGRAL)
  for(p = 0; p < N; ++p)
  {
    if (adm[p]->surface_integration_flg)
      continue;

    Patch_T *patch = adm[p]->patch;
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->Spectral->f = patch->pool[Ind("ADM_integrand_G_U2")];
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
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
  struct PsJs_ADM_S **const adm = obs->items;
  const unsigned N = obs->Nitems;
  double y_CM = 0;
  double x_CM = 0;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  //populate_ADM_momentums_integrand_PdS_GdV(obs);
  
  if (get_parameter("y_CM"))
    y_CM = Pgetd("y_CM");
    
  if (get_parameter("x_CM"))
    x_CM = Pgetd("x_CM");
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    if (!adm[p]->surface_integration_flg)
      continue;
    
    Patch_T *patch = adm[p]->patch;
    unsigned nn = patch->nn;
    unsigned ijk;
    
    Field_T *Py = patch->pool[Ind("ADM_integrand_P_U1")];
    Field_T *Px = patch->pool[Ind("ADM_integrand_P_U0")];
    ADD_AND_ALLOC_FIELD(xPy_adm_integrand);
    ADD_AND_ALLOC_FIELD(yPx_adm_integrand);
    DECLARE_FIELD(xPy_adm_integrand);
    DECLARE_FIELD(yPx_adm_integrand);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x = patch->node[ijk]->x[0]-x_CM;
      double y = patch->node[ijk]->x[1]-y_CM;
      xPy_adm_integrand->v[ijk] = x*Py->v[ijk];
      yPx_adm_integrand->v[ijk] = y*Px->v[ijk];
    }
    
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
    I->Spectral->X_surface = adm[p]->X_surface;
    I->Spectral->Y_surface = adm[p]->Y_surface;
    I->Spectral->Z_surface = adm[p]->Z_surface;
    I->Spectral->I         = adm[p]->I;
    I->Spectral->J         = adm[p]->J;
    I->Spectral->K         = adm[p]->K;
        
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
  if (VOLUME_INTEGRAL)
  for(p = 0; p < N; ++p)
  {
    if (adm[p]->surface_integration_flg)
      continue;

    Patch_T *patch = adm[p]->patch;
    unsigned nn = patch->nn;
    unsigned ijk;
    
    Field_T *Gy = patch->pool[Ind("ADM_integrand_G_U1")];
    Field_T *Gx = patch->pool[Ind("ADM_integrand_G_U0")];
    ADD_AND_ALLOC_FIELD(xGy_adm_integrand);
    ADD_AND_ALLOC_FIELD(yGx_adm_integrand);
    DECLARE_FIELD(xGy_adm_integrand);
    DECLARE_FIELD(yGx_adm_integrand);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x = patch->node[ijk]->x[0]-x_CM;
      double y = patch->node[ijk]->x[1]-y_CM;
      xGy_adm_integrand->v[ijk] = x*Gy->v[ijk];
      yGx_adm_integrand->v[ijk] = y*Gx->v[ijk];
    }
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
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

/* calculating ADM angular momentum in x component */
static double ADM_angular_momentum_x_BBN_CS(Observable_T *const obs)
{
  double Jx = 0;
  struct PsJs_ADM_S **const adm = obs->items;
  const unsigned N = obs->Nitems;
  double y_CM = 0;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  //populate_ADM_momentums_integrand_PdS_GdV(obs);
  
  if (get_parameter("y_CM"))
    y_CM = Pgetd("y_CM");
    
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
   if (!adm[p]->surface_integration_flg)
      continue;
   
    Patch_T *patch = adm[p]->patch;
    unsigned nn = patch->nn;
    unsigned ijk;
    
    Field_T *Py = patch->pool[Ind("ADM_integrand_P_U1")];
    Field_T *Pz = patch->pool[Ind("ADM_integrand_P_U2")];
    
    ADD_AND_ALLOC_FIELD(zPy_adm_integrand);
    ADD_AND_ALLOC_FIELD(yPz_adm_integrand);
    DECLARE_FIELD(zPy_adm_integrand);
    DECLARE_FIELD(yPz_adm_integrand);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double z = patch->node[ijk]->x[2];
      double y = patch->node[ijk]->x[1]-y_CM;
      zPy_adm_integrand->v[ijk] = z*Py->v[ijk];
      yPz_adm_integrand->v[ijk] = y*Pz->v[ijk];
    }
    
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    I->Spectral->X_surface = adm[p]->X_surface;
    I->Spectral->Y_surface = adm[p]->Y_surface;
    I->Spectral->Z_surface = adm[p]->Z_surface;
    I->Spectral->I         = adm[p]->I;
    I->Spectral->J         = adm[p]->J;
    I->Spectral->K         = adm[p]->K;
    
    
    I->Spectral->f = zPy_adm_integrand;
    plan_integration(I);
    Jx -= execute_integration(I);
    
    I->Spectral->f = yPz_adm_integrand;
    plan_integration(I);
    Jx += execute_integration(I);
    
    free_integration(I);
    REMOVE_FIELD(zPy_adm_integrand);
    REMOVE_FIELD(yPz_adm_integrand);
  }
  
  /* volume integration */
  if (VOLUME_INTEGRAL)
  for(p = 0; p < N; ++p)
  {
    if (adm[p]->surface_integration_flg)
      continue;
   
    Patch_T *patch = adm[p]->patch;
    unsigned nn = patch->nn;
    unsigned ijk;
    
    Field_T *Gy = patch->pool[Ind("ADM_integrand_G_U1")];
    Field_T *Gz = patch->pool[Ind("ADM_integrand_G_U2")];
    ADD_AND_ALLOC_FIELD(zGy_adm_integrand);
    ADD_AND_ALLOC_FIELD(yGz_adm_integrand);
    DECLARE_FIELD(zGy_adm_integrand);
    DECLARE_FIELD(yGz_adm_integrand);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double z = patch->node[ijk]->x[2];
      double y = patch->node[ijk]->x[1]-y_CM;
      zGy_adm_integrand->v[ijk] = z*Gy->v[ijk];
      yGz_adm_integrand->v[ijk] = y*Gz->v[ijk];
    }
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
    I->Spectral->f = zGy_adm_integrand;
    plan_integration(I);
    Jx += execute_integration(I);
    
    I->Spectral->f = yGz_adm_integrand;
    plan_integration(I);
    Jx -= execute_integration(I);
    
    free_integration(I);
    REMOVE_FIELD(zGy_adm_integrand);
    REMOVE_FIELD(yGz_adm_integrand);
  }
  
  Jx /= (8*M_PI);
  return Jx;
}

/* calculating ADM angular momentum in y component */
static double ADM_angular_momentum_y_BBN_CS(Observable_T *const obs)
{
  double Jy = 0;
  struct PsJs_ADM_S **const adm = obs->items;
  const unsigned N = obs->Nitems;
  double x_CM = 0;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  //populate_ADM_momentums_integrand_PdS_GdV(obs);
  
  if (get_parameter("x_CM"))
    x_CM = Pgetd("x_CM");
    
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    if (!adm[p]->surface_integration_flg)
      continue;
    
    Patch_T *patch = adm[p]->patch;
    unsigned nn = patch->nn;
    unsigned ijk;
    
    Field_T *Px = patch->pool[Ind("ADM_integrand_P_U0")];
    Field_T *Pz = patch->pool[Ind("ADM_integrand_P_U2")];
    ADD_AND_ALLOC_FIELD(zPx_adm_integrand);
    ADD_AND_ALLOC_FIELD(xPz_adm_integrand);
    DECLARE_FIELD(zPx_adm_integrand);
    DECLARE_FIELD(xPz_adm_integrand);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double z = patch->node[ijk]->x[2];
      double x = patch->node[ijk]->x[0]-x_CM;
      zPx_adm_integrand->v[ijk] = z*Px->v[ijk];
      xPz_adm_integrand->v[ijk] = x*Pz->v[ijk];
    }
    
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
    I->Spectral->X_surface = adm[p]->X_surface;
    I->Spectral->Y_surface = adm[p]->Y_surface;
    I->Spectral->Z_surface = adm[p]->Z_surface;
    I->Spectral->I         = adm[p]->I;
    I->Spectral->J         = adm[p]->J;
    I->Spectral->K         = adm[p]->K;
        
    I->Spectral->f = zPx_adm_integrand;
    plan_integration(I);
    Jy += execute_integration(I);
    
    I->Spectral->f = xPz_adm_integrand;
    plan_integration(I);
    Jy -= execute_integration(I);
    
    free_integration(I);
    REMOVE_FIELD(zPx_adm_integrand);
    REMOVE_FIELD(xPz_adm_integrand);
  }
  
  /* volume integration */
  if (VOLUME_INTEGRAL)
  for(p = 0; p < N; ++p)
  {
    if (adm[p]->surface_integration_flg)
      continue;
   
    Patch_T *patch = adm[p]->patch;
    unsigned nn = patch->nn;
    unsigned ijk;
    
    Field_T *Gx = patch->pool[Ind("ADM_integrand_G_U0")];
    Field_T *Gz = patch->pool[Ind("ADM_integrand_G_U2")];
    ADD_AND_ALLOC_FIELD(zGx_adm_integrand);
    ADD_AND_ALLOC_FIELD(xGz_adm_integrand);
    DECLARE_FIELD(zGx_adm_integrand);
    DECLARE_FIELD(xGz_adm_integrand);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double z = patch->node[ijk]->x[2];
      double x = patch->node[ijk]->x[0]-x_CM;
      zGx_adm_integrand->v[ijk] = z*Gx->v[ijk];
      xGz_adm_integrand->v[ijk] = x*Gz->v[ijk];
    }
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
    I->Spectral->f = zGx_adm_integrand;
    plan_integration(I);
    Jy -= execute_integration(I);
    
    I->Spectral->f = xGz_adm_integrand;
    plan_integration(I);
    Jy += execute_integration(I);
    
    free_integration(I);
    REMOVE_FIELD(zGx_adm_integrand);
    REMOVE_FIELD(xGz_adm_integrand);
  }
  
  Jy /= (8*M_PI);
  return Jy;
}
