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
// obs->quantity = "ADM(P,J)|BBN"  #=> compute P and J ADM for the system 
// obs->quantity = "ADM(P,J)|NS"   #=> compute P and J ADM for single NS 
// obs->quantity = "ADM(P,J)|BH"   #=> compute P and J ADM for single BH
// obs->quantity = "Kommar(M)|BBN" #=> compute Kommar mass for the system 
// obs->quantity = "Kommar(M)|NS"  #=> compute kommar mass for NS 
// obs->quantity = "Kommar(M)|BH"  #=> compute Kommar mass for BH
// obs->quantity = "ADM(M)|BBN"    #=> compute ADM mass for the system 
// obs->quantity = "ADM(M)|NS"     #=> compute ADM mass for NS 
// obs->quantity = "ADM(M)|BH"     #=> compute ADM mass for BH
//
// * plan observable *
// plan_observable(obs);# it finds out the related patches, physical metric etc.
//
// * calculate the observable example:*
// double Px_ADM = obs->Px(obs);# x component
// double Py_ADM = obs->Py(obs);# y component
// double Pz_ADM = obs->Pz(obs);# z component
// double Jx_ADM = obs->Jx(obs);# x component of angular momentum
// double Jy_ADM = obs->Jy(obs);# y component of angular momentum
// double Jz_ADM = obs->Jz(obs);# z component of angular momentum
// double M_ADM  = obs->M(obs) ;# a specifed mass for example ADM mass
// *free*
// free_observable(obs);
*/

#include "obs_header.h"
#define VOLUME_INTEGRAL 1 /* put it to 1 if you want \int{Gdv} */

/* plan and populate items_S sturct and obs struct
// for binary and single objects.
// algorithm: 
// ==========
//
// 1. collect all of the necessary patches.
// 2. populate the required metric for the integrations.
// 3. set flags for surface or volume integral.
// 4. populate normal vectors for surface integrals.
// 5. populate the integrands
// 6. assign the pertinent functions for the calculation.
// */
void obs_plan(Observable_T *obs)
{
  Grid_T *const grid = obs->grid;
  
  if (!strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
    Error0(NO_OPTION);
      
  if (strcmp_i(obs->quantity,"ADM(P,J)|BBN"))
  {  
    const unsigned N_outermost = (unsigned) Pgeti("Number_of_Outermost_Split");
    Patch_T **patches = 0,*patch = 0;
    char stem[1000];
    struct items_S **adm = 0;
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
      
      READ_v(gConf_D2D2)
      READ_v(gConf_D0D2)
      READ_v(gConf_D0D0)
      READ_v(gConf_D0D1)
      READ_v(gConf_D1D2)
      READ_v(gConf_D1D1)
      READ_v(psi);
      
      adm[n]->patch = patch;
      /* populate metric components */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*gConf_D0D0[ijk];
        g01[ijk] = psi4*gConf_D0D1[ijk];
        g02[ijk] = psi4*gConf_D0D2[ijk];
        g11[ijk] = psi4*gConf_D1D1[ijk];
        g12[ijk] = psi4*gConf_D1D2[ijk];
        g22[ijk] = psi4*gConf_D2D2[ijk];
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
        n_physical_metric_surrounding(adm[n],_c_);
      }
    }
    obs->Px = ADM_momentum_x_BHNS_CS;
    obs->Py = ADM_momentum_y_BHNS_CS;
    obs->Pz = ADM_momentum_z_BHNS_CS;
    obs->Jx = ADM_angular_momentum_x_BHNS_CS;
    obs->Jy = ADM_angular_momentum_y_BHNS_CS;
    obs->Jz = ADM_angular_momentum_z_BHNS_CS;
    obs_populate_ADM_integrand_PdS_GdV_binary(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"ADM(P,J)|NS"))
  {  
    Patch_T **patches = 0,*patch = 0;
    struct items_S **adm = 0;
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
      
      READ_v(gConf_D2D2)
      READ_v(gConf_D0D2)
      READ_v(gConf_D0D0)
      READ_v(gConf_D0D1)
      READ_v(gConf_D1D2)
      READ_v(gConf_D1D1)
      READ_v(psi);
      
      adm[n]->patch = patch;
      /* populate metric components */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*gConf_D0D0[ijk];
        g01[ijk] = psi4*gConf_D0D1[ijk];
        g02[ijk] = psi4*gConf_D0D2[ijk];
        g11[ijk] = psi4*gConf_D1D1[ijk];
        g12[ijk] = psi4*gConf_D1D2[ijk];
        g22[ijk] = psi4*gConf_D2D2[ijk];
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
      n_physical_metric_surrounding(adm[n],_c_);
    }
    obs->Px = ADM_momentum_x_BHNS_CS;
    obs->Py = ADM_momentum_y_BHNS_CS;
    obs->Pz = ADM_momentum_z_BHNS_CS;
    obs->Jx = ADM_angular_momentum_x_BHNS_CS;
    obs->Jy = ADM_angular_momentum_y_BHNS_CS;
    obs->Jz = ADM_angular_momentum_z_BHNS_CS;
    obs_populate_ADM_integrand_PdS_GdV_single(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"ADM(P,J)|BH"))
  {  
    Patch_T **patches = 0,*patch = 0;
    struct items_S **adm = 0;
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
      
      READ_v(gConf_D2D2)
      READ_v(gConf_D0D2)
      READ_v(gConf_D0D0)
      READ_v(gConf_D0D1)
      READ_v(gConf_D1D2)
      READ_v(gConf_D1D1)
      READ_v(psi);
      
      adm[n]->patch = patch;
      /* populate metric components */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*gConf_D0D0[ijk];
        g01[ijk] = psi4*gConf_D0D1[ijk];
        g02[ijk] = psi4*gConf_D0D2[ijk];
        g11[ijk] = psi4*gConf_D1D1[ijk];
        g12[ijk] = psi4*gConf_D1D2[ijk];
        g22[ijk] = psi4*gConf_D2D2[ijk];
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
      n_physical_metric_surrounding(adm[n],_c_);
    }
    obs->Px = ADM_momentum_x_BHNS_CS;
    obs->Py = ADM_momentum_y_BHNS_CS;
    obs->Pz = ADM_momentum_z_BHNS_CS;
    obs->Jx = ADM_angular_momentum_x_BHNS_CS;
    obs->Jy = ADM_angular_momentum_y_BHNS_CS;
    obs->Jz = ADM_angular_momentum_z_BHNS_CS;
    obs_populate_ADM_integrand_PdS_GdV_single(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"Kommar(M)|BH"))
  {  
    Patch_T **patches = 0,*patch = 0;
    struct items_S **kommar = 0;
    unsigned p = 0;
    unsigned n,N,ijk,nn;
    
    N = 6/* 6 sides for surroundings */;
        
    patches = calloc(N,sizeof(*patches));
    IsNull(patches);  
    
    /* alloc memory for all patches */
    kommar = calloc(N,sizeof(*kommar));
    IsNull(kommar);
    /* this is where we link to obs struct */
    obs->items = kommar;
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
    
    /* fill kommar struct for each patch */
    for (n = 0; n < N; ++n)
    {
      kommar[n] = calloc(1,sizeof(*kommar[n]));
      IsNull(kommar[n]);
      patch = patches[n];
      nn    = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);
      
      READ_v(gConf_D2D2)
      READ_v(gConf_D0D2)
      READ_v(gConf_D0D0)
      READ_v(gConf_D0D1)
      READ_v(gConf_D1D2)
      READ_v(gConf_D1D1)
      READ_v(psi);
      
      kommar[n]->patch = patch;
      /* populate metric components */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*gConf_D0D0[ijk];
        g01[ijk] = psi4*gConf_D0D1[ijk];
        g02[ijk] = psi4*gConf_D0D2[ijk];
        g11[ijk] = psi4*gConf_D1D1[ijk];
        g12[ijk] = psi4*gConf_D1D2[ijk];
        g22[ijk] = psi4*gConf_D2D2[ijk];
      }
      kommar[n]->g00 = g00;
      kommar[n]->g01 = g01;
      kommar[n]->g02 = g02;
      kommar[n]->g11 = g11;
      kommar[n]->g12 = g12;
      kommar[n]->g22 = g22;
      
      /* surface integral */
      kommar[n]->surface_integration_flg = 1;
      kommar[n]->Z_surface = 1;
      kommar[n]->K = 0;
      n_physical_metric_surrounding(kommar[n],_c_);
    }
    obs->M = obs_Kommar_mass;
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"Kommar(M)|NS"))
  {  
    Patch_T **patches = 0,*patch = 0;
    struct items_S **kommar = 0;
    unsigned p = 0;
    unsigned n,N,ijk,nn;
    
    N = 6/* 6 sides for surroundings */;
        
    patches = calloc(N,sizeof(*patches));
    IsNull(patches);  
    
    /* alloc memory for all patches */
    kommar = calloc(N,sizeof(*kommar));
    IsNull(kommar);
    /* this is where we link to obs struct */
    obs->items = kommar;
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
    
    /* fill kommar struct for each patch */
    for (n = 0; n < N; ++n)
    {
      kommar[n] = calloc(1,sizeof(*kommar[n]));
      IsNull(kommar[n]);
      patch = patches[n];
      nn    = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);
      
      READ_v(gConf_D2D2)
      READ_v(gConf_D0D2)
      READ_v(gConf_D0D0)
      READ_v(gConf_D0D1)
      READ_v(gConf_D1D2)
      READ_v(gConf_D1D1)
      READ_v(psi);
      
      kommar[n]->patch = patch;
      /* populate metric components */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*gConf_D0D0[ijk];
        g01[ijk] = psi4*gConf_D0D1[ijk];
        g02[ijk] = psi4*gConf_D0D2[ijk];
        g11[ijk] = psi4*gConf_D1D1[ijk];
        g12[ijk] = psi4*gConf_D1D2[ijk];
        g22[ijk] = psi4*gConf_D2D2[ijk];
      }
      kommar[n]->g00 = g00;
      kommar[n]->g01 = g01;
      kommar[n]->g02 = g02;
      kommar[n]->g11 = g11;
      kommar[n]->g12 = g12;
      kommar[n]->g22 = g22;
      
      /* surface integral */
      kommar[n]->surface_integration_flg = 1;
      kommar[n]->Z_surface = 1;
      kommar[n]->K = 0;
      n_physical_metric_surrounding(kommar[n],_c_);
    }
    obs->M = obs_Kommar_mass;
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"Kommar(M)|BBN"))
  {  
    Patch_T **patches = 0,*patch = 0;
    struct items_S **kommar = 0;
    unsigned p = 0;
    unsigned n,N,ijk,nn;
    
    N = 10/* 10 sides for surroundings */;
        
    patches = calloc(N,sizeof(*patches));
    IsNull(patches);  
    
    /* alloc memory for all patches */
    kommar = calloc(N,sizeof(*kommar));
    IsNull(kommar);
    /* this is where we link to obs struct */
    obs->items = kommar;
    obs->Nitems = N;
    
    /* first collect all of the patches required */
    p = 0;
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
    
    /* fill kommar struct for each patch */
    for (n = 0; n < N; ++n)
    {
      kommar[n] = calloc(1,sizeof(*kommar[n]));
      IsNull(kommar[n]);
      patch = patches[n];
      nn    = patch->nn;
      
      double *g00 = alloc_double(nn);
      double *g01 = alloc_double(nn);
      double *g02 = alloc_double(nn);
      double *g11 = alloc_double(nn);
      double *g12 = alloc_double(nn);
      double *g22 = alloc_double(nn);
      
      READ_v(gConf_D2D2)
      READ_v(gConf_D0D2)
      READ_v(gConf_D0D0)
      READ_v(gConf_D0D1)
      READ_v(gConf_D1D2)
      READ_v(gConf_D1D1)
      READ_v(psi);
      
      kommar[n]->patch = patch;
      /* populate metric components */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*gConf_D0D0[ijk];
        g01[ijk] = psi4*gConf_D0D1[ijk];
        g02[ijk] = psi4*gConf_D0D2[ijk];
        g11[ijk] = psi4*gConf_D1D1[ijk];
        g12[ijk] = psi4*gConf_D1D2[ijk];
        g22[ijk] = psi4*gConf_D2D2[ijk];
      }
      kommar[n]->g00 = g00;
      kommar[n]->g01 = g01;
      kommar[n]->g02 = g02;
      kommar[n]->g11 = g11;
      kommar[n]->g12 = g12;
      kommar[n]->g22 = g22;
      
      /* surface integral */
      kommar[n]->surface_integration_flg = 1;
      kommar[n]->Z_surface = 1;
      kommar[n]->K = patch->n[2]-1;
      n_physical_metric_surrounding(kommar[n],_c_);
    }
    obs->M = obs_Kommar_mass;
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"ADM(M)|BBN"))
  {  
    const unsigned N_outermost = (unsigned) Pgeti("Number_of_Outermost_Split");
    Patch_T **patches = 0,*patch = 0;
    char stem[1000];
    struct items_S **adm = 0;
    unsigned p = 0,surface_index = 0;
    unsigned n,N,ijk,nn;
    
    if (N_outermost == 0)
      Error0("No outermost patch for integration.\n");
    N = 6*N_outermost/* outermosts */ +
        4/* 4 filling boxes */        +
        19/* box+NS+surroundings for vol integral */ +
        6/* 6 sides for surroundings for area integral */;
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
    
    patches[p++] = GetPatch("left_NS_up",grid);
    patches[p++] = GetPatch("left_NS_down",grid);
    patches[p++] = GetPatch("left_NS_left",grid);
    patches[p++] = GetPatch("left_NS_right",grid);
    patches[p++] = GetPatch("left_NS_back",grid);
    patches[p++] = GetPatch("left_NS_front",grid);
    patches[p++] = GetPatch("left_central_box",grid);
      
    patches[p++] = GetPatch("left_NS_surrounding_up",grid);
    patches[p++] = GetPatch("left_NS_surrounding_down",grid);
    patches[p++] = GetPatch("left_NS_surrounding_left",grid);
    patches[p++] = GetPatch("left_NS_surrounding_right",grid);
    patches[p++] = GetPatch("left_NS_surrounding_back",grid);
    patches[p++] = GetPatch("left_NS_surrounding_front",grid);
    
    patches[p++] = GetPatch("right_BH_surrounding_up",grid);
    patches[p++] = GetPatch("right_BH_surrounding_down",grid);
    patches[p++] = GetPatch("right_BH_surrounding_left",grid);
    patches[p++] = GetPatch("right_BH_surrounding_right",grid);
    patches[p++] = GetPatch("right_BH_surrounding_back",grid);
    patches[p++] = GetPatch("right_BH_surrounding_front",grid);
    
    /* surroundings for surface integrals. NOTE: the order matters */
    surface_index = p;
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
      
      READ_v(gConf_D2D2)
      READ_v(gConf_D0D2)
      READ_v(gConf_D0D0)
      READ_v(gConf_D0D1)
      READ_v(gConf_D1D2)
      READ_v(gConf_D1D1)
      
      adm[n]->patch = patch;
      /* populate metric components, it uses conformal metric */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        g00[ijk] = gConf_D0D0[ijk];
        g01[ijk] = gConf_D0D1[ijk];
        g02[ijk] = gConf_D0D2[ijk];
        g11[ijk] = gConf_D1D1[ijk];
        g12[ijk] = gConf_D1D2[ijk];
        g22[ijk] = gConf_D2D2[ijk];
      }
      adm[n]->g00 = g00;
      adm[n]->g01 = g01;
      adm[n]->g02 = g02;
      adm[n]->g11 = g11;
      adm[n]->g12 = g12;
      adm[n]->g22 = g22;
      
      /* surface integrals params */
      if ( n >= surface_index && 
           regex_search(".+right_BH_surrounding.+",patch->name) )
      {
        adm[n]->surface_integration_flg = 1;
        adm[n]->Z_surface = 1;
        adm[n]->K = 0;
        n_conformal_metric_surrounding(adm[n],_c_);
      }
    }
    obs->M = obs_ADM_mass;
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"ADM(M)|BH"))
  {  
    Patch_T **patches = 0,*patch = 0;
    struct items_S **adm = 0;
    unsigned p = 0;
    unsigned n,N,ijk,nn;
    
    N = 6/* 6 sides for surroundings for area integral */;
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
    /* surroundings for surface integrals. */
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
      
      READ_v(gConf_D2D2)
      READ_v(gConf_D0D2)
      READ_v(gConf_D0D0)
      READ_v(gConf_D0D1)
      READ_v(gConf_D1D2)
      READ_v(gConf_D1D1)
      READ_v(psi);

      adm[n]->patch = patch;
      /* populate metric components, it uses conformal metric */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*gConf_D0D0[ijk];
        g01[ijk] = psi4*gConf_D0D1[ijk];
        g02[ijk] = psi4*gConf_D0D2[ijk];
        g11[ijk] = psi4*gConf_D1D1[ijk];
        g12[ijk] = psi4*gConf_D1D2[ijk];
        g22[ijk] = psi4*gConf_D2D2[ijk];
      }
      adm[n]->g00 = g00;
      adm[n]->g01 = g01;
      adm[n]->g02 = g02;
      adm[n]->g11 = g11;
      adm[n]->g12 = g12;
      adm[n]->g22 = g22;
      
      adm[n]->surface_integration_flg = 1;
      adm[n]->Z_surface = 1;
      adm[n]->K = 0;
      n_physical_metric_surrounding(adm[n],_c_);
    }
    obs->M = obs_BH_ADM_mass;
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"ADM(M)|NS"))
  {  
    Patch_T **patches = 0,*patch = 0;
    struct items_S **adm = 0;
    unsigned p = 0;
    unsigned n,N,ijk,nn;
    
    N = 7/* 6 Ns patches for vol integral */;
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
    /* vol integrals. */
    patches[p++] = GetPatch("left_NS_up",grid);
    patches[p++] = GetPatch("left_NS_down",grid);
    patches[p++] = GetPatch("left_NS_left",grid);
    patches[p++] = GetPatch("left_NS_right",grid);
    patches[p++] = GetPatch("left_NS_back",grid);
    patches[p++] = GetPatch("left_NS_front",grid);
    patches[p++] = GetPatch("left_central_box",grid);
    
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
      
      READ_v(gConf_D2D2)
      READ_v(gConf_D0D2)
      READ_v(gConf_D0D0)
      READ_v(gConf_D0D1)
      READ_v(gConf_D1D2)
      READ_v(gConf_D1D1)
      
      adm[n]->patch = patch;
      /* populate metric components, it uses conformal metric */ 
      for (ijk = 0; ijk < nn; ++ijk)
      {
        g00[ijk] = gConf_D0D0[ijk];
        g01[ijk] = gConf_D0D1[ijk];
        g02[ijk] = gConf_D0D2[ijk];
        g11[ijk] = gConf_D1D1[ijk];
        g12[ijk] = gConf_D1D2[ijk];
        g22[ijk] = gConf_D2D2[ijk];
      }
      adm[n]->g00 = g00;
      adm[n]->g01 = g01;
      adm[n]->g02 = g02;
      adm[n]->g11 = g11;
      adm[n]->g12 = g12;
      adm[n]->g22 = g22;
      
    }
    obs->M = obs_ADM_mass;
    free(patches);
  }
  else
    Error1("There is no such '%s' plan.\n",obs->quantity);
  
}

/* populating normal outward vector for surrounding patches according to the given dir 
// NOTE: the normaliztion is respect to the physical metric gamma_{ij} */
static void n_physical_metric_surrounding(struct items_S *const adm,const Dd_T dir)
{
  Patch_T *const patch = adm->patch;
  const unsigned nn = patch->nn;
  double *n_U0 = alloc_double(nn);
  double *n_U1 = alloc_double(nn);
  double *n_U2 = alloc_double(nn);
  unsigned ijk;
  
  READ_v(gConf_D2D2)
  READ_v(gConf_D0D2)
  READ_v(gConf_D0D0)
  READ_v(gConf_D0D1)
  READ_v(gConf_D1D2)
  READ_v(gConf_D1D1)
  READ_v(psi);
    
  for (ijk = 0; ijk < nn; ++ijk)
  {
    n_U0[ijk] = dq2_dq1(patch,dir,_x_,ijk);
    n_U1[ijk] = dq2_dq1(patch,dir,_y_,ijk);
    n_U2[ijk] = dq2_dq1(patch,dir,_z_,ijk);
    
    /* normalization */
    double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
    double Norm2 = 
psi4*(gConf_D0D0[ijk]*pow(n_U0[ijk], 2) + 2.0*gConf_D0D1[ijk]*
n_U0[ijk]*n_U1[ijk] + 2.0*gConf_D0D2[ijk]*n_U0[ijk]*n_U2[ijk] +
gConf_D1D1[ijk]*pow(n_U1[ijk], 2) + 2.0*gConf_D1D2[ijk]*n_U1[ijk]*
n_U2[ijk] + gConf_D2D2[ijk]*pow(n_U2[ijk], 2));

    double Norm = sqrt(Norm2);
    
    n_U0[ijk] /= Norm;
    n_U1[ijk] /= Norm;
    n_U2[ijk] /= Norm;
    
  }
  adm->n_U0 = n_U0;
  adm->n_U1 = n_U1;
  adm->n_U2 = n_U2;
}

/* populating normal outward vector for surrounding patches according to the given dir 
// NOTE: the normaliztion is respect to the conformal metric gamma_{ij} */
static void n_conformal_metric_surrounding(struct items_S *const adm,const Dd_T dir)
{
  Patch_T *const patch = adm->patch;
  const unsigned nn = patch->nn;
  double *n_U0 = alloc_double(nn);
  double *n_U1 = alloc_double(nn);
  double *n_U2 = alloc_double(nn);
  unsigned ijk;
  
  READ_v(gConf_D2D2)
  READ_v(gConf_D0D2)
  READ_v(gConf_D0D0)
  READ_v(gConf_D0D1)
  READ_v(gConf_D1D2)
  READ_v(gConf_D1D1)
    
  for (ijk = 0; ijk < nn; ++ijk)
  {
    n_U0[ijk] = dq2_dq1(patch,dir,_x_,ijk);
    n_U1[ijk] = dq2_dq1(patch,dir,_y_,ijk);
    n_U2[ijk] = dq2_dq1(patch,dir,_z_,ijk);
    
    /* normalization */
    double Norm2 = 
(gConf_D0D0[ijk]*pow(n_U0[ijk], 2) + 2.0*gConf_D0D1[ijk]*
n_U0[ijk]*n_U1[ijk] + 2.0*gConf_D0D2[ijk]*n_U0[ijk]*n_U2[ijk] +
gConf_D1D1[ijk]*pow(n_U1[ijk], 2) + 2.0*gConf_D1D2[ijk]*n_U1[ijk]*
n_U2[ijk] + gConf_D2D2[ijk]*pow(n_U2[ijk], 2));

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
void obs_free(Observable_T *obs)
{
  if (!obs)
    return;
    
  if (!strcmp_i(obs->grid->kind,"BBN_CubedSpherical_grid"))
    Error0(NO_OPTION);
  
  struct items_S **adm = obs->items;
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
      
      if (_Ind("ADM_integrand_xiP_U0") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiP_U0);
        REMOVE_FIELD(ADM_integrand_xiP_U0);
      }
      if (_Ind("ADM_integrand_xiP_U2") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiP_U2);
        REMOVE_FIELD(ADM_integrand_xiP_U2);
      }
      if (_Ind("ADM_integrand_xiP_U1") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiP_U1);
        REMOVE_FIELD(ADM_integrand_xiP_U1);
      }
      
      if (_Ind("ADM_integrand_xiG_U0") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiG_U0);
        REMOVE_FIELD(ADM_integrand_xiG_U0);
      }
      if (_Ind("ADM_integrand_xiG_U2") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiG_U2);
        REMOVE_FIELD(ADM_integrand_xiG_U2);
      }
      if (_Ind("ADM_integrand_xiG_U1") >= 0)
      {
        DECLARE_FIELD(ADM_integrand_xiG_U1);
        REMOVE_FIELD(ADM_integrand_xiG_U1);
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
static double ADM_momentum_x_BHNS_CS(Observable_T *const obs)
{
  double Px = 0;
  struct items_S **const adm = obs->items;
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
static double ADM_momentum_y_BHNS_CS(Observable_T *const obs)
{
  double Py = 0;
  struct items_S **const adm = obs->items;
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
static double ADM_momentum_z_BHNS_CS(Observable_T *const obs)
{
  double Pz = 0;
  struct items_S **const adm = obs->items;
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
static double ADM_angular_momentum_z_BHNS_CS(Observable_T *const obs)
{
  double Jz = 0;
  struct items_S **const adm = obs->items;
  const unsigned N = obs->Nitems;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    if (!adm[p]->surface_integration_flg)
      continue;
    
    Patch_T *patch = adm[p]->patch;
    Field_T *xiPz = patch->pool[Ind("ADM_integrand_xiP_U2")];
    
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
        
    I->Spectral->f = xiPz;
    plan_integration(I);
    Jz += execute_integration(I);
    
    free_integration(I);
  }
  
  /* volume integration */
  if (VOLUME_INTEGRAL)
  for(p = 0; p < N; ++p)
  {
    if (adm[p]->surface_integration_flg)
      continue;

    Patch_T *patch = adm[p]->patch;
    Field_T *xiGz = patch->pool[Ind("ADM_integrand_xiG_U2")];
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
    I->Spectral->f = xiGz;
    plan_integration(I);
    Jz -= execute_integration(I);
    
    free_integration(I);
  }
  
  Jz /= (8*M_PI);
  return Jz;
}

/* calculating ADM angular momentum in x component */
static double ADM_angular_momentum_x_BHNS_CS(Observable_T *const obs)
{
  double Jx = 0;
  struct items_S **const adm = obs->items;
  const unsigned N = obs->Nitems;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
   if (!adm[p]->surface_integration_flg)
      continue;
   
    Patch_T *patch = adm[p]->patch;
    Field_T *xiPx = patch->pool[Ind("ADM_integrand_xiP_U0")];
    
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
    
    I->Spectral->f = xiPx;
    plan_integration(I);
    Jx += execute_integration(I);
    
    free_integration(I);
  }
  
  /* volume integration */
  if (VOLUME_INTEGRAL)
  for(p = 0; p < N; ++p)
  {
    if (adm[p]->surface_integration_flg)
      continue;
   
    Patch_T *patch = adm[p]->patch;
    Field_T *xiGx = patch->pool[Ind("ADM_integrand_xiG_U0")];
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
    I->Spectral->f = xiGx;
    plan_integration(I);
    Jx -= execute_integration(I);
    
    free_integration(I);
  }
  
  Jx /= (8*M_PI);
  return Jx;
}

/* calculating ADM angular momentum in y component */
static double ADM_angular_momentum_y_BHNS_CS(Observable_T *const obs)
{
  double Jy = 0;
  struct items_S **const adm = obs->items;
  const unsigned N = obs->Nitems;
  Integration_T *I;
  unsigned p;
  assert(N);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    if (!adm[p]->surface_integration_flg)
      continue;
    
    Patch_T *patch = adm[p]->patch;
    Field_T *xiPy = patch->pool[Ind("ADM_integrand_xiP_U1")];
    
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
        
    I->Spectral->f = xiPy;
    plan_integration(I);
    Jy += execute_integration(I);
    
    free_integration(I);
  }
  
  /* volume integration */
  if (VOLUME_INTEGRAL)
  for(p = 0; p < N; ++p)
  {
    if (adm[p]->surface_integration_flg)
      continue;
   
    Patch_T *patch = adm[p]->patch;
    Field_T *xiGy = patch->pool[Ind("ADM_integrand_xiG_U1")];
    
    I  = init_integration();
    I->type = "Integral{f(x)dV},Spectral";
    I->g00 = adm[p]->g00;
    I->g01 = adm[p]->g01;
    I->g02 = adm[p]->g02;
    I->g11 = adm[p]->g11;
    I->g12 = adm[p]->g12;
    I->g22 = adm[p]->g22;
    
    I->Spectral->f = xiGy;
    plan_integration(I);
    Jy -= execute_integration(I);
    
    free_integration(I);
  }
  
  Jy /= (8*M_PI);
  return Jy;
}

/* approximate spin using : S_a = \frac{1}{8\pi}\oint{\xi_{ai} K^{ij}ds^{2}_j} */
void obs_define_spin_integral(double S[3],Grid_T *const grid,const char *const kind)
{
  if (!strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
    Error0(NO_OPTION);
    
  const unsigned N     = 6;  
  double obj_center[3] = {0};
  Patch_T *patches[N];
  unsigned p = 0;
  
  S[0] = S[1] = S[2] = 0;
  
  /* NS spins */
  if (strcmp_i(kind,"NS"))
  {
     /* surroundings for surface integrals */
    patches[p++] = GetPatch("left_NS_surrounding_up",grid);
    patches[p++] = GetPatch("left_NS_surrounding_down",grid);
    patches[p++] = GetPatch("left_NS_surrounding_left",grid);
    patches[p++] = GetPatch("left_NS_surrounding_right",grid);
    patches[p++] = GetPatch("left_NS_surrounding_back",grid);
    patches[p++] = GetPatch("left_NS_surrounding_front",grid);
    obj_center[0]= Pgetd("NS_center_x");
    obj_center[1]= Pgetd("NS_center_y");
    obj_center[2]= Pgetd("NS_center_z");
  }
  /* BH spins */
  if (strcmp_i(kind,"BH"))
  {
    /* surroundings for surface integrals */
    patches[p++] = GetPatch("right_BH_surrounding_up",grid);
    patches[p++] = GetPatch("right_BH_surrounding_down",grid);
    patches[p++] = GetPatch("right_BH_surrounding_left",grid);
    patches[p++] = GetPatch("right_BH_surrounding_right",grid);
    patches[p++] = GetPatch("right_BH_surrounding_back",grid);
    patches[p++] = GetPatch("right_BH_surrounding_front",grid);
    obj_center[0]= Pgetd("BH_center_x");
    obj_center[1]= Pgetd("BH_center_y");
    obj_center[2]= Pgetd("BH_center_z");
  }  
  assert(p==N);
  
  /* carry out the integral for each patch */
  for (p = 0; p < N; ++p)
  {
    Patch_T *patch = patches[p];
    Integration_T *I = 0;
    unsigned nn  = patch->nn;
    struct items_S normal[1] = {0};
    const double *n_comp[3];
    double *g00 = alloc_double(nn);
    double *g01 = alloc_double(nn);
    double *g02 = alloc_double(nn);
    double *g11 = alloc_double(nn);
    double *g12 = alloc_double(nn);
    double *g22 = alloc_double(nn);
    unsigned ijk;
    
    READ_v(gConf_D2D2)
    READ_v(gConf_D0D2)
    READ_v(gConf_D0D0)
    READ_v(gConf_D0D1)
    READ_v(gConf_D1D2)
    READ_v(gConf_D1D1)
    READ_v(psi);
    
    /* populate metric components */ 
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
      g00[ijk] = psi4*gConf_D0D0[ijk];
      g01[ijk] = psi4*gConf_D0D1[ijk];
      g02[ijk] = psi4*gConf_D0D2[ijk];
      g11[ijk] = psi4*gConf_D1D1[ijk];
      g12[ijk] = psi4*gConf_D1D2[ijk];
      g22[ijk] = psi4*gConf_D2D2[ijk];
    }
    
    normal->patch = patch;
    n_physical_metric_surrounding(normal,_c_);
    n_comp[0] = normal->n_U0;
    n_comp[1] = normal->n_U1;
    n_comp[2] = normal->n_U2;
    obs_populate_spin_integrands_Campanelli(patch,obj_center,n_comp);
    
    /* surface integral */
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->g00 = g00;
    I->g01 = g01;
    I->g02 = g02;
    I->g11 = g11;
    I->g12 = g12;
    I->g22 = g22;
    I->Spectral->Z_surface = 1;
    I->Spectral->K         = 0;
    
    I->Spectral->f = patch->pool[Ind("SPIN_integrand_U0")];
    plan_integration(I);
    S[0] += execute_integration(I);
    
    I->Spectral->f = patch->pool[Ind("SPIN_integrand_U1")];
    plan_integration(I);
    S[1] += execute_integration(I);
    
    I->Spectral->f = patch->pool[Ind("SPIN_integrand_U2")];
    plan_integration(I);
    S[2] += execute_integration(I);
    
    /* free */
    DECLARE_FIELD(SPIN_integrand_U0);
    REMOVE_FIELD(SPIN_integrand_U0);
    DECLARE_FIELD(SPIN_integrand_U1);
    REMOVE_FIELD(SPIN_integrand_U1);
    DECLARE_FIELD(SPIN_integrand_U2);
    REMOVE_FIELD(SPIN_integrand_U2);
    free_integration(I);    
    _free(normal->n_U0);
    _free(normal->n_U1);
    _free(normal->n_U2);
    _free(g00);
    _free(g01);
    _free(g02);
    _free(g11);
    _free(g12);
    _free(g22);
  }
  S[0] /= (8*M_PI);
  S[1] /= (8*M_PI);
  S[2] /= (8*M_PI);
}


/* approximate spin using : S_a = \frac{1}{8\pi}\oint{\xi_{ai} K^{ij}ds^{2}_j} */
void 
obs_define_spin_akv
  (
  double S[3]/* spin Sx,Sy,Sz */,
  Grid_T *const grid/* grid */,
  const char *const kind/* "NS" or "BH" */
  )
{
  FUNC_TIC
  
  if (!strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
    Error0(NO_OPTION);
    
  const unsigned N     = 6;  
  Patch_T *patches[N];
  unsigned p = 0;
  
  S[0] = S[1] = S[2] = 0;
  
  /* NS spins */
  if (strcmp_i(kind,"NS"))
  {
     /* surroundings for surface integrals */
    patches[p++] = GetPatch("left_NS_surrounding_up",grid);
    patches[p++] = GetPatch("left_NS_surrounding_down",grid);
    patches[p++] = GetPatch("left_NS_surrounding_left",grid);
    patches[p++] = GetPatch("left_NS_surrounding_right",grid);
    patches[p++] = GetPatch("left_NS_surrounding_back",grid);
    patches[p++] = GetPatch("left_NS_surrounding_front",grid);
  }
  /* BH spins */
  if (strcmp_i(kind,"BH"))
  {
    /* surroundings for surface integrals */
    patches[p++] = GetPatch("right_BH_surrounding_up",grid);
    patches[p++] = GetPatch("right_BH_surrounding_down",grid);
    patches[p++] = GetPatch("right_BH_surrounding_left",grid);
    patches[p++] = GetPatch("right_BH_surrounding_right",grid);
    patches[p++] = GetPatch("right_BH_surrounding_back",grid);
    patches[p++] = GetPatch("right_BH_surrounding_front",grid);
  }  
  assert(p==N);
  
  /* carry out the integral for each patch */
  for (p = 0; p < N; ++p)
  {
    Patch_T *patch = patches[p];
    Integration_T *I = 0;
    unsigned nn  = patch->nn;
    struct items_S normal[1] = {0};
    const double *n_comp[3];
    double *g00 = alloc_double(nn);
    double *g01 = alloc_double(nn);
    double *g02 = alloc_double(nn);
    double *g11 = alloc_double(nn);
    double *g12 = alloc_double(nn);
    double *g22 = alloc_double(nn);
    unsigned ijk;
    
    READ_v(gConf_D2D2)
    READ_v(gConf_D0D2)
    READ_v(gConf_D0D0)
    READ_v(gConf_D0D1)
    READ_v(gConf_D1D2)
    READ_v(gConf_D1D1)
    READ_v(psi);
    
    /* populate metric components */ 
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
      g00[ijk] = psi4*gConf_D0D0[ijk];
      g01[ijk] = psi4*gConf_D0D1[ijk];
      g02[ijk] = psi4*gConf_D0D2[ijk];
      g11[ijk] = psi4*gConf_D1D1[ijk];
      g12[ijk] = psi4*gConf_D1D2[ijk];
      g22[ijk] = psi4*gConf_D2D2[ijk];
    }
    
    normal->patch = patch;
    n_physical_metric_surrounding(normal,_c_);
    n_comp[0] = normal->n_U0;
    n_comp[1] = normal->n_U1;
    n_comp[2] = normal->n_U2;
    obs_populate_spin_integrands_akv(patch,n_comp);
    
    /* surface integral */
    I  = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->g00 = g00;
    I->g01 = g01;
    I->g02 = g02;
    I->g11 = g11;
    I->g12 = g12;
    I->g22 = g22;
    I->Spectral->Z_surface = 1;
    I->Spectral->K         = 0;
    
    I->Spectral->f = patch->pool[Ind("SPIN_integrand_U0")];
    plan_integration(I);
    S[0] += execute_integration(I);
    
    I->Spectral->f = patch->pool[Ind("SPIN_integrand_U1")];
    plan_integration(I);
    S[1] += execute_integration(I);
    
    I->Spectral->f = patch->pool[Ind("SPIN_integrand_U2")];
    plan_integration(I);
    S[2] += execute_integration(I);
    
    /* free */
    DECLARE_FIELD(SPIN_integrand_U0);
    REMOVE_FIELD(SPIN_integrand_U0);
    DECLARE_FIELD(SPIN_integrand_U1);
    REMOVE_FIELD(SPIN_integrand_U1);
    DECLARE_FIELD(SPIN_integrand_U2);
    REMOVE_FIELD(SPIN_integrand_U2);
    free_integration(I);    
    _free(normal->n_U0);
    _free(normal->n_U1);
    _free(normal->n_U2);
    _free(g00);
    _free(g01);
    _free(g02);
    _free(g11);
    _free(g12);
    _free(g22);
  }
  S[0] /= (8*M_PI);
  S[1] /= (8*M_PI);
  S[2] /= (8*M_PI);
  
  FUNC_TOC
}

/* approximate spin using : S = J - RxP */
void obs_define_spin_JRP(double S[3],Grid_T *const grid,const char *const kind)
{
  double J[3] = {0,0,0};
  double R[3] = {0,0,0};
  double P[3] = {0,0,0};
  
  /* NS spins */
  if (strcmp_i(kind,"NS"))
  {
    obs_Rc_NS(R,grid);
    P[0] = Pgetd("NS_Px_ADM");
    P[1] = Pgetd("NS_Py_ADM");
    P[2] = Pgetd("NS_Pz_ADM");
    J[0] = Pgetd("NS_Jx_ADM");
    J[1] = Pgetd("NS_Jy_ADM");
    J[2] = Pgetd("NS_Jz_ADM");
  }
  else if (strcmp_i(kind,"BH"))
  {
    obs_Rc_BH(R,grid);
    P[0] = Pgetd("BH_Px_ADM");
    P[1] = Pgetd("BH_Py_ADM");
    P[2] = Pgetd("BH_Pz_ADM");
    J[0] = Pgetd("BH_Jx_ADM");
    J[1] = Pgetd("BH_Jy_ADM");
    J[2] = Pgetd("BH_Jz_ADM");
  }
  else
    Error0(NO_OPTION);

  S[0] = J[0] - (-P[1]*R[2] + P[2]*R[1]);
  S[1] = J[1] - (P[0]*R[2] - P[2]*R[0]);
  S[2] = J[2] - (-P[0]*R[1] + P[1]*R[0]);
}


