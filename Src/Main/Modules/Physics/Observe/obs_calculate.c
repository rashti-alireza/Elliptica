/*
// Alireza Rashti
// September 2019
*/


#include "obs_calculate.h"
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
void obs_calculate(Observe_T *const obs)
{
  Grid_T *const grid = obs->grid;
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS)
  {
      
  if (strcmp_i(obs->quantity,"ADM(P,J)|BHNS"))
  {  
    Patch_T **patches = 0;
    Patch_T *patch    = 0;
    const char *region = "outermost,filling_box,"
                         "NS_around_OB,BH_around_OB";
    struct items_S **adm = 0;
    Uint n,N,ijk,nn;
    
    /* first collect all of the patches required */
    patches = collect_patches(grid,region,&N);
    
    /* alloc memory for all patches */
    adm = calloc(N,sizeof(*adm));
    IsNull(adm);
    /* this is where we link to obs struct */
    obs->items = adm;
    obs->Nitems = N;
    
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
      if (IsItCovering(patch,"NS_around_OB") ||
          IsItCovering(patch,"BH_around_OB")  )
      {
        adm[n]->surface_integration_flg = 1;
        adm[n]->Z_surface = 1;
        adm[n]->K = patch->n[2]-1;
        n_physical_metric_around(adm[n],_c_);
      }
    }
    
    obs_populate_ADM_integrand_PdS_GdV_binary(obs);
    obs->ret[0] = ADM_momentum_x_BHNS_CS(obs);
    obs->ret[1] = ADM_momentum_y_BHNS_CS(obs);
    obs->ret[2] = ADM_momentum_z_BHNS_CS(obs);
    obs->ret[3] = ADM_angular_momentum_x_BHNS_CS(obs);
    obs->ret[4] = ADM_angular_momentum_y_BHNS_CS(obs);
    obs->ret[5] = ADM_angular_momentum_z_BHNS_CS(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"ADM(P,J)|NS"))
  {  
    Patch_T **patches = 0;
    Patch_T *patch    = 0;
    struct items_S **adm = 0;
    const char *region = "NS_OB";
    Uint n,N,ijk,nn;
    
    /* first collect all of the patches required */
    patches = collect_patches(grid,region,&N);
    
    /* alloc memory for all patches */
    adm = calloc(N,sizeof(*adm));
    IsNull(adm);
    /* this is where we link to obs struct */
    obs->items = adm;
    obs->Nitems = N;
    
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
      adm[n]->K = patch->n[2]-1;
      n_physical_metric_around(adm[n],_c_);
    }
    
    obs_populate_ADM_integrand_PdS_GdV_single(obs);
    obs->ret[0] = ADM_momentum_x_BHNS_CS(obs);
    obs->ret[1] = ADM_momentum_y_BHNS_CS(obs);
    obs->ret[2] = ADM_momentum_z_BHNS_CS(obs);
    obs->ret[3] = ADM_angular_momentum_x_BHNS_CS(obs);
    obs->ret[4] = ADM_angular_momentum_y_BHNS_CS(obs);
    obs->ret[5] = ADM_angular_momentum_z_BHNS_CS(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"ADM(P,J)|BH"))
  {  
    Patch_T **patches = 0;
    Patch_T *patch    = 0;
    struct items_S **adm = 0;
    const char *region = "BH_around_IB";
    Uint n,N,ijk,nn;
    
    /* first collect all of the patches required */
    patches = collect_patches(grid,region,&N);
    
    /* alloc memory for all patches */
    adm = calloc(N,sizeof(*adm));
    IsNull(adm);
    /* this is where we link to obs struct */
    obs->items = adm;
    obs->Nitems = N;
      
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
      n_physical_metric_around(adm[n],_c_);
    }
    obs_populate_ADM_integrand_PdS_GdV_single(obs);
    obs->ret[0] = ADM_momentum_x_BHNS_CS(obs);
    obs->ret[1] = ADM_momentum_y_BHNS_CS(obs);
    obs->ret[2] = ADM_momentum_z_BHNS_CS(obs);
    obs->ret[3] = ADM_angular_momentum_x_BHNS_CS(obs);
    obs->ret[4] = ADM_angular_momentum_y_BHNS_CS(obs);
    obs->ret[5] = ADM_angular_momentum_z_BHNS_CS(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"Kommar(M)|BH"))
  {  
    Patch_T **patches = 0;
    Patch_T *patch    = 0;
    struct items_S **kommar = 0;
    const char *region = "BH_around_IB";
    Uint n,N,ijk,nn;
    
    /* first collect all of the patches required */
    patches = collect_patches(grid,region,&N);
    
    /* alloc memory for all patches */
    kommar = calloc(N,sizeof(*kommar));
    IsNull(kommar);
    /* this is where we link to obs struct */
    obs->items = kommar;
    obs->Nitems = N;
        
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
      n_physical_metric_around(kommar[n],_c_);
    }
    obs->ret[0] = obs_Kommar_mass(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"Kommar(M)|NS"))
  {  
    Patch_T **patches = 0;
    Patch_T *patch    = 0;
    struct items_S **kommar = 0;
    const char *region = "NS_OB";
    Uint n,N,ijk,nn;
    
    /* first collect all of the patches required */
    patches = collect_patches(grid,region,&N);
    
    /* alloc memory for all patches */
    kommar = calloc(N,sizeof(*kommar));
    IsNull(kommar);
    /* this is where we link to obs struct */
    obs->items = kommar;
    obs->Nitems = N;
        
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
      n_physical_metric_around(kommar[n],_c_);
    }
    obs->ret[0] = obs_Kommar_mass(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"Kommar(M)|BHNS"))
  {  
    Patch_T **patches = 0;
    Patch_T *patch    = 0;
    struct items_S **kommar = 0;
    Uint p = 0;
    Uint n,N,ijk,nn;
    
    N = 10/* 10 sides for arounds */;
        
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
    /* arounds for surface integrals */
    patches[p++] = GetPatch("left_NS_around_up",grid);
    patches[p++] = GetPatch("left_NS_around_down",grid);
    patches[p++] = GetPatch("left_NS_around_left",grid);
    patches[p++] = GetPatch("left_NS_around_back",grid);
    patches[p++] = GetPatch("left_NS_around_front",grid);
    
    patches[p++] = GetPatch("right_BH_around_up",grid);
    patches[p++] = GetPatch("right_BH_around_down",grid);
    patches[p++] = GetPatch("right_BH_around_right",grid);
    patches[p++] = GetPatch("right_BH_around_back",grid);
    patches[p++] = GetPatch("right_BH_around_front",grid);

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
      n_physical_metric_around(kommar[n],_c_);
    }
    obs->ret[0] = obs_Kommar_mass(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"ADM(M)|BHNS"))
  {  
    Patch_T **patches1 = 0,**patches2 = 0;
    Patch_T *patch    = 0;
    const char *region;
    struct items_S **adm = 0;
    Uint n,N1,N2,ijk,nn;
    
    /* volume part */
    region   = "outermost,filling_box,NS,NS_around,BH_around";
    patches1 = collect_patches(grid,region,&N1);
    /* surface part */
    region   = "BH_around_IB";
    patches2 = collect_patches(grid,region,&N2);
    
    /* alloc memory for all patches */
    adm = calloc((N1+N2),sizeof(*adm));
    IsNull(adm);
    /* this is where we link to obs struct */
    obs->items = adm;
    obs->Nitems = N1+N2;
    
    /* fill ADM struct for each patch volume part */
    for (n = 0; n < N1; ++n)
    {
      adm[n] = calloc(1,sizeof(*adm[n]));
      IsNull(adm[n]);
      patch = patches1[n];
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
    free(patches1);
    
    for (n = N1; n < N2; ++n)
    {
      adm[n] = calloc(1,sizeof(*adm[n]));
      IsNull(adm[n]);
      patch = patches2[n-N1];
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
      
      adm[n]->surface_integration_flg = 1;
      adm[n]->Z_surface = 1;
      adm[n]->K = 0;
      n_conformal_metric_around(adm[n],_c_);
    }
    obs->ret[0] = obs_ADM_mass(obs);
    free(patches2);
  }
  else if (strcmp_i(obs->quantity,"ADM(M)|BH"))
  {  
    Patch_T **patches = 0;
    Patch_T *patch    = 0;
    struct items_S **adm = 0;
    const char *region = "BH_around_IB";
    Uint n,N,ijk,nn;
    
    /* first collect all of the patches required */
    patches = collect_patches(grid,region,&N);
    
    /* alloc memory for all patches */
    adm = calloc(N,sizeof(*adm));
    IsNull(adm);
    /* this is where we link to obs struct */
    obs->items = adm;
    obs->Nitems = N;
    
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
      n_physical_metric_around(adm[n],_c_);
    }
    obs->ret[0] = obs_BH_ADM_mass(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"ADM(M)|NS"))
  {  
    Patch_T **patches = 0;
    Patch_T *patch    = 0;
    struct items_S **adm = 0;
    const char *region = "NS";
    Uint n,N,ijk,nn;
    
    /* first collect all of the patches required */
    patches = collect_patches(grid,region,&N);
  
    /* alloc memory for all patches */
    adm = calloc(N,sizeof(*adm));
    IsNull(adm);
    /* this is where we link to obs struct */
    obs->items = adm;
    obs->Nitems = N;
    
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
    obs->ret[0] = obs_ADM_mass(obs);
    free(patches);
  }
  else if (strcmp_i(obs->quantity,"CM|BH"))
  {
    Rc_BH(obs);
  }
  else if (strcmp_i(obs->quantity,"CM|NS"))
  {
    obs_Rc_NS(obs);
  }
  else if (strstr_i(obs->quantity,"Spin|JRP|"))
  {
    define_spin_JRP(obs);
  }
  else if (strstr_i(obs->quantity,"Spin|Campanelli|"))
  {
    define_spin_campanelli(obs);
  }
  else if (strstr_i(obs->quantity,"Spin|AKV|"))
  {
    define_spin_akv(obs);
  }
  else if (strcmp_i(obs->quantity,"Irreducible(M)|BH"))
  {
    obs_BH_irreducible_mass_CS(obs);
  }
  else
    Error1("There is no such '%s' plan.\n",obs->quantity);
    
  }/* if (grid->kind == Grid_SplitCubedSpherical_BHNS) */
  else
    Error0(NO_OPTION);
}

/* populating normal outward vector for around patches according to the given dir 
// NOTE: the normaliztion is respect to the physical metric gamma_{ij} */
static void n_physical_metric_around(struct items_S *const adm,const Dd_T dir)
{
  Patch_T *const patch = adm->patch;
  const Uint nn = patch->nn;
  double *n_U0 = alloc_double(nn);
  double *n_U1 = alloc_double(nn);
  double *n_U2 = alloc_double(nn);
  Uint ijk;
  
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

/* populating normal outward vector for around patches according to the given dir 
// NOTE: the normaliztion is respect to the conformal metric gamma_{ij} */
static void n_conformal_metric_around(struct items_S *const adm,const Dd_T dir)
{
  Patch_T *const patch = adm->patch;
  const Uint nn = patch->nn;
  double *n_U0 = alloc_double(nn);
  double *n_U1 = alloc_double(nn);
  double *n_U2 = alloc_double(nn);
  Uint ijk;
  
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

/* calculating ADM momentum in x component */
static double ADM_momentum_x_BHNS_CS(Observe_T *const obs)
{
  double Px = 0;
  struct items_S **const adm = obs->items;
  const Uint N = obs->Nitems;
  Integration_T *I;
  Uint p;
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
    I->Spectral->f = patch->fields[Ind("ADM_integrand_P_U0")];
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
    I->Spectral->f = patch->fields[Ind("ADM_integrand_G_U0")];
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
static double ADM_momentum_y_BHNS_CS(Observe_T *const obs)
{
  double Py = 0;
  struct items_S **const adm = obs->items;
  const Uint N = obs->Nitems;
  Integration_T *I;
  Uint p;
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
    I->Spectral->f = patch->fields[Ind("ADM_integrand_P_U1")];
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
    I->Spectral->f = patch->fields[Ind("ADM_integrand_G_U1")];
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
static double ADM_momentum_z_BHNS_CS(Observe_T *const obs)
{
  double Pz = 0;
  struct items_S **const adm = obs->items;
  const Uint N = obs->Nitems;
  Integration_T *I;
  Uint p;
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
    I->Spectral->f = patch->fields[Ind("ADM_integrand_P_U2")];
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
    I->Spectral->f = patch->fields[Ind("ADM_integrand_G_U2")];
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
static double ADM_angular_momentum_z_BHNS_CS(Observe_T *const obs)
{
  double Jz = 0;
  struct items_S **const adm = obs->items;
  const Uint N = obs->Nitems;
  Integration_T *I;
  Uint p;
  assert(N);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    if (!adm[p]->surface_integration_flg)
      continue;
    
    Patch_T *patch = adm[p]->patch;
    Field_T *xiPz = patch->fields[Ind("ADM_integrand_xiP_U2")];
    
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
    Field_T *xiGz = patch->fields[Ind("ADM_integrand_xiG_U2")];
    
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
static double ADM_angular_momentum_x_BHNS_CS(Observe_T *const obs)
{
  double Jx = 0;
  struct items_S **const adm = obs->items;
  const Uint N = obs->Nitems;
  Integration_T *I;
  Uint p;
  assert(N);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
   if (!adm[p]->surface_integration_flg)
      continue;
   
    Patch_T *patch = adm[p]->patch;
    Field_T *xiPx = patch->fields[Ind("ADM_integrand_xiP_U0")];
    
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
    Field_T *xiGx = patch->fields[Ind("ADM_integrand_xiG_U0")];
    
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
static double ADM_angular_momentum_y_BHNS_CS(Observe_T *const obs)
{
  double Jy = 0;
  struct items_S **const adm = obs->items;
  const Uint N = obs->Nitems;
  Integration_T *I;
  Uint p;
  assert(N);
  
  /* surface integration */
  for(p = 0; p < N; ++p)
  {
    if (!adm[p]->surface_integration_flg)
      continue;
    
    Patch_T *patch = adm[p]->patch;
    Field_T *xiPy = patch->fields[Ind("ADM_integrand_xiP_U1")];
    
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
    Field_T *xiGy = patch->fields[Ind("ADM_integrand_xiG_U1")];
    
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
static void define_spin_campanelli(Observe_T *const obs)
{
  Physics_T *const phys = obs->phys;
  double *const S = obs->ret;
  Patch_T **patches = 0;
  double obj_center[3] = {0};
  const char *region = 0;
  Uint N,p = 0;
  
  S[0] = S[1] = S[2] = 0;
  
  obj_center[0]= Getd("center_x");
  obj_center[1]= Getd("center_y");
  obj_center[2]= Getd("center_z");
  
  /* NS spins */
  if (phys->ctype == NS)
  {
    region = Ftype("NS_around_IB");
  }
  /* BH spins */
  else if (phys->ctype == BH)
  {
    region = Ftype("BH_around_IB");
  }
  else
  {
    Error0(NO_OPTION);
  }
  
  patches = collect_patches(phys->grid,region,&N);
  
  /* carry out the integral for each patch */
  for (p = 0; p < N; ++p)
  {
    Patch_T *patch = patches[p];
    Integration_T *I = 0;
    Uint nn  = patch->nn;
    struct items_S normal[1] = {0};
    const double *n_comp[3];
    double *g00 = alloc_double(nn);
    double *g01 = alloc_double(nn);
    double *g02 = alloc_double(nn);
    double *g11 = alloc_double(nn);
    double *g12 = alloc_double(nn);
    double *g22 = alloc_double(nn);
    Uint ijk;
    
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
    n_physical_metric_around(normal,_c_);
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
    
    I->Spectral->f = patch->fields[Ind("SPIN_integrand_U0")];
    plan_integration(I);
    S[0] += execute_integration(I);
    
    I->Spectral->f = patch->fields[Ind("SPIN_integrand_U1")];
    plan_integration(I);
    S[1] += execute_integration(I);
    
    I->Spectral->f = patch->fields[Ind("SPIN_integrand_U2")];
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
  
  _free(patches);
  
  S[0] /= (8*M_PI);
  S[1] /= (8*M_PI);
  S[2] /= (8*M_PI);
}


/* approximate spin using : S_a = \frac{1}{8\pi}\oint{\xi_{ai} K^{ij}ds^{2}_j} */
static void define_spin_akv(Observe_T *const obs)
{
  Physics_T *const phys = obs->phys;
  double *const S = obs->ret;
  Patch_T **patches = 0;
  const char *region = 0;
  Uint N,p = 0;
  
  S[0] = S[1] = S[2] = 0;
  
  /* NS spins */
  if (phys->ctype == NS)
  {
    region = Ftype("NS_around_IB");
  }
  /* BH spins */
  else if (phys->ctype == BH)
  {
    region = Ftype("BH_around_IB");
  }
  else
  {
    Error0(NO_OPTION);
  }
  
  patches = collect_patches(phys->grid,region,&N);

  
  /* carry out the integral for each patch */
  for (p = 0; p < N; ++p)
  {
    Patch_T *patch = patches[p];
    Integration_T *I = 0;
    Uint nn  = patch->nn;
    struct items_S normal[1] = {0};
    const double *n_comp[3];
    double *g00 = alloc_double(nn);
    double *g01 = alloc_double(nn);
    double *g02 = alloc_double(nn);
    double *g11 = alloc_double(nn);
    double *g12 = alloc_double(nn);
    double *g22 = alloc_double(nn);
    Uint ijk;
    
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
    n_physical_metric_around(normal,_c_);
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
    
    I->Spectral->f = patch->fields[Ind("SPIN_integrand_U0")];
    plan_integration(I);
    S[0] += execute_integration(I);
    
    I->Spectral->f = patch->fields[Ind("SPIN_integrand_U1")];
    plan_integration(I);
    S[1] += execute_integration(I);
    
    I->Spectral->f = patch->fields[Ind("SPIN_integrand_U2")];
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

/* approximate spin using : S = J - RxP */
static void define_spin_JRP(Observe_T *const obs)
{
  Physics_T *const phys = obs->phys;
  double *const S = obs->ret;
  double J[3] = {0,0,0};
  double R[3] = {0,0,0};
  double P[3] = {0,0,0};

  P[0] = Getd("Px_ADM");
  P[1] = Getd("Py_ADM");
  P[2] = Getd("Pz_ADM");
  J[0] = Getd("Jx_ADM");
  J[1] = Getd("Jy_ADM");
  J[2] = Getd("Jz_ADM");
  
  /* NS spins */
  if (phys->ctype == NS)
  {
    obs_Rc_NS(obs);
  }
  else if (phys->ctype == BH)
  {
    Rc_BH(obs);
  }
  else
    Error0(NO_OPTION);

  S[0] = J[0] - (-P[1]*R[2] + P[2]*R[1]);
  S[1] = J[1] - (P[0]*R[2] - P[2]*R[0]);
  S[2] = J[2] - (-P[0]*R[1] + P[1]*R[0]);
}


/* calculating physical center of BH to be used in spin calculations */
static void Rc_BH(Observe_T *const obs)
{
  Physics_T *const phys = obs->phys;
  Grid_T *const grid    = phys->grid;
  const double AH_area  = Getd("AH_area");
  const double x_CM = sysGetd("x_CM");
  const double y_CM = sysGetd("y_CM");
  const double z_CM = sysGetd("z_CM");
  double *const Rc  = obs->ret;
  Uint p;

  Rc[0] = 0;
  Rc[1] = 0;
  Rc[2] = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItCovering(patch,"BH_around_IB"))
      continue;
      
    Uint ijk;
    Uint nn = patch->nn;
    
    READ_v(gConf_D2D2)
    READ_v(gConf_D0D2)
    READ_v(gConf_D0D0)
    READ_v(gConf_D0D1)
    READ_v(gConf_D1D2)
    READ_v(gConf_D1D1)
    READ_v(psi);
    ADD_FIELD(Rc_integrandx)
    ADD_FIELD(Rc_integrandy)
    ADD_FIELD(Rc_integrandz)
    
    double *g00 = alloc_double(nn);
    double *g01 = alloc_double(nn);
    double *g02 = alloc_double(nn);
    double *g11 = alloc_double(nn);
    double *g12 = alloc_double(nn);
    double *g22 = alloc_double(nn);
    
    {/* local variables */
      REALLOC_v_WRITE_v(Rc_integrandx)
      REALLOC_v_WRITE_v(Rc_integrandy)
      REALLOC_v_WRITE_v(Rc_integrandz)

      FOR_ALL_POINTS(ijk,patch)
      {
        double x = patch->node[ijk]->x[0];
        double y = patch->node[ijk]->x[1];
        double z = patch->node[ijk]->x[2];
        double psi4 = Pow2(psi[ijk])*Pow2(psi[ijk]);
        g00[ijk] = psi4*gConf_D0D0[ijk];
        g01[ijk] = psi4*gConf_D0D1[ijk];
        g02[ijk] = psi4*gConf_D0D2[ijk];
        g11[ijk] = psi4*gConf_D1D1[ijk];
        g12[ijk] = psi4*gConf_D1D2[ijk];
        g22[ijk] = psi4*gConf_D2D2[ijk];
        Rc_integrandx[ijk] = x-x_CM;
        Rc_integrandy[ijk] = y-y_CM;
        Rc_integrandz[ijk] = z-z_CM;
      }
    }
    DECLARE_FIELD(Rc_integrandx)
    DECLARE_FIELD(Rc_integrandy)
    DECLARE_FIELD(Rc_integrandz)
    Integration_T *I = init_integration();
    I->type = "Integral{f(x)dS},Spectral";
    I->g00 = g00;
    I->g01 = g01;
    I->g02 = g02;
    I->g11 = g11;
    I->g12 = g12;
    I->g22 = g22;
    I->Spectral->Z_surface = 1;
    I->Spectral->K         = 0;
    
    I->Spectral->f = Rc_integrandx;
    plan_integration(I);
    Rc[0] += execute_integration(I);
    
    I->Spectral->f = Rc_integrandy;
    plan_integration(I);
    Rc[1] += execute_integration(I);
    
    I->Spectral->f = Rc_integrandz;
    plan_integration(I);
    Rc[2] += execute_integration(I);
    
    free_integration(I);
    REMOVE_FIELD(Rc_integrandx)
    REMOVE_FIELD(Rc_integrandy)
    REMOVE_FIELD(Rc_integrandz)
    free(g00);
    free(g01);
    free(g02);
    free(g11);
    free(g12);
    free(g22);
  }
  Rc[0] /= (AH_area);
  Rc[1] /= (AH_area);
  Rc[2] /= (AH_area);
}
