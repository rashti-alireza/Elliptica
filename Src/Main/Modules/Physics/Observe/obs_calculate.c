/*
// Alireza Rashti
// September 2019
// December  2020
*/

/* NOTE: please keep the format of implantation (with all of the ifs elses)
//       the same so later one can add new situations and cases easily. */

#include "obs_calculate.h"

/* plan and populate items_S sturct and obs struct
// for binary and single objects.
// algorithm: 
// ==========
//
// 1. collect all of the necessary patches.
// 2. populate the required metric for the integrations.
// 3. set flags for surface and volume integral.
// 4. populate normal vectors for surface integrals.
// 5. populate the integrands
// 6. call the pertinent functions for the calculation. */
void obs_calculate(Observe_T *const obs)
{
  IFss("ADM(P)|")
  {
    calc_ADM_P(obs);
  }
  else IFss("ADM(J)|")
  {
    calc_ADM_J(obs);
  }
  else IFss("ADM(M)|")
  {
    calc_ADM_mass(obs);
  }
  else IFss("Komar(M)|")
  {
    calc_Kommar_mass(obs);
  }
  else IFss("Irreducible(M)|")
  {
    calc_irreducible_BH_mass(obs);
  }
  else IFss("Baryonic(M)|")
  {
    calc_baryonic_mass(obs);
  }
  else IFss("CM|")
  {
    calc_CM(obs);
  }
  else IFss("Spin|")
  {
    calc_spin(obs);
  }
  else
    Errors("There is no such '%s' plan.\n",obs->quantity);
    
}

/* populating normal outward vector for around patches according to the given dir 
// NOTE: the normalization is respect to the physical metric gamma_{ij} */
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
// NOTE: the normalization is respect to the conformal metric gamma_{ij} */
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
    Free(normal->n_U0);
    Free(normal->n_U1);
    Free(normal->n_U2);
    Free(g00);
    Free(g01);
    Free(g02);
    Free(g11);
    Free(g12);
    Free(g22);
  }
  
  Free(patches);
  
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
    Free(normal->n_U0);
    Free(normal->n_U1);
    Free(normal->n_U2);
    Free(g00);
    Free(g01);
    Free(g02);
    Free(g11);
    Free(g12);
    Free(g22);
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

  R[0] = obs->ret[0];
  R[1] = obs->ret[1];
  R[2] = obs->ret[2];

  S[0] = J[0] - (-P[1]*R[2] + P[2]*R[1]);
  S[1] = J[1] - (P[0]*R[2] - P[2]*R[0]);
  S[2] = J[2] - (-P[0]*R[1] + P[1]*R[0]);
}


/* calculating physical center of BH to be used in spin calculations */
static void Rc_BH(Observe_T *const obs)
{
  Physics_T *const phys = obs->phys;
  Grid_T *const grid    = mygrid(phys,"BH_around_IB");
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


/* calculate adm J for various objects 
// NOTE: the accuracy of this method required to take the integral in
// outermost patches where the metric is conformally flat and trK = 0. */
static void calc_ADM_J(Observe_T *const obs)
{
  SET_MSG
  
  Grid_T *const grid   = obs->grid;  
  Patch_T **patches1   = 0;/* for volume integrals */
  Patch_T **patches2   = 0;/* for surface integrals */
  Patch_T *patch       = 0;
  const char *region   = 0;
  struct items_S **adm = 0;
  Uint N1 = 0;
  Uint N2 = 0;
  Uint n,ijk,nn;
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_SBH)
  {
    IFsc("ADM(J)|BHNS")
    {
      /* NOTE: the "S_inf,default" and "S+V,Ossokine" methods predict the value of system's J 
      // differently from the "S+V,constraint" and "S_obj1+S_obj2,default" methods.
      // the default using the latter for J_ADM measure, and I believed they're the correct ones.
      // however, further investigation in required to find the root cause of this
      // difference. I should emphasize that this only affect the value of the system's J 
      // at the post-processing step and has nothing to do with spins, P adms, or 
      // the ID quality. ALL other quantities are measured accurately. */
      if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S+V,Ossokine"))
      {
        /* this method does not have volume integrals */
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S+V,constraint"))
      {
        /* volume part */
        region   = "outermost,filling_box,NS_around,BH_around";
        patches1 = collect_patches(grid,region,&N1);
        
        /* surface part */
        region   = "NS_around_IB,BH_around_IB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S_obj1+S_obj2,default"))
      {
        /* surface part */
        region   = "NS_OB,BH_around_IB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(J)|NSNS")
    {
      /* NOTE: the "S_inf,default" and "S+V,Ossokine" methods predict the value of system's J 
      // differently from the "S+V,constraint" and "S_obj1+S_obj2,default" methods.
      // the default using the latter for J_ADM measure, and I believed they're the correct ones.
      // however, further investigation in required to find the root cause of this
      // difference. I should emphasize that this only affect the value of the system's J 
      // at the post-processing step and has nothing to do with spins, P adms, or 
      // the ID quality. ALL other quantities are measured accurately. */

      if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S+V,Ossokine"))
      {
        /* this method does not have volume integrals */
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S+V,constraint"))
      {
        /* volume part */
        region   = "outermost,filling_box,NS1_around,NS2_around";
        patches1 = collect_patches(grid,region,&N1);
        
        /* surface part */
        region   = "NS1_around_IB,NS2_around_IB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S_obj1+S_obj2,default"))
      {
        /* NOTE:for maximal slice and conformal flat metric 
        // volume integral is 0. like the case we have for NSNS */
        /* surface part */
        region   = "NS1_OB,NS2_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(J)|NS")
    {
      if (IsIt("S_obj,default"))
      {
        /* surface part */
        region   = "NS_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(J)|NS1")
    {
      if (IsIt("S_obj,default"))
      {
        /* surface part */
        region   = "NS1_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(J)|NS2")
    {
      if (IsIt("S_obj,default"))
      {
        /* surface part */
        region   = "NS2_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(J)|BH")
    {
      if (IsIt("S_obj,default"))
      {
        /* surface part */
        region   = "BH_around_IB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(J)|SBH")
    {
      if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S+V,Ossokine"))
      {
        /* this method does not have volume integrals */
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
  
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
    
    if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_NSNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH)
    {
      IFsc("ADM(J)|BHNS")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(J)|NSNS")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(J)|NS")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(J)|NS1")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(J)|NS2")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(J)|BH")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(J)|SBH")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else
    {
      Error0(NO_OPTION);
    }
  }
  Free(patches1);
 
  /* fill ADM struct for each patch surface part */
  for (n = N1; n < N1+N2; ++n)
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
    
    if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_NSNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH)
    {
      IFsc("ADM(J)|BHNS")
      {
        if (IsIt("S_inf,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,Ossokine"))
        {
          Set_outermost_integral_S_SplitCS(adm)
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,constraint"))
        {
          if (IsItCovering(patch,"NS_around_IB,BH_around_IB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = 0;
            n_physical_metric_around(adm[n],_c_);
          }
          else
          {
            Error0(obs_err_msg);
          }
        }
        else if (IsIt("S_obj1+S_obj2,default"))
        {
          if (IsItCovering(patch,"BH_around_IB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = 0;
            n_physical_metric_around(adm[n],_c_);
          }
          else if (IsItCovering(patch,"NS_OB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = patch->n[2]-1;
            n_physical_metric_around(adm[n],_c_);
          }
          else
          {
            Error0(obs_err_msg);
          }
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(J)|NSNS")
      {
        if (IsIt("S_inf,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,Ossokine"))
        {
          Set_outermost_integral_S_SplitCS(adm)
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,constraint"))
        {
          if (IsItCovering(patch,"NS1_around_IB,NS2_around_IB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = 0;
            n_physical_metric_around(adm[n],_c_);
          }
          else
          {
            Error0(obs_err_msg);
          }
        }
        else if (IsIt("S_obj1+S_obj2,default"))
        {
          if (IsItCovering(patch,"NS1_OB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = patch->n[2]-1;
            n_physical_metric_around(adm[n],_c_);
          }
          else if (IsItCovering(patch,"NS2_OB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = patch->n[2]-1;
            n_physical_metric_around(adm[n],_c_);
          }
          else
          {
            Error0(obs_err_msg);
          }
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(J)|NS")
      {
        if (IsIt("S_obj,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(J)|NS1")
      {
        if (IsIt("S_obj,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(J)|NS2")
      {
        if (IsIt("S_obj,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(J)|BH")
      {
        if (IsIt("S_obj,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = 0;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(J)|SBH")
      {
        if (IsIt("S_inf,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,Ossokine"))
        {
          Set_outermost_integral_S_SplitCS(adm)
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else
    {
      Error0(NO_OPTION);
    }
  }
  Free(patches2);


  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_SBH)
  {
    if (IsIt("S_inf,default"))
    {
      obs_ADM_J_S_default(obs);
    }
    else if (IsIt("S_obj,default"))
    {
      obs_ADM_J_S_default(obs);
    }
    else if (IsIt("S_obj1+S_obj2,default"))
    {
      obs_ADM_J_S_default(obs);
    }
    else if (IsIt("S+V,Ossokine"))
    {
      obs_ADM_J_Stokes_SV_Ossokine(obs);
    }
    else if (IsIt("S+V,constraint"))
    {
      obs_ADM_J_Stokes_SV_constraint(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
}

/* calculate adm P for various objects 
// NOTE: the accuracy of this method required to take the integral in
// outermost patches where the metric is conformally flat and trK = 0. */
static void calc_ADM_P(Observe_T *const obs)
{
  SET_MSG
  
  Grid_T *const grid   = obs->grid;  
  Patch_T **patches1   = 0;/* for volume integrals */
  Patch_T **patches2   = 0;/* for surface integrals */
  Patch_T *patch       = 0;
  const char *region   = 0;
  struct items_S **adm = 0;
  Uint N1 = 0;
  Uint N2 = 0;
  Uint n,ijk,nn;
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_SBH)
  {
    IFsc("ADM(P)|BHNS")
    {
      if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S+V,Ossokine"))
      {
        Collect_outermost_S_V_SplitCS
      }
      else if (IsIt("S+V,Rashti"))
      {
        /* volume part */
        region   = "outermost,filling_box,NS_around,NS,BH_around";
        patches1 = collect_patches(grid,region,&N1);
        
        /* surface part */
        region   = "BH_around_IB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S+V,constraint"))
      {
        /* volume part */
        region   = "outermost,filling_box,NS_around,BH_around";
        patches1 = collect_patches(grid,region,&N1);
        
        /* surface part */
        region   = "NS_around_IB,BH_around_IB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S_obj1+S_obj2,default"))
      {
        /* surface part */
        region   = "NS_OB,BH_around_IB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(P)|NSNS")
    {
      if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S+V,Ossokine"))
      {
        Collect_outermost_S_V_SplitCS
      }
      else if (IsIt("S+V,Rashti"))
      {
        /* volume part */
        region   = "outermost,filling_box,NS1_around,NS1,NS2_around,NS2";
        patches1 = collect_patches(grid,region,&N1);
      }
      else if (IsIt("S+V,constraint"))
      {
        /* volume part */
        region   = "outermost,filling_box,NS1_around,NS1,NS2_around,NS2";
        patches1 = collect_patches(grid,region,&N1);
        
        /* surface part */
        region   = "NS1_around_IB,NS2_around_IB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S_obj1+S_obj2,default"))
      {
        /* surface part */
        region   = "NS1_OB,NS2_OB,";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(P)|NS")
    {
      if (IsIt("S_obj,default"))
      {
        /* surface part */
        region   = "NS_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(P)|NS1")
    {
      if (IsIt("S_obj,default"))
      {
        /* surface part */
        region   = "NS1_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(P)|NS2")
    {
      if (IsIt("S_obj,default"))
      {
        /* surface part */
        region   = "NS2_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(P)|BH")
    {
      if (IsIt("S_obj,default"))
      {
        /* surface part */
        region   = "BH_around_IB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(P)|SBH")
    {
      if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else if (IsIt("S+V,Ossokine"))
      {
        Collect_outermost_S_V_SplitCS
      }
      else if (IsIt("S+V,Rashti"))
      {
        /* volume part */
        region   = "outermost,BH_around";
        patches1 = collect_patches(grid,region,&N1);
        
        /* surface part */
        region   = "BH_around_IB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
  
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
    
    if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_NSNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH)
    {
      IFsc("ADM(P)|BHNS")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(P)|NSNS")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(P)|NS")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(P)|NS1")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(P)|NS2")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(P)|BH")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else IFsc("ADM(P)|SBH")
      {
        ;/* nothing, to keep all options this empty if stays here */
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else
    {
      Error0(NO_OPTION);
    }
  }
  Free(patches1);
 
  /* fill ADM struct for each patch surface part */
  for (n = N1; n < N1+N2; ++n)
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
    
    if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_NSNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH)
    {
      IFsc("ADM(P)|BHNS")
      {
        if (IsIt("S_inf,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,Ossokine"))
        {
          Set_outermost_integral_S_SplitCS(adm)
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,Rashti"))
        {
          if (IsItCovering(patch,"BH_around_IB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = 0;
            n_physical_metric_around(adm[n],_c_);
          }
          else
          {
            Error0(obs_err_msg);
          }
        }
        else if (IsIt("S+V,constraint"))
        {
          if (IsItCovering(patch,"NS_around_IB,BH_around_IB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = 0;
            n_physical_metric_around(adm[n],_c_);
          }
          else
          {
            Error0(obs_err_msg);
          }
        }
        else if (IsIt("S_obj1+S_obj2,default"))
        {
          if (IsItCovering(patch,"BH_around_IB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = 0;
            n_physical_metric_around(adm[n],_c_);
          }
          else if (IsItCovering(patch,"NS_OB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = patch->n[2]-1;
            n_physical_metric_around(adm[n],_c_);
          }
          else
          {
            Error0(obs_err_msg);
          }
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(P)|NSNS")
      {
        if (IsIt("S_inf,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,Ossokine"))
        {
          Set_outermost_integral_S_SplitCS(adm)
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,Rashti"))
        {
          ;
        }
        else if (IsIt("S+V,constraint"))
        {
          if (IsItCovering(patch,"NS1_around_IB,NS2_around_IB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = 0;
            n_physical_metric_around(adm[n],_c_);
          }
          else
          {
            Error0(obs_err_msg);
          }
        }
        else if (IsIt("S_obj1+S_obj2,default"))
        {
          if (IsItCovering(patch,"NS1_OB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = patch->n[2]-1;
            n_physical_metric_around(adm[n],_c_);
          }
          else if (IsItCovering(patch,"NS2_OB"))
          {
            adm[n]->surface_integration_flg = 1;
            adm[n]->Z_surface = 1;
            adm[n]->K = patch->n[2]-1;
            n_physical_metric_around(adm[n],_c_);
          }
          else
          {
            Error0(obs_err_msg);
          }
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(P)|NS")
      {
        if (IsIt("S_obj,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(P)|NS1")
      {
        if (IsIt("S_obj,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(P)|NS2")
      {
        if (IsIt("S_obj,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(P)|BH")
      {
        if (IsIt("S_obj,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = 0;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(P)|SBH")
      {
        if (IsIt("S_inf,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,Ossokine"))
        {
          Set_outermost_integral_S_SplitCS(adm)
          n_physical_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,Rashti"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = 0;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else
    {
      Error0(NO_OPTION);
    }  
  }
  Free(patches2);
  
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_SBH)
  {
    if (IsIt("S_inf,default"))
    {
      obs_ADM_P_S_default(obs);
    }
    else if (IsIt("S_obj,default"))
    {
      obs_ADM_P_S_default(obs);
    }
    else if (IsIt("S_obj1+S_obj2,default"))
    {
      obs_ADM_P_S_default(obs);
    }
    else if (IsIt("S+V,Ossokine"))
    {
      obs_ADM_P_Stokes_SV_Ossokine(obs);
    }
    else if (IsIt("S+V,Rashti"))
    {
      obs_ADM_P_Stokes_SV_Rashti(obs);
    }
    else if (IsIt("S+V,constraint"))
    {
      obs_ADM_P_Stokes_SV_constraint(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
  
}
  
/* calculate Komar mass for various objects */
static void calc_Kommar_mass(Observe_T *const obs)
{
  SET_MSG
  
  Grid_T *const grid    = obs->grid;
  Patch_T **patches1   = 0;/* for volume integrals */
  Patch_T **patches2   = 0;/* for surface integrals */
  Patch_T *patch       = 0;
  const char *region   = 0;
  struct items_S **Komar = 0;
  Uint N1 = 0;
  Uint N2 = 0;
  Uint n,ijk,nn;
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_SBH)
  {
    IFsc("Komar(M)|BHNS")
    {
      if (IsIt("S+V,default"))
      {
        /* volume part */
        region   = "NS";
        patches1 = collect_patches(grid,region,&N1);
        /* surface part */
        region   = "BH_around_IB";
        patches2 = collect_patches(grid,region,&N2); 
      }
      else if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);  
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("Komar(M)|NSNS")
    {
      if (IsIt("S+V,default"))
      {
        /* volume part */
        region   = "NS1,NS2";
        patches1 = collect_patches(grid,region,&N1);
      }
      else if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);  
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("Komar(M)|NS")
    {
      if (IsIt("V_obj,default"))
      {
        /* volume part */
        region = "NS";
        patches1 = collect_patches(grid,region,&N1);
      }
      else if (IsIt("S_obj,default"))
      {
        /* surface part */
        region = "NS_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("Komar(M)|NS1")
    {
      if (IsIt("V_obj,default"))
      {
        /* volume part */
        region = "NS1";
        patches1 = collect_patches(grid,region,&N1);
      }
      else if (IsIt("S_obj,default"))
      {
        /* surface part */
        region = "NS1_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("Komar(M)|NS2")
    {
      if (IsIt("V_obj,default"))
      {
        /* volume part */
        region = "NS2";
        patches1 = collect_patches(grid,region,&N1);
      }
      else if (IsIt("S_obj,default"))
      {
        /* surface part */
        region = "NS2_OB";
        patches2 = collect_patches(grid,region,&N2);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("Komar(M)|BH")
    {
      /* surface part */
      region   = "BH_around_IB";
      patches2 = collect_patches(grid,region,&N2); 
    }
    else IFsc("Komar(M)|SBH")
    {
      if (IsIt("S_obj,default"))
      {
        /* surface part */
        region   = "BH_around_IB";
        patches2 = collect_patches(grid,region,&N2); 
      }
      else if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2); 
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
  
  /* alloc memory for all patches */
  Komar = calloc(N1+N2,sizeof(*Komar));
  IsNull(Komar);
  /* this is where we link to obs struct */
  obs->items = Komar;
  obs->Nitems = N1+N2;
      
  /* fill Komar struct for each patch volume part */
  for (n = 0; n < N1; ++n)
  {
    Komar[n] = calloc(1,sizeof(*Komar[n]));
    IsNull(Komar[n]);
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
    READ_v(psi);
    
    Komar[n]->patch = patch;
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
    Komar[n]->g00 = g00;
    Komar[n]->g01 = g01;
    Komar[n]->g02 = g02;
    Komar[n]->g11 = g11;
    Komar[n]->g12 = g12;
    Komar[n]->g22 = g22;
  }
  Free(patches1);
  
  /* fill Komar struct for each patch surface part */
  for (n = N1; n < N1+N2; ++n)
  {
    Komar[n] = calloc(1,sizeof(*Komar[n]));
    IsNull(Komar[n]);
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
    READ_v(psi);
    
    Komar[n]->patch = patch;
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
    Komar[n]->g00 = g00;
    Komar[n]->g01 = g01;
    Komar[n]->g02 = g02;
    Komar[n]->g11 = g11;
    Komar[n]->g12 = g12;
    Komar[n]->g22 = g22;
    
    if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_NSNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH)
    {
      IFsc("Komar(M)|BHNS")
      {
        if (IsIt("S+V,default"))
        {
          /* surface integral */
          Komar[n]->surface_integration_flg = 1;
          Komar[n]->Z_surface = 1;
          Komar[n]->K = 0;
          n_physical_metric_around(Komar[n],_c_);
        }
        else if (IsIt("S_inf,default"))
        {
          /* NOTE: we can use a closer surface to the objects
          // since Komar is independent of surface, so: */
          Set_outermost_integral_S_SplitCS(Komar)
          n_physical_metric_around(Komar[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("Komar(M)|NSNS")
      {
        if (IsIt("S+V,default"))
        {
          ;
        }
        else if (IsIt("S_inf,default"))
        {
          /* NOTE: we can use a closer surface to the objects
          // since Komar is independent of surface, so: */
          Set_outermost_integral_S_SplitCS(Komar)
          n_physical_metric_around(Komar[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("Komar(M)|NS")
      {
        if (IsIt("V_obj,default"))
        {
          ;
        }
        else if (IsIt("S_obj,default"))
        {
          /* surface integral */
          Komar[n]->surface_integration_flg = 1;
          Komar[n]->Z_surface = 1;
          Komar[n]->K = patch->n[2]-1;
          n_physical_metric_around(Komar[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("Komar(M)|NS1")
      {
        if (IsIt("V_obj,default"))
        {
          ;
        }
        else if (IsIt("S_obj,default"))
        {
          /* surface integral */
          Komar[n]->surface_integration_flg = 1;
          Komar[n]->Z_surface = 1;
          Komar[n]->K = patch->n[2]-1;
          n_physical_metric_around(Komar[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("Komar(M)|NS2")
      {
        if (IsIt("V_obj,default"))
        {
          ;
        }
        else if (IsIt("S_obj,default"))
        {
          /* surface integral */
          Komar[n]->surface_integration_flg = 1;
          Komar[n]->Z_surface = 1;
          Komar[n]->K = patch->n[2]-1;
          n_physical_metric_around(Komar[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("Komar(M)|BH")
      {
        /* surface integral */
        Komar[n]->surface_integration_flg = 1;
        Komar[n]->Z_surface = 1;
        Komar[n]->K = 0;
        n_physical_metric_around(Komar[n],_c_);
      }
      else IFsc("Komar(M)|SBH")
      {
        if (IsIt("S_obj,default"))
        {
          /* surface integral */
          Komar[n]->surface_integration_flg = 1;
          Komar[n]->Z_surface = 1;
          Komar[n]->K = 0;
          n_physical_metric_around(Komar[n],_c_);
        }
        else if (IsIt("S_inf,default"))
        {
          /* NOTE: we can use a closer surface to the objects
          // since Komar is independent of surface, so: */
          Set_outermost_integral_S_SplitCS(Komar)
          n_physical_metric_around(Komar[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else
    {
      Error0(NO_OPTION);
    }
  }
  Free(patches2);
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_SBH)
  {
    if (IsIt("S+V,default"))
    {
      obs->ret[0] = obs_Komar_mass(obs);
    }
    else if (IsIt("S_obj,default"))
    {
      obs->ret[0] = obs_Komar_mass(obs);
    }
    else if (IsIt("S_inf,default"))
    {
      obs->ret[0] = obs_Komar_mass(obs);
    }
    else if (IsIt("V_obj,default"))
    {
      obs->ret[0] = obs_Komar_mass(obs);
    }
    else if (IsIt("S+V,default"))
    {
      obs->ret[0] = obs_Komar_mass(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
}

/* calculate ADM mass for various objects */
static void calc_ADM_mass(Observe_T *const obs)
{
  SET_MSG
  
  /* in these cases use gConf and not g */
  const int IsConf   = (
                        IsIt("S+V,default")   ||
                        IsIt("V_obj,default") ||
                        IsIt("S+V,conformal")
                       );
  Grid_T *const grid = obs->grid;
  Patch_T **patches1 = 0;/* for volume integrals */
  Patch_T **patches2 = 0;/* for surface integrals */
  Patch_T *patch     = 0;
  struct items_S **adm = 0;
  const char *region   = 0;
  Uint N1 = 0;
  Uint N2 = 0;
  Uint n,ijk,nn;
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_SBH)
  {
    IFsc("ADM(M)|BHNS")
    {
      if (IsIt("S+V,default"))
      {
        /* volume part */
        region   = "outermost,filling_box,NS,NS_around,BH_around";
        patches1 = collect_patches(grid,region,&N1);
        /* surface part */
        region   = "BH_around_IB";
        patches2 = collect_patches(grid,region,&N2); 
      }
      else if (IsIt("S+V,conformal"))
      {
        /* volume part */
        region   = "outermost,filling_box,NS,NS_around,BH_around";
        patches1 = collect_patches(grid,region,&N1);
        /* surface part */
        region   = "BH_around_IB";
        patches2 = collect_patches(grid,region,&N2); 
      }
      else if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);  
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(M)|NSNS")
    {
      if (IsIt("S+V,default"))
      {
        /* volume part */
        region   = "outermost,filling_box,NS1,NS1_around,NS2,NS2_around";
        patches1 = collect_patches(grid,region,&N1);
      }
      else if (IsIt("S+V,conformal"))
      {
        /* volume part */
        region   = "outermost,filling_box,NS1,NS1_around,NS2,NS2_around";
        patches1 = collect_patches(grid,region,&N1);
      }
      else if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);  
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(M)|NS")
    {
      if (IsIt("V_obj,default"))
      {
        region   = "NS";
        patches1 = collect_patches(grid,region,&N1);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(M)|NS1")
    {
      if (IsIt("V_obj,default"))
      {
        region   = "NS1";
        patches1 = collect_patches(grid,region,&N1);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(M)|NS2")
    {
      if (IsIt("V_obj,default"))
      {
        region   = "NS2";
        patches1 = collect_patches(grid,region,&N1);
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else IFsc("ADM(M)|SBH")
    {
      if (IsIt("S+V,default"))
      {
        /* volume part */
        region   = "outermost,BH_around";
        patches1 = collect_patches(grid,region,&N1);
        /* surface part */
        region   = "BH_around_IB";
        patches2 = collect_patches(grid,region,&N2); 
      }
      else if (IsIt("S_inf,default"))
      {
        /* surface part */
        region   = "outermost_OB";
        patches2 = collect_patches(grid,region,&N2);  
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(NO_OPTION);
  }

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
    READ_v(psi);
    
    adm[n]->patch = patch;
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double psi4 = (IsConf ? 1. : Pow2(psi[ijk])*Pow2(psi[ijk]));
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
  }
  Free(patches1);
  
  /* for surface part */
  for (n = N1; n < N1+N2; ++n)
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
    READ_v(psi);
    
    adm[n]->patch = patch;
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double psi4 = (IsConf ? 1. : Pow2(psi[ijk])*Pow2(psi[ijk]));
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
    
    if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_NSNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH)
    {
      IFsc("ADM(M)|BHNS")
      {
        if (IsIt("S+V,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = 0;
          n_conformal_metric_around(adm[n],_c_);
        }
        else if (IsIt("S+V,conformal"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = 0;
          n_conformal_metric_around(adm[n],_c_);
        }
        else if (IsIt("S_inf,default"))
        {
          /* surface integral */
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(M)|NSNS")
      {
        if (IsIt("S+V,default"))
        {
          ;
        }
        else if (IsIt("S+V,conformal"))
        {
          ;
        }
        else if (IsIt("S_inf,default"))
        {
          /* surface integral */
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(M)|NS")
      {
        if (IsIt("V_obj,default"))
        {
          ;
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(M)|NS1")
      {
        if (IsIt("V_obj,default"))
        {
          ;
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(M)|NS2")
      {
        if (IsIt("V_obj,default"))
        {
          ;
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else IFsc("ADM(M)|SBH")
      {
        if (IsIt("S+V,default"))
        {
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = 0;
          n_conformal_metric_around(adm[n],_c_);
        }
        else if (IsIt("S_inf,default"))
        {
          /* surface integral */
          adm[n]->surface_integration_flg = 1;
          adm[n]->Z_surface = 1;
          adm[n]->K = patch->n[2]-1;
          n_physical_metric_around(adm[n],_c_);
        }
        else
        {
          Error0(obs_err_msg);
        }
      }
      else
      {
        Error0(obs_err_msg);
      }
    }
    else
      Error0(NO_OPTION);
  }
  Free(patches2);
  
  if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
      grid->kind == Grid_SplitCubedSpherical_NSNS ||
      grid->kind == Grid_SplitCubedSpherical_SBH)
  {
    if (IsIt("S+V,default"))
    {
      obs->ret[0] = obs_ADM_mass_SV_isotropic(obs);
    }
    else if (IsIt("S+V,conformal"))
    {
      obs->ret[0] = obs_ADM_mass_SV_conformal(obs);
    }
    else if (IsIt("S_inf,default"))
    {
      obs->ret[0] = obs_ADM_mass_S2(obs);
    }
    else if (IsIt("V_obj,default"))
    {
      obs->ret[0] = obs_ADM_mass_SV_isotropic(obs);
    }
    else 
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
}

/* ->: (+/-)integral{S ds} + (+/-)integral{V dv}.
// Carring out integration of ADM P or J. */
double obs_integral_SV (Observe_T *const obs,
                        const char *const sS/* integrand for S */,
                        const char *const sV/* intergrand for V */,
                        const char sign_sS/* [+/-] integral of S */,
                        const char sign_sV/* [+/-] integral of V */)
{
  double ret = 0;
  struct items_S **const adm = obs->items;
  const Uint N = obs->Nitems;
  const double Sign[2] = {-1.,1.};
  Uint p;
  
  /* some checks */
  assert(N);
  assert(sign_sS == '+' || sign_sS == '-');
  assert(sign_sV == '+' || sign_sV == '-');
  
  for(p = 0; p < N; ++p)
  {
    Patch_T *patch = adm[p]->patch;
    double sum_s = 0.;
    double sum_v = 0.;
    
    if (adm[p]->surface_integration_flg)
    {
      Integration_T *I = 0;
      I  = init_integration();
      I->type = "Integral{f(x)dS},Spectral";
      I->Spectral->f = patch->fields[Ind(sS)];
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
      sum_s = execute_integration(I);
      free_integration(I);
    }
    else
    {
      Integration_T *I = 0;
      I  = init_integration();
      I->type = "Integral{f(x)dV},Spectral";
      I->Spectral->f = patch->fields[Ind(sV)];
      I->g00 = adm[p]->g00;
      I->g01 = adm[p]->g01;
      I->g02 = adm[p]->g02;
      I->g11 = adm[p]->g11;
      I->g12 = adm[p]->g12;
      I->g22 = adm[p]->g22;
      
      plan_integration(I);
      sum_v = execute_integration(I);
      free_integration(I);
    }
    ret += (Sign[sign_sS == '+']*sum_s+Sign[sign_sV == '+']*sum_v);
  }
  
  return ret;
}


/* calculate irreducible mass of BH */
static void calc_irreducible_BH_mass(Observe_T *const obs)
{
  SET_MSG
  
  Grid_T *const grid = obs->grid;
  
  if (IsIt("S_obj,default"))
  {
    if (grid->kind == Grid_SplitCubedSpherical_BHNS ||
        grid->kind == Grid_SplitCubedSpherical_SBH)
    {
      obs_BH_irreducible_mass_CS(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(obs_err_msg);
  }

}

/* center mass 
// NOTE: CM of each object measured with respect to system CM. */
static void calc_CM(Observe_T *const obs)
{
  SET_MSG
  
  IFsc("CM|BH")
  {
    if (IsIt("S_obj,default"))
    {
      Rc_BH(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else IFsc("CM|NS")
  {
    if (IsIt("V_obj,default"))
    {
      obs_Rc_NS(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else IFsc("CM|NS1")
  {
    if (IsIt("V_obj,default"))
    {
      obs_Rc_NS(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else IFsc("CM|NS2")
  {
    if (IsIt("V_obj,default"))
    {
      obs_Rc_NS(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(obs_err_msg);
  }
}

/* calculate spin */
static void calc_spin(Observe_T *const obs)
{
  SET_MSG
  
  IFsc("Spin|BH")
  {
    if (IsIt("S_obj,JRP"))
    {
      define_spin_JRP(obs);
    }
    else if (IsIt("S_obj,Campanelli"))
    {
      define_spin_campanelli(obs);
    }
    else if (IsIt("S_obj,AKV"))
    {
      define_spin_akv(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else IFsc("Spin|NS")
  {
    if (IsIt("S_obj,JRP"))
    {
      define_spin_JRP(obs);
    }
    else if (IsIt("S_obj,Campanelli"))
    {
      define_spin_campanelli(obs);
    }
    else if (IsIt("S_obj,AKV"))
    {
      define_spin_akv(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else IFsc("Spin|NS1")
  {
    if (IsIt("S_obj,JRP"))
    {
      define_spin_JRP(obs);
    }
    else if (IsIt("S_obj,Campanelli"))
    {
      define_spin_campanelli(obs);
    }
    else if (IsIt("S_obj,AKV"))
    {
      define_spin_akv(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else IFsc("Spin|NS2")
  {
    if (IsIt("S_obj,JRP"))
    {
      define_spin_JRP(obs);
    }
    else if (IsIt("S_obj,Campanelli"))
    {
      define_spin_campanelli(obs);
    }
    else if (IsIt("S_obj,AKV"))
    {
      define_spin_akv(obs);
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(obs_err_msg);
  }
}


/* calculate Baryonic mass */
static void calc_baryonic_mass(Observe_T *const obs)
{
  SET_MSG
  
  Physics_T *const phys = obs->phys;
  
  IFsc("Baryonic(M)|NS")
  {
    if (IsIt("V_obj,default"))
    {
      obs->ret[0] = star_NS_baryonic_gConf_mass
          (phys,Getd("Euler_equation_constant"));
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else IFsc("Baryonic(M)|NS1")
  {
    if (IsIt("V_obj,default"))
    {
      obs->ret[0] = star_NS_baryonic_gConf_mass
          (phys,Getd("Euler_equation_constant"));
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else IFsc("Baryonic(M)|NS2")
  {
    if (IsIt("V_obj,default"))
    {
      obs->ret[0] = star_NS_baryonic_gConf_mass
          (phys,Getd("Euler_equation_constant"));
    }
    else
    {
      Error0(obs_err_msg);
    }
  }
  else
  {
    Error0(obs_err_msg);
  }
}
