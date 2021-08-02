/*
// Alireza Rashti
// November 2020
*/


/* functions for updates and computations pertinent to adm */


#include "adm_updates.h"

/* measuring Hamiltonian constraint.
// method:
// =======
// from_scratch: using beta,psi,alphaPsi and trK to make Kij, 
//               then using psi, gConf to make g
//               and then finally make R and using sources 
//               to calucalte constraints.
// from_identities: using AConfIJ and various other identities
//                  to calculate constraints.
// */
void adm_compute_constraints(Physics_T *const phys,
                                  const char *const region,
                                  const char *const method,
                                  const char *const ham,
                                  const char *const mom)
{
  if (strcmp_i(method,"from_identities"))
  {
    printf(Pretty0"method: from_identities.\n");
    fflush(stdout);
    
    Grid_T *const grid = mygrid(phys,region);
    Uint p;
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *const patch = grid->patch[p];
      
      adm_ham_and_mom_from_identities(patch,ham,mom);
      
    }
  }
  else if (strcmp_i(method,"from_scratch"))
  {
    Grid_T *const grid = mygrid(phys,region);
    Uint p;
    
    printf(Pretty0"method: from_scratch.\n");
    fflush(stdout);
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *const patch = grid->patch[p];
      
      adm_ham_and_mom_from_scratch(patch,ham,mom);
      
    }
  }
  else
    Error0(NO_OPTION);
  
}

/* update AConf^{ij} */
void adm_update_AConfIJ(Physics_T *const phys,const char *const region)
{
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *const patch = grid->patch[p];
    
    adm_update_AConfIJ_patch(patch);
    
  }
}

/* update adm_K^{ij} */
void adm_update_adm_KIJ(Physics_T *const phys,const char *const region)
{
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *const patch = grid->patch[p];
    
    adm_update_adm_KIJ_patch(patch);
    
  }
}

/* update adm_K_{ij} */
void adm_update_adm_Kij(Physics_T *const phys,const char *const region)
{
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *const patch = grid->patch[p];
    
    adm_update_adm_Kij_patch(patch);
    
  }
}

/* update adm_g_{ij} */
void adm_update_adm_gij(Physics_T *const phys,const char *const region)
{
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *const patch = grid->patch[p];
    
    adm_update_adm_g_patch(patch);
    
  }
}

/* handy macro for beta update */
#define UPDATE_beta(x)       \
    beta_U##x[ijk]       = B0_U##x[ijk]+B1_U##x[ijk];
    
#define UPDATE_dbeta(x,y)    \
   dbeta_U##x##D##y[ijk] = dB0_U##x##D##y[ijk]+dB1_U##x##D##y[ijk];

#define UPDATE_ddbeta(x,y,z) \
   ddbeta_U##x##D##y##D##z[ijk] = ddB0_U##x##D##y##D##z[ijk]+ddB1_U##x##D##y##D##z[ijk];

#define PREP_beta(x)       \
  REALLOC_v_WRITE_v(beta_U##x);READ_v(B0_U##x);READ_v(B1_U##x);
    
#define PREP_dbeta(x,y)    \
  REALLOC_v_WRITE_v(dbeta_U##x##D##y); READ_v(dB0_U##x##D##y);READ_v(dB1_U##x##D##y);
  
#define PREP_ddbeta(x,y,z) \
  REALLOC_v_WRITE_v(ddbeta_U##x##D##y##D##z); READ_v(ddB0_U##x##D##y##D##z);READ_v(ddB1_U##x##D##y##D##z);

/* update beta and its derivatives.
// note: it assumes B0 and B1 are already updated.
// note: the reasion we have a dedicated function for beta update
// is B1 part, since this might not get resolved well in outermost patches,
// thus, we should calculated it analytically. */
void adm_update_beta(Physics_T *const phys,const char *const region)
{
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *const patch = grid->patch[p];
    
    adm_update_beta_U0(patch);
    adm_update_beta_U1(patch);
    adm_update_beta_U2(patch);
  }
}

/* update beta_U0 patch wise. */
void adm_update_beta_U0(Patch_T *const patch)
{
  PREP_beta(0);
  PREP_dbeta(0,0);
  PREP_dbeta(0,1);
  PREP_dbeta(0,2);
  PREP_ddbeta(0,0,0);
  PREP_ddbeta(0,0,1);
  PREP_ddbeta(0,0,2);
  PREP_ddbeta(0,1,1);
  PREP_ddbeta(0,1,2);
  PREP_ddbeta(0,2,2);

  FOR_ALL_ijk
  {
    UPDATE_beta(0);
    UPDATE_dbeta(0,0);
    UPDATE_dbeta(0,1);
    UPDATE_dbeta(0,2);
    UPDATE_ddbeta(0,0,0);
    UPDATE_ddbeta(0,0,1);
    UPDATE_ddbeta(0,0,2);
    UPDATE_ddbeta(0,1,1);
    UPDATE_ddbeta(0,1,2);
    UPDATE_ddbeta(0,2,2);
  }
}

/* update beta_U1 patch wise. */
void adm_update_beta_U1(Patch_T *const patch)
{
  PREP_beta(1);
  PREP_dbeta(1,0);
  PREP_dbeta(1,1);
  PREP_dbeta(1,2);
  PREP_ddbeta(1,0,0);
  PREP_ddbeta(1,0,1);
  PREP_ddbeta(1,0,2);
  PREP_ddbeta(1,1,1);
  PREP_ddbeta(1,1,2);
  PREP_ddbeta(1,2,2);

  FOR_ALL_ijk
  {
    UPDATE_beta(1);
    UPDATE_dbeta(1,0);
    UPDATE_dbeta(1,1);
    UPDATE_dbeta(1,2);
    UPDATE_ddbeta(1,0,0);
    UPDATE_ddbeta(1,0,1);
    UPDATE_ddbeta(1,0,2);
    UPDATE_ddbeta(1,1,1);
    UPDATE_ddbeta(1,1,2);
    UPDATE_ddbeta(1,2,2);
  }
}

/* update beta_U2 patch wise. */
void adm_update_beta_U2(Patch_T *const patch)
{
  PREP_beta(2);
  PREP_dbeta(2,0);
  PREP_dbeta(2,1);
  PREP_dbeta(2,2);
  PREP_ddbeta(2,0,0);
  PREP_ddbeta(2,0,1);
  PREP_ddbeta(2,0,2);
  PREP_ddbeta(2,1,1);
  PREP_ddbeta(2,1,2);
  PREP_ddbeta(2,2,2);

  FOR_ALL_ijk
  {
    UPDATE_beta(2);
    UPDATE_dbeta(2,0);
    UPDATE_dbeta(2,1);
    UPDATE_dbeta(2,2);
    UPDATE_ddbeta(2,0,0);
    UPDATE_ddbeta(2,0,1);
    UPDATE_ddbeta(2,0,2);
    UPDATE_ddbeta(2,1,1);
    UPDATE_ddbeta(2,1,2);
    UPDATE_ddbeta(2,2,2);
  }
}

#undef UPDATE_beta
#undef UPDATE_dbeta
#undef UPDATE_ddbeta
#undef PREP_beta
#undef PREP_dbeta
#undef PREP_ddbeta

/* update B1^i and its derivatives */
void adm_update_adm_B1I(Physics_T *const phys,const char *const region)
{
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  if(Pcmps(P_"B1I_form","inspiral"))
  {
    struct General_Arg_S param[1] = {0};
    
    param->omega = sysGetd("angular_velocity");
    param->Vr    = sysGetd("infall_velocity");
    param->CM[0] = sysGetd("x_CM");
    param->CM[1] = sysGetd("y_CM");
    param->CM[2] = sysGetd("z_CM");
    param->D     = sysGetd("separation");
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *const patch = grid->patch[p];
      
      adm_update_B1I_patch(patch,param);
      
    }
  }
  else if(Pcmps(P_"B1I_form","zero"))
  {
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *const patch = grid->patch[p];
      
      adm_update_B1I_patch(patch,0);
      
    }
  }
  else
    Error0(NO_OPTION);
  
}

/* B1^i = omega x (r-r_CM) + Vr/D*(r-r_CM) */
void adm_update_B1I_inspiral(Patch_T *const patch,void *params)
{
  const struct General_Arg_S *const par = params;
  const double Omega = par->omega;
  const double Vr    = par->Vr;
  const double D     = par->D;
  const double x_CM  = par->CM[0];
  const double y_CM  = par->CM[1];
  const double z_CM  = par->CM[2];

  /* B^1 */
  REALLOC_v_WRITE_v(B1_U0)
  REALLOC_v_WRITE_v(B1_U1)
  REALLOC_v_WRITE_v(B1_U2)
  
  FOR_ALL_ijk
  {
    double x     = patch->node[ijk]->x[0];
    double y     = patch->node[ijk]->x[1];
    double z     = patch->node[ijk]->x[2];

    B1_U0[ijk] = Omega*(-y+y_CM)+Vr*(x-x_CM)/D;
    B1_U1[ijk] = Omega*(x-x_CM) +Vr*(y-y_CM)/D;
    B1_U2[ijk] =                 Vr*(z-z_CM)/D;
  }

  /* 1st derivatives */
  REALLOC_v_WRITE_v(dB1_U0D2)
  REALLOC_v_WRITE_v(dB1_U0D1)
  REALLOC_v_WRITE_v(dB1_U0D0)
  
  REALLOC_v_WRITE_v(dB1_U1D2)
  REALLOC_v_WRITE_v(dB1_U1D1)
  REALLOC_v_WRITE_v(dB1_U1D0)
  
  REALLOC_v_WRITE_v(dB1_U2D2)
  set_field_to_zero(dB1_U2D1)
  set_field_to_zero(dB1_U2D0)
  
  FOR_ALL_ijk
  {
    dB1_U0D2[ijk] = 0;
    dB1_U0D1[ijk] = -Omega;
    dB1_U0D0[ijk] = Vr/D;
    
    dB1_U1D2[ijk] = 0;
    dB1_U1D1[ijk] = Vr/D;
    dB1_U1D0[ijk] = Omega;
    
    /* dB1_U2D[01] is zero. */
    dB1_U2D2[ijk] = Vr/D;
  }

  /* 2nd derivatives => all zero. */
  set_field_to_zero(ddB1_U0D0D0)
  set_field_to_zero(ddB1_U0D1D2)
  set_field_to_zero(ddB1_U0D1D1)
  set_field_to_zero(ddB1_U0D0D1)
  set_field_to_zero(ddB1_U0D2D2)
  set_field_to_zero(ddB1_U0D0D2)
  
  set_field_to_zero(ddB1_U1D0D0)
  set_field_to_zero(ddB1_U1D1D2)
  set_field_to_zero(ddB1_U1D1D1)
  set_field_to_zero(ddB1_U1D0D1)
  set_field_to_zero(ddB1_U1D2D2)
  set_field_to_zero(ddB1_U1D0D2)
  
  set_field_to_zero(ddB1_U2D0D0)
  set_field_to_zero(ddB1_U2D1D2)
  set_field_to_zero(ddB1_U2D1D1)
  set_field_to_zero(ddB1_U2D0D1)
  set_field_to_zero(ddB1_U2D2D2)
  set_field_to_zero(ddB1_U2D0D2)
}

/* identically B1^i = 0 */
void adm_update_B1I_zero(Patch_T *const patch,void *params)
{
  set_field_to_zero(B1_U0)
  set_field_to_zero(B1_U1)
  set_field_to_zero(B1_U2)

  /* 1st derivatives */
  set_field_to_zero(dB1_U0D2)
  set_field_to_zero(dB1_U0D1)
  set_field_to_zero(dB1_U0D0)
  
  set_field_to_zero(dB1_U1D2)
  set_field_to_zero(dB1_U1D1)
  set_field_to_zero(dB1_U1D0)
  
  set_field_to_zero(dB1_U2D2)
  set_field_to_zero(dB1_U2D1)
  set_field_to_zero(dB1_U2D0)
  
  /* 2nd derivatives => all zero. */
  set_field_to_zero(ddB1_U0D0D0)
  set_field_to_zero(ddB1_U0D1D2)
  set_field_to_zero(ddB1_U0D1D1)
  set_field_to_zero(ddB1_U0D0D1)
  set_field_to_zero(ddB1_U0D2D2)
  set_field_to_zero(ddB1_U0D0D2)
  
  set_field_to_zero(ddB1_U1D0D0)
  set_field_to_zero(ddB1_U1D1D2)
  set_field_to_zero(ddB1_U1D1D1)
  set_field_to_zero(ddB1_U1D0D1)
  set_field_to_zero(ddB1_U1D2D2)
  set_field_to_zero(ddB1_U1D0D2)
  
  set_field_to_zero(ddB1_U2D0D0)
  set_field_to_zero(ddB1_U2D1D2)
  set_field_to_zero(ddB1_U2D1D1)
  set_field_to_zero(ddB1_U2D0D1)
  set_field_to_zero(ddB1_U2D2D2)
  set_field_to_zero(ddB1_U2D0D2) 
  
  UNUSED(params);
}

