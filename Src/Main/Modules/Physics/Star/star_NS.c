/*
// Alireza Rashti
// November 2020
*/


/* collection of NS star affairs */

#include "star_NS.h"

/* find Euler equation constant to meet NS baryonic mass */
int star_NS_idealfluid_gConf_find_Euler_const(Physics_T *const phys)
{
  FUNC_TIC
  
  Root_Finder_T *root = init_root_finder(1);
  const double W1  = Getd("Euler_const_update_weight");
  const double W2  = 1-W1;
  double *Euler_const = 0;
  double guess[1] = {Getd("Euler_equation_constant")};
  const double RESIDUAL = sqrt(Getd("RootFinder_Tolerance"));
  struct NS_Euler_eq_const_RootFinder_S params[1] = {0};
  double bar_mass,adm_mass,kommar_mass;
  
  bar_mass = star_NS_baryonic_gConf_mass(phys,guess[0]);
  observe(phys,"ADM(M)|NS",&adm_mass);
  observe(phys,"Kommar(M)|NS",&kommar_mass);

  printf(Pretty0"current NS baryonic mass = %e\n",bar_mass);
  printf(Pretty0"current NS ADM mass      = %e\n",adm_mass);
  printf(Pretty0"current NS Kommar mass   = %e\n",kommar_mass);
  
  Setd("baryonic_mass_current",bar_mass);
  Setd("ADM_mass",adm_mass);
  Setd("Kommar_mass",kommar_mass);
  
  params->phys = phys;
  params->NS_baryonic_mass = Getd("baryonic_mass");
  
  root->type        = Gets("RootFinder_Method");
  root->tolerance   = Getd("RootFinder_Tolerance");
  root->MaxIter     = (unsigned)Geti("RootFinder_Max_Number_of_Iteration");
  root->x_gss       = guess;
  root->params      = params;
  root->f[0]        = Euler_eq_const_gConf_rootfinder_eq;
  root->verbose     = 1;
  plan_root_finder(root);
  
  Euler_const       = execute_root_finder(root);
  /* if root finder is not OK for some reason */
  if (GRT(root->residual,RESIDUAL))
    Euler_const[0] = guess[0];/* don't update */
  else
    Euler_const[0] = W1*Euler_const[0]+W2*guess[0];
  
  Setd("Euler_equation_constant",Euler_const[0]);
  free(Euler_const);
  free_root_finder(root);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* root finder eqution for Euler equation constant */
static double Euler_eq_const_gConf_rootfinder_eq(void *params,const double *const x)
{
  struct NS_Euler_eq_const_RootFinder_S *const par = params; 
  
  return star_NS_baryonic_gConf_mass(par->phys,x[0]) - par->NS_baryonic_mass;
}

/* extrapolate matter fields out NS */
int star_NS_idealfluid_extrapolate_matter_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  if(Pcmps("star_extrapolate_matter_fields","poly2"))
  {
    /* make phi, W => enthalpy */
    const char *fields_name[] = {"phi","enthalpy",0};
    star_extrapolate(phys,fields_name,"poly2");
    phys->region = Ftype("NS_around");
    star_W_spin_vector_idealfluid_update(phys);
  }
  else
  {
    Error0(NO_OPTION);
  }
  
  FUNC_TOC
  return EXIT_SUCCESS;
}


/* calculate W = Omega_NS x (r-C_NS) */
void star_W_spin_vector_idealfluid_update(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid = phys->grid;
  double Omega_NS[3],C_NS[3];
  unsigned p;
  
  Omega_NS[0] = Getd("Omega_U0");
  Omega_NS[1] = Getd("Omega_U1");
  Omega_NS[2] = Getd("Omega_U2");
  C_NS[0]     = Getd("Center_x");
  C_NS[1]     = Getd("Center_y");
  C_NS[2]     = Getd("Center_z");
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if_not_cover(patch,phys) continue;
    
    W_spin_vector_idealfluid(patch,Omega_NS,C_NS);
    
  }
  
  FUNC_TOC
}

/* calculate W = Omega_NS x (r-C_NS) */
void W_spin_vector_idealfluid(Patch_T *const patch,const double Omega_NS[3],const double C_NS[3])
{
  REALLOC_v_WRITE_v(W_U0)
  REALLOC_v_WRITE_v(W_U1)
  REALLOC_v_WRITE_v(W_U2)
  unsigned nn = patch->nn;
  unsigned ijk;
  
  for (ijk = 0; ijk < nn; ++ijk)
  {
    double x = patch->node[ijk]->x[0]-C_NS[0];
    double y = patch->node[ijk]->x[1]-C_NS[1];
    double z = patch->node[ijk]->x[2]-C_NS[2];
    
    /* spin part */
    W_U0[ijk] = Omega_NS[1]*z-Omega_NS[2]*y;
    W_U1[ijk] = Omega_NS[2]*x-Omega_NS[0]*z;
    W_U2[ijk] = Omega_NS[0]*y-Omega_NS[1]*x;
  }
  
}

