/*
// Alireza Rashti
// November 2020
*/


/* collection of NS star affairs */

#include "star_NS.h"

/* find Euler equation constant to meet NS baryonic mass */
void star_idealfluid_NS_nonflat_find_Euler_const(Obj_Man_T *const obj)
{
  FUNC_TIC
  
  Grid_T *const grid = obj->grid;
  Root_Finder_T *root = init_root_finder(1);
  char opar[OPAR_LEN] = {'\0'};
  const double W1  = Getd("Euler_const_update_weight");
  const double W2  = 1-W1;
  double *Euler_const = 0;
  double guess[1] = {Getd("Euler_equation_constant")};
  const double RESIDUAL = sqrt(Getd("RootFinder_Tolerance"));
  struct NS_Euler_eq_const_RootFinder_S params[1] = {0};
  Observable_T *obs = 0;
  double bar_mass,adm_mass,kommar_mass;
  
  bar_mass = star_NS_baryonic_nonflat_mass(obj,guess[0]);
  obs = init_observable(grid,"ADM(M)|NS");
  adm_mass = obs->M(obs);
  free_observable(obs);
  
  obs = init_observable(grid,"Kommar(M)|NS");
  kommar_mass = obs->M(obs);
  free_observable(obs);

  printf(Pretty0"current NS baryonic mass = %e\n",bar_mass);
  printf(Pretty0"current NS ADM mass      = %e\n",adm_mass);
  printf(Pretty0"current NS Kommar mass   = %e\n",kommar_mass);
  
  Setd("baryonic_mass_current",bar_mass);
  Setd("ADM_mass",adm_mass);
  Setd("Kommar_mass",kommar_mass);
  
  params->obj = obj;
  params->NS_baryonic_mass = Getd("baryonic_mass");
  
  root->type        = Gets("RootFinder_Method");
  root->tolerance   = Getd("RootFinder_Tolerance");
  root->MaxIter     = (unsigned)Geti("RootFinder_Max_Number_of_Iteration");
  root->x_gss       = guess;
  root->params      = params;
  root->f[0]        = Euler_eq_const_nonflat_rootfinder_eq;
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
}

/* root finder eqution for Euler equation constant */
static double Euler_eq_const_nonflat_rootfinder_eq(void *params,const double *const x)
{
  struct NS_Euler_eq_const_RootFinder_S *const par = params; 
  
  return star_NS_baryonic_nonflat_mass(par->obj,x[0]) - par->NS_baryonic_mass;
}

