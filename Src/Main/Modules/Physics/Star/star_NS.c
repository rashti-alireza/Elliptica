/*
// Alireza Rashti
// November 2020
*/


/* collection of NS star affairs */

#include "star_NS.h"

/* adjust force balance */
int star_NS_idealfluid_gConf_force_balance(Physics_T *const phys)
{
  FUNC_TIC
  
  Patch_T *patch    = 0;
  Patch_T **patches = 0;
  const double NS_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  const char *const par = Gets("force_balance_equation");
  Interpolation_T *interp_s = init_interpolation();
  char *adjust[3];
  double dh1[3] = {0},dh2[3] = {0}, X[3] = {0};
  unsigned Np;
  char reg[99];
  
  sprintf(reg,"%s_cent.*",phys->spos);
  patches = regex_collect_patches(phys->grid,reg,&Np);
  
  /* initial values before adjustments */
  patch = x_in_which_patch(NS_center,patches,Np);
  assert(patch);
  assert(X_of_x(X,NS_center,patch));
  
  DECLARE_FIELD(denthalpy_D0);
  DECLARE_FIELD(denthalpy_D1);
  DECLARE_FIELD(denthalpy_D2);
  
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  
  interp_s->field = denthalpy_D0;
  plan_interpolation(interp_s);
  dh1[0] = execute_interpolation(interp_s);
  
  interp_s->field = denthalpy_D1;
  plan_interpolation(interp_s);
  dh1[1] = execute_interpolation(interp_s);
  
  interp_s->field = denthalpy_D2;
  plan_interpolation(interp_s);
  dh1[2] = execute_interpolation(interp_s);
  
  /* show method */  
  printf(Pretty0"method = %s\n",par);
  
  /* set functions method */
  parse_adjust_parameter(par,adjust);
  
  void (*force_balance_0)(Physics_T *const phys) = 
            get_func_force_balance_adjustment(adjust[0]);
            
  void (*force_balance_1)(Physics_T *const phys) = 
            get_func_force_balance_adjustment(adjust[1]);
            
  void (*force_balance_2)(Physics_T *const phys) = 
            get_func_force_balance_adjustment(adjust[2]);
  
  if (force_balance_0)
    force_balance_0(phys);
  
  if (force_balance_1)
    force_balance_1(phys);
    
  if (force_balance_2)
    force_balance_2(phys);
    
  _free(adjust[0]);
  _free(adjust[1]);
  _free(adjust[2]);
  
  /* update STRESS_ENERGY and related */
  Seti("enthalpy_neat",0);
  Physics(phys,STRESS_ENERGY);
  
  interp_s->field = denthalpy_D0;
  plan_interpolation(interp_s);
  dh2[0] = execute_interpolation(interp_s);
  
  interp_s->field = denthalpy_D1;
  plan_interpolation(interp_s);
  dh2[1] = execute_interpolation(interp_s);
  
  interp_s->field = denthalpy_D2;
  plan_interpolation(interp_s);
  dh2[2] = execute_interpolation(interp_s);
  
  /* print initial values before adjustments */
  printf(Pretty0"Enthalpy derivatives at NS center before force balance eq.:\n");
  printf(Pretty0"dh(%g,%g,%g)/dx = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh1[0]);
  printf(Pretty0"dh(%g,%g,%g)/dy = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh1[1]);
  printf(Pretty0"dh(%g,%g,%g)/dz = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh1[2]);
    
  /* print initial values after adjustments */
  printf(Pretty0"Enthalpy derivatives at NS center after force balance eq.:\n");
  printf(Pretty0"dh(%g,%g,%g)/dx = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh2[0]);
  printf(Pretty0"dh(%g,%g,%g)/dy = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh2[1]);
  printf(Pretty0"dh(%g,%g,%g)/dz = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh2[2]);
    
  printf(Pretty0"Changes in enthalpy derivatives after force balance eq.:\n");
  printf(Pretty0"dh2/dx-dh1/dx = %+g\n",dh2[0]-dh1[0]);
  printf(Pretty0"dh2/dy-dh1/dy = %+g\n",dh2[1]-dh1[1]);
  printf(Pretty0"dh2/dz-dh1/dz = %+g\n",dh2[2]-dh1[2]);
 
  free_interpolation(interp_s);
  _free(patches);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

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
  
  if(Pcmps("star_NS_extrapolate_matter_fields","poly2"))
  {
    /* make phi, W => enthalpy */
    const char *fields_name[] = {"phi","enthalpy",0};
    
    star_NS_extrapolate(phys,fields_name,"poly2");
    
    phys->region = Ftype("NS_around");
    star_W_spin_vector_idealfluid_update(phys);
    phys_set_region(phys);
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
    
    IF_not_cover(patch,phys) continue;
    
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

/* force_balance_equation : adjust x_CM at direction d/dx */
static void force_balance_ddx_x_CM(Physics_T *const phys)
{
  force_balance_eq_root_finders(phys,0,"x_CM");
}

/* force_balance_equation : adjust x_CM at direction d/dy */
static void force_balance_ddy_x_CM(Physics_T *const phys)
{
  force_balance_eq_root_finders(phys,1,"x_CM");
}

/* force_balance_equation : adjust x_CM at direction d/dz */
static void force_balance_ddz_x_CM(Physics_T *const phys)
{
  force_balance_eq_root_finders(phys,2,"x_CM");
}

/* force_balance_equation : adjust y_CM at direction d/dx */
static void force_balance_ddx_y_CM(Physics_T *const phys)
{
  force_balance_eq_root_finders(phys,0,"y_CM");
}

/* force_balance_equation : adjust y_CM at direction d/dy */
static void force_balance_ddy_y_CM(Physics_T *const phys)
{
  force_balance_eq_root_finders(phys,1,"y_CM");
}

/* force_balance_equation : adjust y_CM at direction d/dz */
static void force_balance_ddz_y_CM(Physics_T *const phys)
{
  force_balance_eq_root_finders(phys,2,"y_CM");
}

/* force_balance_equation : adjust Omega at direction d/dx */
static void force_balance_ddx_Omega(Physics_T *const phys)
{
  force_balance_eq_root_finders(phys,0,"angular_velocity");
}

/* force_balance_equation : adjust Omega at direction d/dy */
static void force_balance_ddy_Omega(Physics_T *const phys)
{
  force_balance_eq_root_finders(phys,1,"angular_velocity");
}

/* force_balance_equation : adjust Omega at direction d/dz */
static void force_balance_ddz_Omega(Physics_T *const phys)
{
  force_balance_eq_root_finders(phys,2,"angular_velocity");
}

/* force_balance_equation : adjust Omega along the V2CM = R_cm-R_NS. */
static void force_balance_ddCM_Omega(Physics_T *const phys)
{
  force_balance_eq_root_finders(phys,-1,"angular_velocity");
}


/* find parameter par using force balance equation in direction dir */
static void force_balance_eq_root_finders(Physics_T *const phys,const int dir, const char *const par)
{
  const double D            = sysGetd("separation");
  const double Vr           = sysGetd("infall_velocity");
  const double NS_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  const double RESIDUAL     = sqrt(Getd("RootFinder_Tolerance"));
  const double Omega_BHNS   = sysGetd("angular_velocity");
  const double x_CM         = sysGetd("x_CM");
  const double y_CM         = sysGetd("y_CM");
  const double Scale        = 0.1;/* scale the weight */
  const double Rel_Change   = 0.1;/* relative change */
  double W1  = Getd("Force_Balance_Update_Weight");
  double W2  = 1-W1;
  double *new_par,old_par = 0;
  double guess[1],X[3],V2CM[2]={0,0};
  struct Force_Balance_RootFinder_S params[1] = {0};
  Patch_T *patch    = 0;
  Patch_T **patches = 0;
  char desc[1000] = {'\0'};
  unsigned Np;
  char reg[99];
  
  sprintf(reg,"%s_cent.*",phys->spos);
  patches = regex_collect_patches(phys->grid,reg,&Np);
  
  /* initial values before adjustments */
  patch = x_in_which_patch(NS_center,patches,Np);
  assert(patch);
  assert(X_of_x(X,NS_center,patch));
  
  params->patch      = patch;
  params->X          = X;
  params->Vr         = Vr;
  params->D          = D;
  params->V2CM       = V2CM;
  params->find_x_CM  = 0;
  params->find_y_CM  = 0;
  params->find_Omega = 0;
  params->dir        = dir;
  
  /* which par */
  if (strcmp_i("angular_velocity",par))
  {
    params->find_Omega = 1;
    params->x_CM       = x_CM;
    params->y_CM       = y_CM;
    guess[0]           = Omega_BHNS;
    old_par            = Omega_BHNS;
  }
  else if (strcmp_i("x_CM",par))
  {
    params->find_x_CM  = 1;
    params->Omega_BHNS = Omega_BHNS;
    params->y_CM       = y_CM;
    guess[0]           = x_CM;
    old_par            = x_CM;
  }
  else if (strcmp_i("y_CM",par))
  {
    params->find_y_CM  = 1;
    params->Omega_BHNS = Omega_BHNS;
    params->x_CM       = x_CM;
    guess[0]           = y_CM;
    old_par            = y_CM;
  }
  else
    Error0(NO_OPTION);
  
  /* which direction */
  if (dir == -1)
  {
    V2CM[0] = x_CM-NS_center[0];
    V2CM[1] = y_CM-NS_center[1];
    double V2CM_norm = root_square(2,V2CM,0);
    assert(!EQL(V2CM_norm,0));
    V2CM[0] /= V2CM_norm;
    V2CM[1] /= V2CM_norm;
    params->dLnGamma = V2CM[0]*star_NS_idealfluid_gConf_dLnGamma_force_bal(patch,X,0) + 
                       V2CM[1]*star_NS_idealfluid_gConf_dLnGamma_force_bal(patch,X,1);
  }
  else if (dir == 0)
  {
    params->dLnGamma = star_NS_idealfluid_gConf_dLnGamma_force_bal(patch,X,0);
  }
  else if (dir == 1)
  {
    params->dLnGamma = star_NS_idealfluid_gConf_dLnGamma_force_bal(patch,X,1);
  }
  else if (dir == 2)
  {
    params->dLnGamma = star_NS_idealfluid_gConf_dLnGamma_force_bal(patch,X,2);
  }
  else
    Error0(NO_OPTION);
  
  sprintf(desc,"Solving Force Balance Eq. for '%s' at direction 'x^%d'",par,dir);
  
  Root_Finder_T *root = init_root_finder(1);
  root->description   = desc;
  root->verbose       = 1;
  root->type          = Gets("RootFinder_Method");
  root->tolerance     = Getd("RootFinder_Tolerance");
  root->MaxIter       = (unsigned)Geti("RootFinder_Max_Number_of_Iteration");
  root->x_gss         = guess;
  root->params        = params;
  
  Error0("WHAT TO DO for B?");//root->f[0]          = force_balance_root_finder_eq;
  
  plan_root_finder(root);
  
  new_par = execute_root_finder(root);
  
  if (root->exit_status != ROOT_FINDER_OK && GRT(root->residual,RESIDUAL))
  {
    print_root_finder_exit_status(root);
  }
  
  /* if change is big try to soften it */
  if (fabs(1-new_par[0]/old_par) > Rel_Change)
  {
    W1 *= Scale;
    W2  = 1-W1;
  }
  new_par[0] = W1*new_par[0]+W2*old_par;
  Psetd(par,new_par[0]);
  
  /* since B1 has been changed let's update the pertinent fields */
  
  Error0("WHAT TO DO for B?");//update_B1_dB1_Beta_dBete_Aij_dAij(grid);
  
  free_root_finder(root);
  free(new_par);
  _free(patches);
}

/* getting adjustment str, returns the relevant function. */
static fAdjustment_t *get_func_force_balance_adjustment(const char *const adjust)
{
  fAdjustment_t *f = 0;
  
  if (!adjust)
  {
    f = 0;
  }
  else if (strcmp_i(adjust,"none"))
  {
    f = 0;
  }
  else if (strcmp_i(adjust,"d/dx:x_CM"))
  {
    f = force_balance_ddx_x_CM;
  }
  else if (strcmp_i(adjust,"d/dCM:Omega"))
  {
    f = force_balance_ddCM_Omega;
  }
  else if (strcmp_i(adjust,"d/dy:x_CM"))
  {
    f = force_balance_ddy_x_CM;
  }
  else if (strcmp_i(adjust,"d/dz:x_CM"))
  {
    f = force_balance_ddz_x_CM;
  }
  else if (strcmp_i(adjust,"d/dx:y_CM"))
  {
    f = force_balance_ddx_y_CM;
  }
  else if (strcmp_i(adjust,"d/dy:y_CM"))
  {
    f = force_balance_ddy_y_CM;
  }
  else if (strcmp_i(adjust,"d/dz:y_CM"))
  {
    f = force_balance_ddz_y_CM;
  }
  else if (strcmp_i(adjust,"d/dx:Omega"))
  {
    f = force_balance_ddx_Omega;
  }
  else if (strcmp_i(adjust,"d/dy:Omega"))
  {
    f = force_balance_ddy_Omega;
  }
  else if (strcmp_i(adjust,"d/dz:Omega"))
  {
    f = force_balance_ddz_Omega;
  }
  else
    Error0(NO_OPTION);
  
  return f;
}

/* parsing adjust parameter value and fill components consequently */
static void parse_adjust_parameter(const char *const par,char *adjust[3])
{
  if (!strstr_i(par,"adjust(") && !strstr_i(par,"none"))
    Error1("Syntax error for '%s'.\n",par);
  
  /* if it is none */  
  if (strcmp_i(par,"none"))
  {
    adjust[0] = 0;
    adjust[1] = 0;
    adjust[2] = 0;
    
    return;
  }
  
  /* parse if not none */
  char *str = dup_s(par);
  char *save;
  char *main_str = sub_s(str,'(',')',&save);
  char *tok  = tok_s(main_str,',',&save);
  
  adjust[0] = dup_s(tok);
  tok = tok_s(0,',',&save);
  adjust[1] = dup_s(tok);
  tok = tok_s(0,',',&save);
  adjust[2] = dup_s(tok);
  
  free(str);
}
