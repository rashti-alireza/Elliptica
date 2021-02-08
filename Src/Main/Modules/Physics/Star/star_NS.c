/*
// Alireza Rashti
// November 2020
*/


/* collection of NS star affairs */

#include "star_NS.h"

/* adjust NS center to be fixed at a specific location otherwise
// the star might drift away. */
int star_NS_keep_center_fixed(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,"NS");
  Patch_T *patch     = 0;
  const double NS_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  Interpolation_T *interp_s = init_interpolation();
  double dh1[3] = {0},dh2[3] = {0}, X[3] = {0};
  
  /* initial values before adjustments */
  patch = x_in_which_patch(NS_center,grid->patch,grid->np);
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

  IF_sval("adjust_center_method","Interpolation")
  {
    adjust_NS_center_interpolation(phys);
  }
  else IF_sval("adjust_center_method","Taylor_expansion")
  {
    adjust_NS_center_Taylor_expansion(phys);
  }
  else
    Error0(NO_OPTION);

  /* note: enthalpy is already updated! */
  
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
  printf(Pretty0"Enthalpy derivatives at NS center before adjustment:\n");
  printf(Pretty0"dh(%g,%g,%g)/dx = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh1[0]);
  printf(Pretty0"dh(%g,%g,%g)/dy = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh1[1]);
  printf(Pretty0"dh(%g,%g,%g)/dz = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh1[2]);
    
  /* print initial values after adjustments */
  printf(Pretty0"Enthalpy derivatives at NS center after adjustment:\n");
  printf(Pretty0"dh(%g,%g,%g)/dx = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh2[0]);
  printf(Pretty0"dh(%g,%g,%g)/dy = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh2[1]);
  printf(Pretty0"dh(%g,%g,%g)/dz = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh2[2]);
    
  printf(Pretty0"Changes in enthalpy derivatives after adjustment:\n");
  printf(Pretty0"dh2/dx-dh1/dx = %+g\n",dh2[0]-dh1[0]);
  printf(Pretty0"dh2/dy-dh1/dy = %+g\n",dh2[1]-dh1[1]);
  printf(Pretty0"dh2/dz-dh1/dz = %+g\n",dh2[2]-dh1[2]);
 
  free_interpolation(interp_s);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}  

/* adjust force balance */
int star_NS_idealfluid_gConf_force_balance(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,"NS");
  const double NS_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  const char *const par = Gets("force_balance_equation");
  char *adjust[3] = {0};
  double dh1[3] = {0};
  double dh2[3] = {0}; 
  double X[3]   = {0};
  
  /* show info */  
  printf(Pretty0"method        = %s\n",par);
  printf(Pretty0"update weight = %e\n",
         Getd("force_balance_update_weight"));
  
  /* set functions method */
  parse_adjust_parameter(par,adjust);
  
  void (*force_balance_0)(Physics_T *const phys) = 
            get_func_force_balance_adjustment(adjust[0]);
            
  void (*force_balance_1)(Physics_T *const phys) = 
            get_func_force_balance_adjustment(adjust[1]);
            
  void (*force_balance_2)(Physics_T *const phys) = 
            get_func_force_balance_adjustment(adjust[2]);
  
  /* initial values before adjustments */
  if (0)//force_balance_0 || force_balance_1 || force_balance_2)
  {
    Interpolation_T *interp_s = init_interpolation();
    Patch_T *patch = x_in_which_patch(NS_center,grid->patch,grid->np);
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
    
    free_interpolation(interp_s);
  }
  
  /* execute force balance */
  if (force_balance_0)
    force_balance_0(phys);
  
  if (force_balance_1)
    force_balance_1(phys);
    
  if (force_balance_2)
    force_balance_2(phys);
    
  Free(adjust[0]);
  Free(adjust[1]);
  Free(adjust[2]);
  
  /* update stress energy tensor and related and show difference */
  if (0)//force_balance_0 || force_balance_1 || force_balance_2)
  {
    Sets("enthalpy_neat","no");
    physics(phys,STRESS_ENERGY_UPDATE);
    
    Interpolation_T *interp_s = init_interpolation();
    Patch_T *patch = x_in_which_patch(NS_center,grid->patch,grid->np);
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
    dh2[0] = execute_interpolation(interp_s);
    
    interp_s->field = denthalpy_D1;
    plan_interpolation(interp_s);
    dh2[1] = execute_interpolation(interp_s);
    
    interp_s->field = denthalpy_D2;
    plan_interpolation(interp_s);
    dh2[2] = execute_interpolation(interp_s);
    
    /* print initial values before adjustments */
    printf(Pretty0"Enthalpy derivatives at NS center before "
                  "force balance eq.:\n");
    printf(Pretty0"dh(%g,%g,%g)/dx = %+g\n",
      NS_center[0],NS_center[1],NS_center[2],dh1[0]);
    printf(Pretty0"dh(%g,%g,%g)/dy = %+g\n",
      NS_center[0],NS_center[1],NS_center[2],dh1[1]);
    printf(Pretty0"dh(%g,%g,%g)/dz = %+g\n",
      NS_center[0],NS_center[1],NS_center[2],dh1[2]);
      
    /* print initial values after adjustments */
    printf(Pretty0"Enthalpy derivatives at NS center after force "
                  "balance eq.:\n");
    printf(Pretty0"dh(%g,%g,%g)/dx = %+g\n",
      NS_center[0],NS_center[1],NS_center[2],dh2[0]);
    printf(Pretty0"dh(%g,%g,%g)/dy = %+g\n",
      NS_center[0],NS_center[1],NS_center[2],dh2[1]);
    printf(Pretty0"dh(%g,%g,%g)/dz = %+g\n",
      NS_center[0],NS_center[1],NS_center[2],dh2[2]);
      
    printf(Pretty0"Changes in enthalpy derivatives after "
                  "force balance eq.:\n");
    printf(Pretty0"dh2/dx-dh1/dx = %+g\n",dh2[0]-dh1[0]);
    printf(Pretty0"dh2/dy-dh1/dy = %+g\n",dh2[1]-dh1[1]);
    printf(Pretty0"dh2/dz-dh1/dz = %+g\n",dh2[2]-dh1[2]);
  }
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* find Euler equation constant to meet NS baryonic mass */
int star_NS_ifluid_gConf_find_EulerC_fix_baryon_mass(Physics_T *const phys)
{
  FUNC_TIC
  
  Root_Finder_T *root = init_root_finder(1);
  const double W1  = Getd("Euler_const_update_weight");
  const double W2  = 1-W1;
  double *Euler_const = 0;
  double guess[1] = {Getd("Euler_equation_constant")};
  const double Residual = sqrt(Getd("RootFinder_Tolerance"));
  struct NS_Euler_eq_const_RootFinder_S params[1] = {0};
  double bar_mass,adm_mass,Komar_mass;
  
  bar_mass = star_NS_baryonic_gConf_mass(phys,guess[0]);
  observe(phys,"ADM(M)"  ,Gets("Observe_ADM_M")  ,&adm_mass);
  observe(phys,"Komar(M)",Gets("Observe_Komar_M"),&Komar_mass);

  printf(Pretty0"current NS baryonic mass = %e\n",bar_mass);
  printf(Pretty0"current NS ADM mass      = %e\n",adm_mass);
  printf(Pretty0"current NS Komar mass    = %e\n",Komar_mass);
  
  Setd("baryonic_mass_current",bar_mass);
  Setd("ADM_mass",adm_mass);
  Setd("Komar_mass",Komar_mass);
  
  params->phys = phys;
  params->NS_baryonic_mass = Getd("baryonic_mass");
  
  root->type        = Gets("RootFinder_Method");
  root->tolerance   = Getd("RootFinder_Tolerance");
  root->MaxIter     = (Uint)Geti("RootFinder_Iteration");
  root->x_gss       = guess;
  root->params      = params;
  root->f[0]        = Euler_eq_const_gConf_rootfinder_eq;
  root->verbose     = strstr_i(Gets("RootFinder_verbose"),"yes");
  plan_root_finder(root);
  Euler_const       = execute_root_finder(root);
  
  printf(Pretty0"update weight = %e\n",W1);
  /* if root finder is not OK for some reason */
  if (GRT(root->residual,Residual))
    Euler_const[0] = guess[0];/* don't update */
  else
    Euler_const[0] = W1*Euler_const[0]+W2*guess[0];
  
  Setd("Euler_equation_constant",Euler_const[0]);
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
  
  IF_sval("extrapolate_matter_fields","inverse_r_expmr")
  {
    /* make phi, W => enthalpy */
    const char *fields_name[] = {"phi","enthalpy",0};
    
    star_NS_extrapolate(phys,fields_name,"inverse_r_expmr");
    
    star_W_spin_vector_idealfluid_update(phys,"NS_around");
  }
  else IF_sval("extrapolate_matter_fields","inverse_r_expmAr")
  {
    /* make phi, W => enthalpy */
    const char *fields_name[] = {"phi","enthalpy",0};
    
    star_NS_extrapolate(phys,fields_name,"inverse_r_expmAr");
    
    star_W_spin_vector_idealfluid_update(phys,"NS_around");
  }
  else IF_sval("extrapolate_matter_fields","inverse_r2_expmAr")
  {
    /* make phi, W => enthalpy */
    const char *fields_name[] = {"phi","enthalpy",0};
    
    star_NS_extrapolate(phys,fields_name,"inverse_r2_expmAr");
    
    star_W_spin_vector_idealfluid_update(phys,"NS_around");
  }
  else IF_sval("extrapolate_matter_fields","enthalpy_expmr_phi_inverse_r2")
  {
    const char *fields_name[2] = {0,0};
    
    fields_name[0] = "enthalpy";
    star_NS_extrapolate(phys,fields_name,"expmr");
    
    fields_name[0] = "phi";
    star_NS_extrapolate(phys,fields_name,"inverse_r2");
    
    star_W_spin_vector_idealfluid_update(phys,"NS_around");
  }
  else IF_sval("extrapolate_matter_fields","exp2")
  {
    /* make phi, W => enthalpy */
    const char *fields_name[] = {"phi","enthalpy",0};
    
    star_NS_extrapolate(phys,fields_name,"exp2");
    
    star_W_spin_vector_idealfluid_update(phys,"NS_around");
  }
  else IF_sval("extrapolate_matter_fields","poly2")
  {
    /* make phi, W => enthalpy */
    const char *fields_name[] = {"phi","enthalpy",0};
    
    star_NS_extrapolate(phys,fields_name,"poly2");
    
    star_W_spin_vector_idealfluid_update(phys,"NS_around");
  }
  else IF_sval("extrapolate_matter_fields","inverse_r2")
  {
    /* make phi, W => enthalpy */
    const char *fields_name[] = {"phi","enthalpy",0};
    
    star_NS_extrapolate(phys,fields_name,"inverse_r2");
    
    star_W_spin_vector_idealfluid_update(phys,"NS_around");
  }
  else IF_sval("extrapolate_matter_fields","inverse_r2_expmr")
  {
    /* make phi, W => enthalpy */
    const char *fields_name[] = {"phi","enthalpy",0};
    
    star_NS_extrapolate(phys,fields_name,"inverse_r2_expmr");
    
    star_W_spin_vector_idealfluid_update(phys,"NS_around");
  }
  else
  {
    Error0(NO_OPTION);
  }
  
  FUNC_TOC
  return EXIT_SUCCESS;
}


/* calculate W = Omega_NS x (r-C_NS) */
void star_W_spin_vector_idealfluid_update(Physics_T *const phys,
                                          const char *const region)
{
  FUNC_TIC
  AssureType(phys->ctype == NS);
  
  Grid_T *const grid = mygrid(phys,region);
  double Omega_NS[3],C_NS[3];
  
  Omega_NS[0] = Getd("Omega_x");
  Omega_NS[1] = Getd("Omega_y");
  Omega_NS[2] = Getd("Omega_z");
  C_NS[0]     = Getd("center_x");
  C_NS[1]     = Getd("center_y");
  C_NS[2]     = Getd("center_z");
  
  FOR_ALL_p(grid->np)
  {
    Patch_T *patch = grid->patch[p];
    
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
  Uint nn = patch->nn;
  Uint ijk;
  
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
  Grid_T *const grid        = mygrid(phys,"NS");
  const double D            = sysGetd("separation");
  const double Vr           = sysGetd("infall_velocity");
  const double NS_center[3] = {Getd("center_x"),
                               Getd("center_y"),
                               Getd("center_z")};
  const double Residual     = sqrt(Getd("RootFinder_Tolerance"));
  const double Omega_sys    = sysGetd("angular_velocity");
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
  char s[999] = {'\0'};
  
  printf(Pretty0"Solving Force Balance Eq. for '%s' at direction 'x^%d'\n",
         par,dir);
  
  /* initial values before adjustments */
  patch = x_in_which_patch(NS_center,grid->patch,grid->np);
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
    guess[0]           = Omega_sys;
    old_par            = Omega_sys;
  }
  else if (strcmp_i("x_CM",par))
  {
    params->find_x_CM  = 1;
    params->Omega      = Omega_sys;
    params->y_CM       = y_CM;
    guess[0]           = x_CM;
    old_par            = x_CM;
  }
  else if (strcmp_i("y_CM",par))
  {
    params->find_y_CM  = 1;
    params->Omega      = Omega_sys;
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
  
  Root_Finder_T *root = init_root_finder(1);
  root->type          = Gets("RootFinder_Method");
  root->tolerance     = Getd("RootFinder_Tolerance");
  root->MaxIter       = (Uint)Geti("RootFinder_Iteration");
  root->x_gss         = guess;
  root->params        = params;
  root->f[0]          = star_NS_idealfluid_gConf_root_force_bal;
  root->verbose       = strstr_i(Gets("RootFinder_verbose"),"yes");
  plan_root_finder(root);
  new_par = execute_root_finder(root);
  
  if (root->exit_status != ROOT_FINDER_OK && GRT(root->residual,Residual))
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
  
  /* update parameter */
  sprintf(s,"%s_%s",phys->ssys,par);
  Psetd(s,new_par[0]);
  
  /* update */
  //physics(phys,ADM_UPDATE_B1I);
  //physics(phys,ADM_UPDATE_beta);
  
  free_root_finder(root);
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
    Errors("Syntax error for '%s'.\n",par);
  
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

/* adjust the center of NS at the designated point, in case it moved.
// we need only to draw enthalpy to NS center.
// to do so, we demand shifted_enthalpy(r) = enthalpy(dR+r), in which dR 
// is the amount the center is displaced from NS center (dR = R-r_cen);
// thus f_new(center) = f_old(center+dR) => f_new(r) = f_old(r+dR).
// finally update the enthalpy and its derivatives. */
static void adjust_NS_center_interpolation(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid  = mygrid(phys,"NS,NS_around");
  const double W1     = Getd("adjust_center_update_weight");
  const double W2     = 1.-W1;
  double NS_center[3] = {Getd("center_x"),
                         Getd("center_y"),
                         Getd("center_z")};
  double R[3]  = {0};
  double dR[3] = {0};
  Uint p;
  
  star_NS_find_where_denthalpy_is_0(phys,R);
  
  dR[0] = R[0]-NS_center[0];
  dR[1] = R[1]-NS_center[1];
  dR[2] = R[2]-NS_center[2];
  
  /* if it is already located at the designted point */
  if (EQL(dR[0],0.) && EQL(dR[1],0.) && EQL(dR[2],0.))
  {
    FUNC_TOC
    return;
  }
    
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    const double *h_old = 0;
    double x[3],Xp[3];
    
    /* now shift enthalpy */
    Patch_T *patchp = 0;
    DECLARE_FIELD(enthalpy);
    ADD_FIELD(shifted_enthalpy);
    REALLOC_v_WRITE_v(shifted_enthalpy);
    make_coeffs_3d(enthalpy);
    h_old = enthalpy->v;
    
    FOR_ALL_ijk
    {
      x[0] = patch->node[ijk]->x[0]+dR[0];
      x[1] = patch->node[ijk]->x[1]+dR[1];
      x[2] = patch->node[ijk]->x[2]+dR[2];
      patchp = x_in_which_patch(x,grid->patch,grid->np);
      
      if (patchp) 
      {
        X_of_x(Xp,x,patchp);
        shifted_enthalpy[ijk] = W2*h_old[ijk]+W1*f_of_X("enthalpy",Xp,patchp);
      }
      else
      {
        /* try */
        patchp = x_in_which_patch_force(x,grid->patch,grid->np,Xp);
        if (patchp)
          shifted_enthalpy[ijk] = W2*h_old[ijk]+W1*f_of_X("enthalpy",Xp,patchp);
        /* if point x located outside of NS around */
        else
          shifted_enthalpy[ijk] = h_old[ijk];
      }
    }
  }
  
  /* now clean enthalpy and copy new value and remove extras */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
 
    DECLARE_FIELD(enthalpy);
    DECLARE_FIELD(shifted_enthalpy);
    free_coeffs(enthalpy);
    free(enthalpy->v);
    enthalpy->v = shifted_enthalpy->v;
    shifted_enthalpy->v = 0;
    REMOVE_FIELD(shifted_enthalpy);
  }
  
  /* update enthalpy derivatives */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    dField_di(denthalpy_D0);
    dField_di(denthalpy_D1);
    dField_di(denthalpy_D2);
  }
  
  FUNC_TOC
}

/* find the NS center (coords where denthalpy is 0) 
// using d(enthalpy)/dx^i = 0 and put the coords location into xh0. */
void star_NS_find_where_denthalpy_is_0(Physics_T *const phys,double xdh0[3])
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,"NS");
  double *NS_center;
  struct NC_Center_RootFinder_S params[1] = {0};
  const double Residual = sqrt(Getd("RootFinder_Tolerance"));
  double guess[3] = {Getd("center_x"),
                     Getd("center_y"),
                     Getd("center_z")};
  
  Root_Finder_T *root = init_root_finder(3);                   
  params->root_finder = root;
  root->type        = Gets("RootFinder_Method");
  root->tolerance   = Getd("RootFinder_Tolerance");
  root->MaxIter     = (Uint)Geti("RootFinder_Iteration");
  root->x_gss       = guess;
  root->params      = params;
  root->f[0]        = dh_dx0_root_finder_eq;
  root->f[1]        = dh_dx1_root_finder_eq;
  root->f[2]        = dh_dx2_root_finder_eq;    
  params->patches   = grid->patch;
  params->Np        = grid->np;
  root->verbose     = strstr_i(Gets("RootFinder_verbose"),"yes");

  plan_root_finder(root);
  NS_center = execute_root_finder(root);
    
  /* if root finder was successful */
  if (root->exit_status == ROOT_FINDER_OK || LSS(root->residual,Residual))
  {
    xdh0[0] = NS_center[0];
    xdh0[1] = NS_center[1];
    xdh0[2] = NS_center[2];
    
    printf(Pretty0"Current NS center found at (%g,%g,%g)\n",
                  NS_center[0],NS_center[1],NS_center[2]);
    printf(Pretty0"Change in x direction = %+g\n",NS_center[0]-guess[0]);
    printf(Pretty0"Change in y direction = %+g\n",NS_center[1]-guess[1]);
    printf(Pretty0"Change in z direction = %+g\n",NS_center[2]-guess[2]);
  }
  else
  {
    print_root_finder_exit_status(root);
    xdh0[0] = guess[0];
    xdh0[1] = guess[1];
    xdh0[2] = guess[2];
  }
  
  free_root_finder(root);
  
  FUNC_TOC
}

/* dh/dx^0 = 0 */
static double dh_dx0_root_finder_eq(void *params,const double *const x)
{
  struct NC_Center_RootFinder_S *const par = params;
  Patch_T **const patches = par->patches;
  const Uint Np           = par->Np;
  Patch_T *const patch    = x_in_which_patch(x,patches,Np);
  Interpolation_T *interp_s;
  double interp,X[3];

  /* if this point is out of this patches, exit */
  if (!patch)
  {
    par->root_finder->interrupt = 1;
    return 0;
  }
  
  DECLARE_FIELD(denthalpy_D0);
  
  X_of_x(X,x,patch);
  interp_s = init_interpolation();
  interp_s->field = denthalpy_D0;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}

/* dh/dx^1 = 0 */
static double dh_dx1_root_finder_eq(void *params,const double *const x)
{
  struct NC_Center_RootFinder_S *const par = params;
  Patch_T **const patches = par->patches;
  const Uint Np           = par->Np;
  Patch_T *const patch    = x_in_which_patch(x,patches,Np);
  Interpolation_T *interp_s;
  double interp,X[3];

  /* if this point is out of this patches, exit */
  if (!patch)
  {
    par->root_finder->interrupt = 1;
    return 0;
  }
  
  DECLARE_FIELD(denthalpy_D1);
  
  X_of_x(X,x,patch);
  interp_s = init_interpolation();
  interp_s->field = denthalpy_D1;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}

/* dh/dx^2 = 0 */
static double dh_dx2_root_finder_eq(void *params,const double *const x)
{
  struct NC_Center_RootFinder_S *const par = params;
  Patch_T **const patches = par->patches;
  const Uint Np           = par->Np;
  Patch_T *const patch    = x_in_which_patch(x,patches,Np);
  Interpolation_T *interp_s;
  double interp,X[3];

  /* if this point is out of this patches, exit */
  if (!patch)
  {
    par->root_finder->interrupt = 1;
    return 0;
  }
  
  DECLARE_FIELD(denthalpy_D2);
  
  X_of_x(X,x,patch);
  interp_s = init_interpolation();
  interp_s->field = denthalpy_D2;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}


/* adjust the center of NS at the designated point, in case it moved. 
// in this method, we tune enthalpy values such that the derivative of 
// the enthalpy be 0 at NS center, using Taylor expansion. */
static void adjust_NS_center_Taylor_expansion(Physics_T *const phys)
{
  FUNC_TIC
  
  Grid_T *const grid  = mygrid(phys,"NS,NS_around");
  const double W      = Getd("adjust_center_update_weight");
  double NS_center[3] = {Getd("center_x"),
                         Getd("center_y"),
                         Getd("center_z")};
  Patch_T *const patch_center = 
        x_in_which_patch(NS_center,grid->patch,grid->np);
  double R[3]  = {0};
  double dR[3] = {0};
  double dh[3] = {0};
  double Xc[3] = {0};
  Uint p;
  
  star_NS_find_where_denthalpy_is_0(phys,R);
  
  dR[0] = R[0]-NS_center[0];
  dR[1] = R[1]-NS_center[1];
  dR[2] = R[2]-NS_center[2];
  
  /* if it is already located at the designted point */
  if (EQL(dR[0],0.) && EQL(dR[1],0.) && EQL(dR[2],0.))
  {
    FUNC_TOC
    return;
  }
  
  if (!patch_center)
    Error0("could not find the pertinent patch!");
  
  /* denthalpy/d? at the center */
  X_of_x(Xc,NS_center,patch_center);
  dh[0] = f_of_X("denthalpy_D0",Xc,patch_center);
  dh[1] = f_of_X("denthalpy_D1",Xc,patch_center);
  dh[2] = f_of_X("denthalpy_D2",Xc,patch_center);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    WRITE_v(enthalpy)
  
    FOR_ALL_ijk
    {
      double x = patch->node[ijk]->x[0];
      double y = patch->node[ijk]->x[1];
      double z = patch->node[ijk]->x[2];
      
      enthalpy[ijk] -= W*(dh[0]*(x-NS_center[0])+
                          dh[1]*(y-NS_center[1])+
                          dh[2]*(z-NS_center[2]));
    }
    
    dField_di(denthalpy_D0);
    dField_di(denthalpy_D1);
    dField_di(denthalpy_D2);
  }
  
  FUNC_TOC
}

/* populate psi, alphaPsi, enthalpy, phi and W fields using
// TOV solution (stems must be given) */
void star_populate_psi_alphaPsi_matter_fields_TOV
      (Physics_T *const phys,const char *const region,
      const char *const Psi,const char *const AlphaPsi,
      const char *const Enthalpy,const char *const Rho0,
      const char *const Phi,const char *const W)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  const double NSx   = Getd("center_x");
  const double NSy   = Getd("center_y");
  const double NSz   = Getd("center_z");
  const double O_x   = Getd("omega_x");
  const double O_y   = Getd("omega_y");
  const double O_z   = Getd("omega_z");
  const double Omega = sysGetd("angular_velocity");
  const double y_CM  = sysGetd("y_CM");
  
  /* find TOV solution */
  TOV_T *tov = TOV_init();
  tov->phys  = phys;
  tov->bar_m = Getd("baryonic_mass");
  tov = TOV_solution(tov);
  const double R_Schwar = tov->r[tov->N-1];/* NS's Schwarzchild radius */
  const double M_NS  = tov->ADM_m;/* NS adm mass */
  
  FOR_ALL_p(grid->np)
  {
    Patch_T *patch = grid->patch[p];
    
    if (IsItCovering(patch,"NS"))
    {
      REALLOC_v_WRITE_v_STEM(psi,Psi)
      REALLOC_v_WRITE_v_STEM(alphaPsi,AlphaPsi)
      REALLOC_v_WRITE_v_STEM(enthalpy,Enthalpy)
      REALLOC_v_WRITE_v_STEM(rho0,Rho0)
      REALLOC_v_WRITE_v_STEM(phi,Phi)
      REALLOC_v_WRITE_v_STEM(w_U0,W)
      REALLOC_v_WRITE_v_STEM(w_U1,W)
      REALLOC_v_WRITE_v_STEM(w_U2,W)
      
      Interpolation_T *interp_psi = init_interpolation();
      interp_psi->method          = "Natural_Cubic_Spline_1D";
      interp_psi->N_cubic_spline_1d->f = tov->psi;
      interp_psi->N_cubic_spline_1d->x = tov->rbar;
      interp_psi->N_cubic_spline_1d->N = tov->N;
      plan_interpolation(interp_psi);

      Interpolation_T *interp_h = init_interpolation();
      interp_h->method          = "Natural_Cubic_Spline_1D";
      interp_h->N_cubic_spline_1d->f = tov->h;
      interp_h->N_cubic_spline_1d->x = tov->rbar;
      interp_h->N_cubic_spline_1d->N = tov->N;
      plan_interpolation(interp_h);
      
      EoS_T *eos = init_EoS(phys);

      FOR_ALL_ijk
      {
        double x = patch->node[ijk]->x[0]-NSx;
        double y = patch->node[ijk]->x[1]-NSy;
        double z = patch->node[ijk]->x[2]-NSz;
        double r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
        double alpha;
        double enthalpy_h;

        interp_psi->N_cubic_spline_1d->h = r;
        interp_h->N_cubic_spline_1d->h   = r;

        /* psi */
        psi[ijk] = execute_interpolation(interp_psi);

        /* alphaPsi */
        enthalpy_h = execute_interpolation(interp_h);
        alpha = sqrt(1-2*M_NS/R_Schwar)/enthalpy_h;
        alphaPsi[ijk] = psi[ijk]*alpha;

        /* enthalpy */
        enthalpy[ijk] = enthalpy_h;

        /* rho0 */
        eos->h = enthalpy_h;
        rho0[ijk] = eos->rest_mass_density(eos);

        /* phi Newtonian approximation */
        phi[ijk] = -Omega*(NSy-y_CM)*x;

        /* spin part */
        w_U0[ijk] = O_y*z-O_z*y;
        w_U1[ijk] = O_z*x-O_x*z;
        w_U2[ijk] = O_x*y-O_y*x;
      }
      free_interpolation(interp_psi);
      free_interpolation(interp_h);
      free_EoS(eos);
    }
    else
    {
      REALLOC_v_WRITE_v_STEM(psi,Psi)
      REALLOC_v_WRITE_v_STEM(alphaPsi,AlphaPsi)

      FOR_ALL_ijk
      {
        double x    = patch->node[ijk]->x[0]-NSx;
        double y    = patch->node[ijk]->x[1]-NSy;
        double z    = patch->node[ijk]->x[2]-NSz;
        double r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
        double alpha;

        /* psi */
        psi[ijk] = 1+0.5*M_NS/r;

        /* alphaPsi */
        alpha = (1-0.5*M_NS/r)/(1+0.5*M_NS/r);
        alphaPsi[ijk] = psi[ijk]*alpha;
      }
    }
  }
  
  TOV_free(tov);
  
  FUNC_TOC
}
