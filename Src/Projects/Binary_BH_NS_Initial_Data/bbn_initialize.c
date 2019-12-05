/*
// Alireza Rashti
// June 2019
*/

#include "bbn_initialize.h"

/* initialize this system according to the previous grid. 
// ->return value: the next grid as a result of this initialization. */
Grid_T *bbn_initialize_next_grid(Grid_T *const grid_prev)
{
  Grid_T *grid_next = 0;
  
  if (!grid_prev)/* if grid is empty come up with an approximation */
  {
    /* if we use TOV and Kerr-Schil black hole approximation */
    if (strcmp_i(GetParameterS_E("BH_NS_initialization"),"TOV_KerrSchild"))
      grid_next = TOV_KerrSchild_approximation();
    else
      abortEr(NO_OPTION);
  }
  else/* use previous grid to make the next one with new adjustments */
  {
    grid_next = make_next_grid_using_previous_grid(grid_prev);
  }
  
  return grid_next;   
}

/* finding different quantities and then make the next grid using previous grid
// first find the values of the following parameters and some adjustment:
// . orbital angular velocity
// . center of rotation (center of mass)
// . Euler equation constant.
// . NS center
// . drag NS to the center
// . find NS surface
// . BH_radius
// . Omega_BH
// ->return value: the next grid called 'grid_next' */
static Grid_T *make_next_grid_using_previous_grid(Grid_T *const grid_prev)
{
  Grid_T *grid_next = 0;
  struct Grid_Params_S *GridParams = init_GridParams();/* adjust some pars for construction of next grid */
  
  /* NOTE: the order of function calls are important */

  /* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S */
  bbn_update_stress_energy_tensor(grid_prev,0);
  
  /* find Euler equation constant to meet NS baryonic mass */
  find_Euler_eq_const(grid_prev);
  
  /* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S */
  bbn_update_stress_energy_tensor(grid_prev,0);
  
  /* find y_CM or orbital_angular_velocity using force balance equation */
  force_balance_eq(grid_prev);
  
  /* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S */
  bbn_update_stress_energy_tensor(grid_prev,0);
  
  /* extrapolate fluid fields outside NS */
  extrapolate_fluid_fields_outsideNS(grid_prev);
  
  /* keep NS center fixed at (0,C,0) by adjusting enthalpy */
  keep_NS_center_fixed(grid_prev);
  
  /* find NS surface using h = 1 */
  find_NS_surface(grid_prev,GridParams);
  
  /* P_ADM control */
  P_ADM_control(grid_prev);
  
  /* find the apparent horizon radius to acquire the desired BH mass */
  find_AH_radius(grid_prev,GridParams);
  
  /* find the Omega_BH to acquire the desired BH spin */
  //find_BH_Omega(grid_prev,GridParams);
  
  /* make new grid with new parameters */
  const double bh_chi  = GetParameterD_E("BH_X_U2");
  const double bh_mass = GetParameterD_E("BH_mass");
  //const double bh_R    = bh_mass*(1+sqrt(1-SQR(bh_chi)));
  //GridParams->R_BH_r = bh_R;
  GridParams->a_BH   = bh_chi*bh_mass;
  GridParams->BH_R_type = "PerfectSphere";
  grid_next = creat_bbn_grid_CS(GridParams);
  
  /* fields: */
  /* creating all of the fields needed for construction of Initial Data */
  bbn_allocate_fields(grid_next);
  
  /* populating the free data part of initial data that we chose ourself */
  bbn_populate_free_data(grid_next);

  /* use previous grid to interpolate values of the fields for 
  // the next grid  and initialzing some other fields */
  interpolate_and_initialize_to_next_grid(grid_next,grid_prev);
  
  /* taking partial derivatives of the fields needed for equations */
  bbn_partial_derivatives_fields(grid_next);
  
  /* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S */
  bbn_update_stress_energy_tensor(grid_next,1);
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  bbn_update_Aij(grid_next);
  
  /* make normal vectorn on BH horizon */
  make_normal_vector_on_BH_horizon(grid_next);
  
  /* freeing */
  free_Grid_Params_S(GridParams);
  
  return grid_next;
}

/* keep NS center fixed at (0,C,0) by adjusting enthalpy */
static void keep_NS_center_fixed(Grid_T *const grid)
{
  const double C    = -0.5*GetParameterD_E("BH_NS_separation");
  struct NC_Center_RootFinder_S par[1] = {0};
  Root_Finder_T root_finder[1] = {0};
  const double x_center[3] = {0,C,0};
  
  par->patch = GetPatch("left_centeral_box",grid);
  par->root_finder = root_finder;
  
  /* print initial values before adjustments */
  printf("dh/dx before adjustment:\n");
  printf("dh/dx(%g,%g,%g)|NS center = %g\n",
    x_center[0],x_center[1],x_center[2],
    dh_dx0_root_finder_eq(par,x_center));
  printf("dh/dy(%g,%g,%g)|NS center = %g\n",
    x_center[0],x_center[1],x_center[2],
    dh_dx1_root_finder_eq(par,x_center));
  printf("dh/dz(%g,%g,%g)|NS center = %g\n",
    x_center[0],x_center[1],x_center[2],
    dh_dx2_root_finder_eq(par,x_center));
  
  if (strcmp_i(GetParameterS_E("NS_adjust_center_method"),"draw_enthalpy"))
  {
    adjust_NS_center_draw_enthalpy(grid);
  }
  else if (strcmp_i(GetParameterS_E("NS_adjust_center_method"),"tune_enthalpy"))
  {
    adjust_NS_center_tune_enthalpy(grid);
  }
  else
    abortEr(NO_OPTION);
  
  /* print initial values after adjustments */
  printf("dh/dx after adjustment:\n");
  printf("dh/dx(%g,%g,%g)|NS center = %g\n",
    x_center[0],x_center[1],x_center[2],
    dh_dx0_root_finder_eq(par,x_center));
  printf("dh/dy(%g,%g,%g)|NS center = %g\n",
    x_center[0],x_center[1],x_center[2],
    dh_dx1_root_finder_eq(par,x_center));
  printf("dh/dz(%g,%g,%g)|NS center = %g\n",
    x_center[0],x_center[1],x_center[2],
    dh_dx2_root_finder_eq(par,x_center));
}

/* controlling P_ADM */
static void P_ADM_control(Grid_T *const grid)
{
  char *adjust[3];
  const char *const par = GetParameterS_E("P_ADM_control_method");
  
  parse_adjust_parameter(par,adjust);
  
  void (*P_ADM_control_0)(Grid_T *const grid) =
                              get_func_P_ADM_adjustment(adjust[0]);
  
  void (*P_ADM_control_1)(Grid_T *const grid) =
                              get_func_P_ADM_adjustment(adjust[1]);
                              
  void (*P_ADM_control_2)(Grid_T *const grid) =
                              get_func_P_ADM_adjustment(adjust[2]);
  
  /* update P_ADM momentum parameters */
  Observable_T *obs = init_observable(grid);
  double p1[3] = {0};
  double p2[3] = {0};
  obs->quantity = "ADM_momentums";
  plan_observable(obs);
  
  /* get previous P_ADMs */
  p1[0] = GetParameterD_E("P_ADM_x");
  p1[1] = GetParameterD_E("P_ADM_y");
  p1[2] = GetParameterD_E("P_ADM_z");
  
  /* get the current P_ADMs */
  p2[0] = obs->Px_ADM(obs);
  p2[1] = obs->Py_ADM(obs);
  p2[2] = obs->Pz_ADM(obs);
  
  update_parameter_double_format("P_ADM_x",p2[0]);
  update_parameter_double_format("P_ADM_y",p2[1]);
  update_parameter_double_format("P_ADM_z",p2[2]);
  
  update_parameter_double_format("P_ADM_x_prev",p1[0]);
  update_parameter_double_format("P_ADM_y_prev",p1[1]);
  update_parameter_double_format("P_ADM_z_prev",p1[2]);
  
  free_observable(obs);
  
  if (P_ADM_control_0)
    P_ADM_control_0(grid);
  
  if (P_ADM_control_1)
    P_ADM_control_1(grid);
    
  if (P_ADM_control_2)
    P_ADM_control_2(grid);
    
  _free(adjust[0]);
  _free(adjust[1]);
  _free(adjust[2]);
}

/* getting adjustment str, returns the relevant function. */
fAdjustment_t *get_func_P_ADM_adjustment(const char *const adjust)
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
  else if (strcmp_i(adjust,"x_CM"))
  {
    f = Py_ADM_is0_by_x_CM;
  }
  else if (strcmp_i(adjust,"y_CM"))
  {
    f = Px_ADM_is0_by_y_CM;
  }
  else if (strcmp_i(adjust,"boost_Vx"))
  {
    f = Px_ADM_is0_by_x_boost;
  }
  else if (strcmp_i(adjust,"boost_Vy"))
  {
    f = Py_ADM_is0_by_y_boost;
  }
  else if (strcmp_i(adjust,"boost_Vz"))
  {
    f = Pz_ADM_is0_by_z_boost;
  }
  else
    abortEr(NO_OPTION);
  
  return f;
}

/* getting adjustment str, returns the relevant function. */
fAdjustment_t *get_func_force_balance_adjustment(const char *const adjust)
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
    abortEr(NO_OPTION);
  
  return f;
}

/* parsing adjust parameter value and fill components consequently */
static void parse_adjust_parameter(const char *const par,char *adjust[3])
{
  if (!strstr_i(par,"adjust(") && !strstr_i(par,"none"))
    abortEr_s("Syntax error for '%s'.\n",par);
  
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

/* adjust various quantities according to parameter:
// "force_balance_equation" */
static void force_balance_eq(Grid_T *const grid)
{
  char *adjust[3];
  const char *const par = GetParameterS_E("force_balance_equation");
  
  parse_adjust_parameter(par,adjust);
  
  void (*force_balance_0)(Grid_T *const grid) = 
            get_func_force_balance_adjustment(adjust[0]);
            
  void (*force_balance_1)(Grid_T *const grid) = 
            get_func_force_balance_adjustment(adjust[1]);
            
  void (*force_balance_2)(Grid_T *const grid) = 
            get_func_force_balance_adjustment(adjust[2]);
  
  
  if (force_balance_0)
    force_balance_0(grid);
  
  if (force_balance_1)
    force_balance_1(grid);
    
  if (force_balance_2)
    force_balance_2(grid);
    
  _free(adjust[0]);
  _free(adjust[1]);
  _free(adjust[2]);
}

/* adjust the boot velocity at the outer boundary to diminish P_ADM
// it only makes changes in the specified direction x. */
static void Px_ADM_is0_by_x_boost(Grid_T *const grid)
{
  const double SMALL_FAC = 1E-2;
  const double dP   = GetParameterD_E("P_ADM_control_tolerance");
  const double W    = GetParameterD_E("Solving_Field_Update_Weight");
  double p1[3] = {0};
  double p2[3] = {0};
  double v1[3] = {0};
  double v2[3] = {0};
  double v0[3] = {0};
  double  v[3] = {0};
  static unsigned iter = 0;
  
  /* get previous P_ADMs */
  p1[0] = GetParameterD_E("P_ADM_x_prev");
  
  /* get the current P_ADMs */
  p2[0] = GetParameterD_E("P_ADM_x");
  
  if (iter == 0)
  {
    update_parameter_double_format("v1_boost_x",p2[0]*SMALL_FAC);
  }
  else if (iter == 1)
  {
    update_parameter_double_format("v2_boost_x",p2[0]*SMALL_FAC);
  }
  else if (iter > 1)
  {
    /* get the penultimate boost velocity */
    v1[0] = GetParameterD_E("v1_boost_x");
    
    /* get the ultimate boost velocity */
    v2[0] = GetParameterD_E("v2_boost_x");
    
    /* get the boost velocity */
    v0[0] = GetParameterD_E("v*_boost_x");
    
    /* calculate the new boost velocity */
    v[0] = (v2[0]*p1[0]-v1[0]*p2[0])/(p1[0]-p2[0]);
    
    /* take cure of 0 denominator */
    if (EQL(p1[0],p2[0])) v[0] = 0;
    
    /* change the boost velocity relaxed */
    v[0] = W*v[0]+(1-W)*v0[0];
    
    /* update parameters */
    update_parameter_double_format("v1_boost_x",v2[0]);
    
    update_parameter_double_format("v2_boost_x",v[0]);
    
    const double dPx_Px = (p2[0]-p1[0])/fabs(p2[0]);
    printf("dPx/|Px| = %+e\n",dPx_Px);
    
    /* if change in momentum is big */
    if (GRT(fabs(dPx_Px),dP))
    {
      update_parameter_double_format("v*_boost_x",v[0]);
      printf("-->boost velocity_x = %e -> %e \n",v2[0],v[0]);
    }
    else
      printf("-->boost velocity_x = %e -> no update.\n",
        GetParameterD_E("v*_boost_x"));
      
  }
  
  iter++;
  UNUSED(grid);
} 

/* adjust the boot velocity at the outer boundary to diminish P_ADM
// it only makes changes in the specified direction y. */
static void Py_ADM_is0_by_y_boost(Grid_T *const grid)
{
  const double SMALL_FAC = 1E-2;
  const double dP   = GetParameterD_E("P_ADM_control_tolerance");
  const double W    = GetParameterD_E("Solving_Field_Update_Weight");
  double p1[3] = {0};
  double p2[3] = {0};
  double v1[3] = {0};
  double v2[3] = {0};
  double v0[3] = {0};
  double  v[3] = {0};
  static unsigned iter = 0;
  
  /* get previous P_ADMs */
  p1[1] = GetParameterD_E("P_ADM_y_prev");
  
  /* get the current P_ADMs */
  p2[1] = GetParameterD_E("P_ADM_y");
  
  if (iter == 0)
  {
    update_parameter_double_format("v1_boost_y",p2[1]*SMALL_FAC);
  }
  else if (iter == 1)
  {
    update_parameter_double_format("v2_boost_y",p2[1]*SMALL_FAC);
  }
  else if (iter > 1)
  {
    /* get the penultimate boost velocity */
    v1[1] = GetParameterD_E("v1_boost_y");
    
    /* get the ultimate boost velocity */
    v2[1] = GetParameterD_E("v2_boost_y");
    
    /* get the boost velocity */
    v0[1] = GetParameterD_E("v*_boost_y");
    
    /* calculate the new boost velocity */
    v[1] = (v2[1]*p1[1]-v1[1]*p2[1])/(p1[1]-p2[1]);
    
    /* take cure of 0 denominator */
    if (EQL(p1[1],p2[1])) v[1] = 0;
    
    /* change the boost velocity relaxed */
    v[1] = W*v[1]+(1-W)*v0[1];
    
    /* update parameters */
    update_parameter_double_format("v1_boost_y",v2[1]);
    
    update_parameter_double_format("v2_boost_y",v[1]);
    
    const double dPy_Py = (p2[1]-p1[1])/fabs(p2[1]);
    printf("dPy/|Py| = %+e\n",dPy_Py);
    
    if (GRT(fabs(dPy_Py),dP))
    {
      update_parameter_double_format("v*_boost_y",v[1]);
      printf("-->boost velocity_y = %e -> %e \n",v2[1],v[1]);
    }
    else
      printf("-->boost velocity_y = %e -> no update.\n",
        GetParameterD_E("v*_boost_y"));
      
  }
  
  iter++;
  UNUSED(grid);
} 

/* adjust the boot velocity at the outer boundary to diminish P_ADM
// it only makes changes in the specified direction z. */
static void Pz_ADM_is0_by_z_boost(Grid_T *const grid)
{
  const double SMALL_FAC = 1E-2;
  const double dP   = GetParameterD_E("P_ADM_control_tolerance");
  const double W    = GetParameterD_E("Solving_Field_Update_Weight");
  double p1[3] = {0};
  double p2[3] = {0};
  double v1[3] = {0};
  double v2[3] = {0};
  double v0[3] = {0};
  double  v[3] = {0};
  static unsigned iter = 0;
  
  /* get previous P_ADMs */
  p1[2] = GetParameterD_E("P_ADM_z_prev");
  
  /* get the current P_ADMs */
  p2[2] = GetParameterD_E("P_ADM_z");
  
  if (iter == 0)
  {
    update_parameter_double_format("v1_boost_z",p2[2]*SMALL_FAC);
  }
  else if (iter == 1)
  {
    update_parameter_double_format("v2_boost_z",p2[2]*SMALL_FAC);
  }
  else if (iter > 1)
  {
    /* get the penultimate boost velocity */
    v1[2] = GetParameterD_E("v1_boost_z");
    
    /* get the ultimate boost velocity */
    v2[2] = GetParameterD_E("v2_boost_z");
    
    /* get the boost velocity */
    v0[2] = GetParameterD_E("v*_boost_z");
    
    /* calculate the new boost velocity */
    v[2] = (v2[2]*p1[2]-v1[2]*p2[2])/(p1[2]-p2[2]);
    
    /* take cure of 0 denominator */
    if (EQL(p1[2],p2[2])) v[2] = 0;
    
    /* change the boost velocity relaxed */
    v[2] = W*v[2]+(1-W)*v0[2];
    
    /* update parameters */
    update_parameter_double_format("v1_boost_z",v2[2]);
    
    update_parameter_double_format("v2_boost_z",v[2]);
    
    const double dPz_Pz = (p2[2]-p1[2])/fabs(p2[2]);
    printf("dPz/|Pz| = %+e\n",dPz_Pz);
      
    if (GRT(fabs(dPz_Pz),dP))
    {
      update_parameter_double_format("v*_boost_z",v[2]);
      printf("-->boost velocity_z = %e -> %e \n",v2[2],v[2]);
    }
    else
      printf("-->boost velocity_z = %e -> no update.\n",
        GetParameterD_E("v*_boost_z"));
  }
  
  iter++;
  UNUSED(grid);
} 

/* adjust the center of NS at the designated point, in case it moved. 
// in this method, we tune enthalpy values such that the derivative of 
// the enthalpy be 0 at (0,NS_C,0). */
static void adjust_NS_center_tune_enthalpy(Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn    = patch->nn;
    unsigned ijk;
    
    if (!IsItNSPatch(patch) && !IsItNSSurroundingPatch(patch))
      continue;
    
    {/* local variables */
      PREP_FIELD(enthalpy)
      GET_FIELD(denthalpy_D2)
      GET_FIELD(denthalpy_D0)
    
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double x = patch->node[ijk]->x[0];
        double z = patch->node[ijk]->x[2];
        
        enthalpy[ijk] = enthalpy[ijk]*(1-denthalpy_D0[ijk]*x-denthalpy_D2[ijk]*z);
      }
    }
    {/* local variables */
      DECLARE_FIELD(enthalpy)

      DECLARE_AND_EMPTY_FIELD(denthalpy_D2)
      DECLARE_AND_EMPTY_FIELD(denthalpy_D1)
      DECLARE_AND_EMPTY_FIELD(denthalpy_D0)
      
      denthalpy_D2->v = Partial_Derivative(enthalpy,"z");
      denthalpy_D1->v = Partial_Derivative(enthalpy,"y");
      denthalpy_D0->v = Partial_Derivative(enthalpy,"x");
    }
  }
  
}

/* adjust the center of NS at the designated point, in case it moved.
// we need only to draw enthalpy to (0,NS_C,0).
// to do so, we demand shifted_enthalpy(r) = enthalpy(R+r), in which R 
// is the amount the center is displaced from (0,-D/2,0);
// and finally update the enthalpy and its derivatives. */
static void adjust_NS_center_draw_enthalpy(Grid_T *const grid)
{
  /* find the NS center */
  find_NS_center(grid);
  
  char par_name[1000];
  double *NS_center = 0;
  const double C    = -0.5*GetParameterD_E("BH_NS_separation");
  sprintf(par_name,"grid%u_NS_center_x",grid->gn);
  NS_center = GetParameterArrayF_E(par_name);
  const double R[3] = {NS_center[0],NS_center[1]-C,NS_center[2]};
  unsigned p,ijk;
  
  /* if it is already located at the designted point */
  if (EQL(0,R[0]) && EQL(0,R[1]) && EQL(0,R[2]))
    return;
    
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItNSPatch(patch) && !IsItNSSurroundingPatch(patch))
      continue;
    
    /* now shift enthalpy */
    Patch_T *patchp = 0;
    char *stem, hint[1000];
    DECLARE_FIELD(enthalpy);
    ADD_FIELD(shifted_enthalpy);
    DECLARE_FIELD(shifted_enthalpy);
    
    make_coeffs_3d(enthalpy);
    
    stem = strstr(patch->name,"_");/* the patch->name convention is grid\d?_root */
    assert(stem);
    stem++;
    sprintf(hint,"%s",stem);
    unsigned nn = patch->nn;
    double x[3],Xp[3];
    for (ijk = 0; ijk < nn; ++ijk)
    {
      x[0] = patch->node[ijk]->x[0]+R[0];
      x[1] = patch->node[ijk]->x[1]+R[1];
      x[2] = patch->node[ijk]->x[2]+R[2];
      find_X_and_patch(x,hint,grid,Xp,&patchp);
      
      /* if point x located outside of NS surrounding */
      if (LookUpField("enthalpy",patchp) < 0)
        shifted_enthalpy->v[ijk] = 1;
      else
        shifted_enthalpy->v[ijk] = interpolate_from_patch_prim("enthalpy",Xp,patchp);
    }
    
  }
  
  /* now clean enthalpy and copy new value and remove extras */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItNSPatch(patch) && !IsItNSSurroundingPatch(patch))
      continue;
    
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
    if(!IsItNSPatch(patch) && !IsItNSSurroundingPatch(patch))
      continue;
      
    DECLARE_FIELD(enthalpy)

    DECLARE_AND_EMPTY_FIELD(denthalpy_D2)
    DECLARE_AND_EMPTY_FIELD(denthalpy_D1)
    DECLARE_AND_EMPTY_FIELD(denthalpy_D0)
    
    denthalpy_D2->v = Partial_Derivative(enthalpy,"z");
    denthalpy_D1->v = Partial_Derivative(enthalpy,"y");
    denthalpy_D0->v = Partial_Derivative(enthalpy,"x");
  }
  
}

/* force_balance_equation : adjust x_CM at direction d/dx */
static void force_balance_ddx_x_CM(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,0,"x_CM");
}

/* force_balance_equation : adjust x_CM at direction d/dy */
static void force_balance_ddy_x_CM(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,1,"x_CM");
}

/* force_balance_equation : adjust x_CM at direction d/dz */
static void force_balance_ddz_x_CM(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,2,"x_CM");
}

/* force_balance_equation : adjust y_CM at direction d/dx */
static void force_balance_ddx_y_CM(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,0,"y_CM");
}

/* force_balance_equation : adjust y_CM at direction d/dy */
static void force_balance_ddy_y_CM(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,1,"y_CM");
}

/* force_balance_equation : adjust y_CM at direction d/dz */
static void force_balance_ddz_y_CM(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,2,"y_CM");
}

/* force_balance_equation : adjust Omega at direction d/dx */
static void force_balance_ddx_Omega(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,0,"BH_NS_orbital_angular_velocity");
}

/* force_balance_equation : adjust Omega at direction d/dy */
static void force_balance_ddy_Omega(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,1,"BH_NS_orbital_angular_velocity");
}

/* force_balance_equation : adjust Omega at direction d/dz */
static void force_balance_ddz_Omega(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,2,"BH_NS_orbital_angular_velocity");
}

/* find parameter par using force balance equation in direction dir */
static void force_balance_eq_root_finders(Grid_T *const grid,const int dir, const char *const par)
{
  const double D            = GetParameterD_E("BH_NS_separation");
  const double Vr           = GetParameterD_E("BH_NS_infall_velocity");
  const double NS_center[3] = {0,-D/2,0};/* since we keep the NS center always here */
  const double RESIDUAL     = sqrt(GetParameterD_E("RootFinder_Tolerance"));
  const double Omega_BHNS   = GetParameterD_E("BH_NS_orbital_angular_velocity");
  const double y_CM         = GetParameterD_E("y_CM");
  const double x_CM         = GetParameterD_E("x_CM");
  const double W1  = GetParameterD_E("Solving_Field_Update_Weight");
  const double W2  = 1-W1;
  double *new_par,old_par;
  double guess[1],X[3];
  struct Force_Balance_RootFinder_S params[1];
  Patch_T *patch;
  char desc[1000] = {'\0'};
  
  patch = GetPatch("left_centeral_box",grid);
  X_of_x(X,NS_center,patch);
  
  params->patch      = patch;
  params->X          = X;
  params->Vr         = Vr;
  params->D          = D;
  params->find_x_CM  = 0;
  params->find_y_CM  = 0;
  params->find_Omega = 0;
  params->dir        = dir;
  
  /* which par */
  if (strcmp_i("BH_NS_orbital_angular_velocity",par))
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
    abortEr(NO_OPTION);
  
  /* which direction */
  if (dir == 0)
  {
    params->dLnGamma = dLnGamma_in_force_balance_eq(patch,X,0);
  }
  else if (dir == 1)
  {
    params->dLnGamma = dLnGamma_in_force_balance_eq(patch,X,1);
  }
  else if (dir == 2)
  {
    params->dLnGamma = dLnGamma_in_force_balance_eq(patch,X,2);
  }
  else
    abortEr(NO_OPTION);
  
  sprintf(desc,"Solving Force Balance Eq. for '%s' at direction 'x^%d'",par,dir);
  
  Root_Finder_T *root = init_root_finder(1);
  root->description   = desc;
  root->verbose       = 1;
  root->type          = GetParameterS_E("RootFinder_Method");
  root->tolerance     = GetParameterD_E("RootFinder_Tolerance");
  root->MaxIter       = (unsigned)GetParameterI_E("RootFinder_Max_Number_of_Iteration");
  root->x_gss         = guess;
  root->params        = params;
  root->f[0]          = force_balance_root_finder_eq;
  plan_root_finder(root);
  
  new_par = execute_root_finder(root);
  
  if (root->exit_status != ROOT_FINDER_OK && GRT(root->residual,RESIDUAL))
  {
    print_root_finder_exit_status(root);
  }
  
  new_par[0] = W1*new_par[0]+W2*old_par;
  update_parameter_double_format(par,new_par[0]);
  printf("Updating %s: %g -> %g\n",par,old_par,new_par[0]);
  
  /* since B1 has been changed let's update the pertinent fields */
  update_B1_dB1_Beta_dBete_Aij_dAij(grid);
  
  free_root_finder(root);
  free(new_par);
}

/* find the NS center using d(enthalpy)/dx^i = 0 */
static void find_NS_center(Grid_T *const grid)
{
  double *NS_center;
  Root_Finder_T *root = init_root_finder(3);
  struct NC_Center_RootFinder_S params[1];
  const double RESIDUAL = sqrt(GetParameterD_E("RootFinder_Tolerance"));
  Flag_T success_f = NONE;
  double guess[3];/* initial guess for root finder */
  char par_name[1000];
  unsigned p;
  
  guess[0] = guess[2] = 0;
  guess[1] = -0.5*GetParameterD_E("BH_NS_separation");
  params->root_finder = root;
  root->description = "Finding NS center:";
  root->type        = GetParameterS_E("RootFinder_Method");
  root->tolerance   = GetParameterD_E("RootFinder_Tolerance");
  root->MaxIter     = (unsigned)GetParameterI_E("RootFinder_Max_Number_of_Iteration");
  root->x_gss       = guess;
  root->params      = params;
  root->verbose     = 1;
  root->f[0]        = dh_dx0_root_finder_eq;
  root->f[1]        = dh_dx1_root_finder_eq;
  root->f[2]        = dh_dx2_root_finder_eq;    
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItNSPatch(patch))
      continue;

    root->interrupt = 0;
    params->patch = patch;
    plan_root_finder(root);
    NS_center     = execute_root_finder(root);
    
    /* if root finder was successful */
    if (root->exit_status == ROOT_FINDER_OK || LSS(root->residual,RESIDUAL))
    {
      success_f = YES;
      /* save the position of NS center */
      sprintf(par_name,"grid%u_NS_center_x",grid->gn);
      add_parameter_array(par_name,NS_center,3);
      /* save the patch stem where the NS center takes place */
      sprintf(par_name,"grid%u_NS_center_patch",grid->gn);
      char *stem = strstr(patch->name,"_");
      assert(stem);
      stem++;
      add_parameter(par_name,stem);
      printf("NS center found at (%g,%g,%g).\n",NS_center[0],NS_center[1],NS_center[2]);
      free(NS_center);
      break;
    }
    free(NS_center);
  }
  if (success_f != YES)
  {
    print_root_finder_exit_status(root);
    abortEr("NS center could not be found.\n");
  }
  free_root_finder(root);
}

/* dh/dx^0 = 0 */
static double dh_dx0_root_finder_eq(void *params,const double *const x)
{
  struct NC_Center_RootFinder_S *const par = params;
  Patch_T *const patch = par->patch;
  DECLARE_FIELD(denthalpy_D0);
  Interpolation_T *interp_s;
  Needle_T *needle = alloc_needle();
  unsigned flg;
  double interp,X[3];
  
  needle->grid = patch->grid;
  needle->x    = x;
  needle_in(needle,patch);
  point_finder(needle);
  flg = needle->Nans;
  free_needle(needle);
  
  /* if this point is out of this patch, exit */
  if (flg == 0)
  {
    par->root_finder->interrupt = 1;
    return 0;
  }
  
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
  Patch_T *const patch = par->patch;
  DECLARE_FIELD(denthalpy_D1);
  Interpolation_T *interp_s;
  Needle_T *needle = alloc_needle();
  unsigned flg;
  double interp,X[3];
  
  needle->grid = patch->grid;
  needle->x    = x;
  needle_in(needle,patch);
  point_finder(needle);
  flg = needle->Nans;
  free_needle(needle);
  
  /* if this point is out of this patch, exit */
  if (flg == 0)
  {
    par->root_finder->interrupt = 1;
    return 0;
  }
  
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
  Patch_T *const patch = par->patch;
  DECLARE_FIELD(denthalpy_D2);
  Interpolation_T *interp_s;
  Needle_T *needle = alloc_needle();
  unsigned flg;
  double interp,X[3];
  
  needle->grid = patch->grid;
  needle->x    = x;
  needle_in(needle,patch);
  point_finder(needle);
  flg = needle->Nans;
  free_needle(needle);
  
  /* if this point is out of this patch, exit */
  if (flg == 0)
  {
    par->root_finder->interrupt = 1;
    return 0;
  }
  
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

/* find Euler equation constant to meet NS baryonic mass */
static void find_Euler_eq_const(Grid_T *const grid)
{
  Root_Finder_T *root = init_root_finder(1);
  double *Euler_const = 0;
  double guess[1];/* initial guess for Euler const */
  struct Euler_eq_const_RootFinder_S params[1];
  
  params->grid = grid;
  params->NS_baryonic_mass = GetParameterD_E("NS_baryonic_mass");
  guess[0] = GetParameterD_E("Euler_equation_constant");
  
  root->description = "Finding Euler equation constant using NS baryonic mass";
  root->type        = GetParameterS_E("RootFinder_Method");
  root->tolerance   = GetParameterD_E("RootFinder_Tolerance");
  root->MaxIter     = (unsigned)GetParameterI_E("RootFinder_Max_Number_of_Iteration");
  root->x_gss       = guess;
  root->params      = params;
  root->f[0]        = Euler_eq_const_rootfinder_eq;
  root->verbose     = 1;
  plan_root_finder(root);
  Euler_const       = execute_root_finder(root);
  update_parameter_double_format("Euler_equation_constant",Euler_const[0]);
  printf("Euler Equation const. updated: %g -> %g\n",guess[0],Euler_const[0]);
  free(Euler_const);
  free_root_finder(root);
}

/* find y_CM by demanding Px_ADM = 0 */
static void Px_ADM_is0_by_y_CM(Grid_T *const grid)
{
  double dy_CM = 0,px0,y_CM_new,p[3]={0};
  const double dP   = GetParameterD_E("P_ADM_control_tolerance");
  const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity");
  const double y_CM0 = GetParameterD_E("y_CM0");
  const double y_CM  = GetParameterD_E("y_CM");
  const double M_NS  = GetParameterD_E("NS_mass");
  const double M_BH  = GetParameterD_E("BH_mass");
  
  /* get P_ADM */
  p[0] = GetParameterD_E("P_ADM_x");
  p[1] = GetParameterD_E("P_ADM_y");
  p[2] = GetParameterD_E("P_ADM_z");
  px0  = GetParameterD_E("P_ADM_x_prev");
  
  printf("ADM momentums before y_CM update:\n");
  printf("P_ADM = (%e,%e,%e).\n",p[0],p[1],p[2]);
  
  /* changing center of mass */
  dy_CM    = -p[0]/(Omega_BHNS*(M_NS+M_BH));
  y_CM_new = y_CM0+dy_CM;
  
  const double dPx_Px = (px0-p[0])/fabs(p[0]);
  printf("dPx/|Px| = %+e\n",dPx_Px);
  
  /* having found new x_CM now update */
  if (GRT(fabs(dPx_Px),dP))
  {
    update_parameter_double_format("y_CM",y_CM_new);
    update_B1_dB1_Beta_dBete_Aij_dAij(grid);
    printf("Update Center of Rotation: %g -> %g.\n",y_CM,y_CM_new);
  }
  else
    printf("Update Center of Rotation: no update.\n");
    
}

/* find x_CM by demanding Py_ADM = 0 */
static void Py_ADM_is0_by_x_CM(Grid_T *const grid)
{
  double  dx_CM = 0,py0,x_CM_new,p[3]={0};
  const double dP   = GetParameterD_E("P_ADM_control_tolerance");
  const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity");
  const double x_CM0 = GetParameterD_E("x_CM0");
  const double x_CM  = GetParameterD_E("x_CM");
  const double M_NS  = GetParameterD_E("NS_mass");
  const double M_BH  = GetParameterD_E("BH_mass");
  
  /* get P_ADM */
  p[0] = GetParameterD_E("P_ADM_x");
  p[1] = GetParameterD_E("P_ADM_y");
  p[2] = GetParameterD_E("P_ADM_z");
  py0  = GetParameterD_E("P_ADM_y_prev");
  
  printf("ADM momentums before x_CM update:\n");
  printf("P_ADM = (%e,%e,%e).\n",p[0],p[1],p[2]);
  
  /* changing center of mass */
  dx_CM    = p[1]/(Omega_BHNS*(M_NS+M_BH));
  x_CM_new = x_CM0+dx_CM;
  
  const double dPy_Py = (py0-p[1])/fabs(p[1]);
  printf("dPy/|Py| = %+e\n",dPy_Py);
  
  /* having found new x_CM now update */
  if (GRT(fabs(dPy_Py),dP))
  {
    update_parameter_double_format("x_CM",x_CM_new);
    update_B1_dB1_Beta_dBete_Aij_dAij(grid);
    printf("Update Center of Rotation: %g -> %g.\n",x_CM,x_CM_new);
  }
  else
    printf("Update Center of Rotation: no update.\n");
    
}

/* update: (B1, dB1),(Beta,dBeta),(Aij and dAij)
// over the whole grid excluding inside of BH */
static void update_B1_dB1_Beta_dBete_Aij_dAij(Grid_T *const grid)
{
  const unsigned np = grid->np;
  unsigned p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    if (IsItInsideBHPatch(patch))
      continue;
    
    bbn_update_B1_U012(patch);
    bbn_update_derivative_B1_U0(patch);
    bbn_update_derivative_B1_U1(patch);
    bbn_update_derivative_B1_U2(patch);
    
    bbn_update_Beta_U0(patch);
    bbn_update_Beta_U1(patch);
    bbn_update_Beta_U2(patch);
    
    bbn_update_derivative_Beta_U0(patch);
    bbn_update_derivative_Beta_U1(patch);
    bbn_update_derivative_Beta_U2(patch);
    
    bbn_update_psi10A_UiUj(patch);
  }
}

/* find the apparent horizon radius to acquire the desired BH mass */
static void find_AH_radius(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  const double target_bh_mass  = GetParameterD_E("BH_mass");
  const double current_bh_mass = bbn_BH_Kommar_mass(grid);
  const double current_r_excision = GetParameterD_E("r_excision");
  const double W1  = GetParameterD_E("Solving_Field_Update_Weight");
  //const double W2  = 1-W1;
  const double dr  = -current_r_excision*(current_bh_mass/target_bh_mass-1);
  const double r_excision = current_r_excision + W1*dr;
  
  GridParams->R_BH_r = r_excision;
  update_parameter_double_format("r_excision",r_excision);
  
  printf("-->BH Kommar's mass = %e\n",current_bh_mass);
  printf("-->r_excision = %e -> %e \n",current_r_excision,r_excision);
}

/* find the Omega_BH to acquire the desired BH spin */
//static void find_BH_Omega(Grid_T *const grid,struct Grid_Params_S *const GridParams);

/* root finder eqution for Euler equation constant */
static double Euler_eq_const_rootfinder_eq(void *params,const double *const x)
{
  struct Euler_eq_const_RootFinder_S *const par = params;
  
  return bbn_NS_baryonic_mass(par->grid,x[0]) - par->NS_baryonic_mass;
}

/* use previous grid to interpolate values of the fields that will be solved for the next grid */
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev)
{
  const unsigned np = grid_next->np;
  unsigned p;
 
  /* the following fields are interpolated: */
  /* B0_U[0-2],psi,eta,phi,enthalpy */
  
  /* to avoid race condition between threads write all coeffs */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid_prev->np; ++p)/* note: grid_prev has more patches! */
  {
    Patch_T *patch = grid_prev->patch[p];
    
    DECLARE_FIELD(B0_U0)
    DECLARE_FIELD(B0_U1)
    DECLARE_FIELD(B0_U2)
    DECLARE_FIELD(psi)
    DECLARE_FIELD(eta)
    make_coeffs_3d(B0_U0);
    make_coeffs_3d(B0_U1);
    make_coeffs_3d(B0_U2);
    make_coeffs_3d(psi);
    make_coeffs_3d(eta);
    
    if (IsItNSPatch(patch))
    {
      DECLARE_FIELD(phi)
      DECLARE_FIELD(enthalpy)
      make_coeffs_3d(phi);
      make_coeffs_3d(enthalpy);
    }
  }
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < np; ++p)
  {
    Patch_T *patch = grid_next->patch[p];
    unsigned nn = patch->nn;
    char hint[100],*root_name;
    unsigned ijk;
    
    root_name = strstr(patch->name,"_");/* the patch->name convention is grid\d?_root */
    assert(root_name);
    root_name++;
    sprintf(hint,"%s",root_name);
    
    PREP_FIELD(B0_U0)
    PREP_FIELD(B0_U1)
    PREP_FIELD(B0_U2)
    PREP_FIELD(psi)
    PREP_FIELD(eta)
    
    if (IsItNSPatch(patch))
    {
      PREP_FIELD(phi)
      PREP_FIELD(enthalpy)
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double x[3] = { patch->node[ijk]->x[0],
                        patch->node[ijk]->x[1],
                        patch->node[ijk]->x[2] };
        double Xp[3] = {0};/* (X,Y,Z)(x,y,z) in grid_prev */
        Patch_T *patchp = 0;/* patch in grid_prev contains (x,y,z) */
        
        /* finding X and patch in grid_prev, associated to x */
        find_X_and_patch(x,hint,grid_prev,Xp,&patchp);
        
        B0_U0[ijk]    = interpolate_from_patch_prim("B0_U0",Xp,patchp);
        B0_U1[ijk]    = interpolate_from_patch_prim("B0_U1",Xp,patchp);
        B0_U2[ijk]    = interpolate_from_patch_prim("B0_U2",Xp,patchp);
        psi[ijk]      = interpolate_from_patch_prim("psi",Xp,patchp);
        eta[ijk]      = interpolate_from_patch_prim("eta",Xp,patchp);
        phi[ijk]      = interpolate_from_patch_prim("phi",Xp,patchp);
        enthalpy[ijk] = interpolate_from_patch_prim("enthalpy",Xp,patchp);
      }
    }
    else
    {
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double x[3] = { patch->node[ijk]->x[0],
                        patch->node[ijk]->x[1],
                        patch->node[ijk]->x[2] };
        double Xp[3] = {0};/* (X,Y,Z)(x,y,z) in grid_prev */
        Patch_T *patchp = 0;/* patch in grid_prev contains (x,y,z) */
        
        /* finding X and patch in grid_prev, associated to x */
        find_X_and_patch(x,hint,grid_prev,Xp,&patchp);
        
        B0_U0[ijk] = interpolate_from_patch_prim("B0_U0",Xp,patchp);
        B0_U1[ijk] = interpolate_from_patch_prim("B0_U1",Xp,patchp);
        B0_U2[ijk] = interpolate_from_patch_prim("B0_U2",Xp,patchp);
        psi[ijk]   = interpolate_from_patch_prim("psi",Xp,patchp);
        eta[ijk]   = interpolate_from_patch_prim("eta",Xp,patchp);
      }
    }
  }/* end of for (p = 0; p < np; ++p) */
  
  /* initializing some other fields: */
  /* W_U[0-2],Beta_U[0-2],B1_U[0-2] */
  const double Omega_NS_x = GetParameterD_E("NS_Omega_U0");
  const double Omega_NS_y = GetParameterD_E("NS_Omega_U1");
  const double Omega_NS_z = GetParameterD_E("NS_Omega_U2");
  const double C_NS = GetParameterD_E("NS_Center");
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < np; ++p)
  {
    Patch_T *patch = grid_next->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    bbn_update_B1_U012(patch);
    bbn_update_Beta_U0(patch);
    bbn_update_Beta_U1(patch);
    bbn_update_Beta_U2(patch);
    
    if (IsItNSPatch(patch))
    {
      PREP_FIELD(W_U0)
      PREP_FIELD(W_U1)
      PREP_FIELD(W_U2)
      
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double x = patch->node[ijk]->x[0];
        double y = patch->node[ijk]->x[1]-C_NS;
        double z = patch->node[ijk]->x[2];
        
        /* spin part */
        W_U0[ijk] = Omega_NS_y*z-Omega_NS_z*y;
        W_U1[ijk] = Omega_NS_z*x-Omega_NS_x*z;
        W_U2[ijk] = Omega_NS_x*y-Omega_NS_y*x;
      }
    }/* end of if (IsItNSPatch(patch)) */
    
  }/* end of for (p = 0; p < np; ++p) */
  
}

/* given field name, X and patch, finds the value of the field in X  
// using interpolation.
// ->return value: f(X) */
static double interpolate_from_patch_prim(const char *const field,const double *const X,Patch_T *const patch)
{
  double interp;
  Interpolation_T *interp_s = init_interpolation();
  
  interp_s->field = patch->pool[Ind(field)];
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}

/* given a cartesian point x on the grid, it finds the corresponding X and patch 
// on which this x takes place. 
// hint, is the name of the patch that potentially has the given x */
static void find_X_and_patch(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch)
{
  Interface_T **face;
  SubFace_T *subf;
  Needle_T *needle = alloc_needle();
  const double LOW_RES_ERR = 1E-9;
  unsigned *found;
  unsigned p,f,sf;
  
  needle->grid = grid;
  needle->x    = x;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!strstr(patch->name,hint))
    {
      continue;
    }
    else
    {
      needle_in(needle,patch);
      
      /* find all neighbors of this patch */
      face = patch->interface;
      /* for all faces */
      FOR_ALL(f,face)
      {
        /* for all subfaces */
        for (sf = 0; sf < face[f]->ns; ++sf)
        {
          subf = face[f]->subface[sf];
          if (subf->outerB || subf->innerB)
            continue;
            
          needle_in(needle,grid->patch[subf->adjPatch]);
        }
      }
      break;
    }/* end of else */
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  point_finder(needle);
  found = needle->ans;
  
  /* if it could not find X in neither the given hint patch nor its neighbors */
  if (needle->Nans)
  {
    *ppatch = grid->patch[found[0]];
    X_of_x(X,x,*ppatch);
  }
  else/* if no patch found let's find it in the other patches */
  {
    needle->ex   = needle->in;
    needle->Nex  = needle->Nin;
    needle->in   = 0;
    needle->Nin  = 0;
    
    point_finder(needle);
    found = needle->ans;
    
    /* if not found */
    if (!needle->Nans)
    {
      /* at the horizon of BH at low resulotion some points 
      // might not be found, let's fix this by hand. */
      if (strstr(hint,"right_BH"))
      {
        *ppatch = GetPatch(hint,grid);
        X_of_x(X,x,*ppatch);
        
        /* NOTE: This is assumed Cubed Spherical coords */
        if (LSS(fabs(X[2]),LOW_RES_ERR))
        {
          X[2] = 0;/* the culprit is at low resolution, 
                   // X[2] won't be found very close to 0! */
          /* make sure the X falls in the interval */
          assert(LSSEQL(X[0],1) && GRTEQL(X[0],-1));
          assert(LSSEQL(X[1],1) && GRTEQL(X[1],-1));
        }
        else
        {
          fprintf(stderr,"The point (%g,%g,%g) could not be found!\n",x[0],x[1],x[2]);
          abortEr("Point not found!\n");
        }
      }
      else
      {
        fprintf(stderr,"The point (%g,%g,%g) could not be found!\n",x[0],x[1],x[2]);
        abortEr("Point not found!\n");
      }
    }
    else
    {  
      *ppatch = grid->patch[found[0]];
      X_of_x(X,x,*ppatch);
    }
  }
  free_needle(needle);
}

#define ij(i,j) ((j)+Nphi*(i))
/* given the grid, find the NS surface on Ylm points 
// i.e (theta,phi) collocations are = (Legendre,EquiSpaced),
// using the fact that at the surface enthalpy = 1.
// it fills also NS radius attributes:
//   GridParams->Max_R_NS_l;
//   GridParams->NS_R_Ylm->realClm;
//   GridParams->NS_R_Ylm->imagClm;
//   GridParams->NS_R_Ylm->Lmax;
//
// we assumed cubed spherical grid */
static void find_NS_surface_Ylm_method_CS(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  pr_line_custom('=');
  printf("{ Finding the surface of NS, Ylm method ...\n");
  
  /* the stucture for the root finder */
  struct NS_surface_RootFinder_S par[1];
  unsigned Ntheta,Nphi;/* total number of theta and phi points */
  const unsigned lmax = (unsigned)GetParameterI_E("NS_surface_Ylm_expansion_max_l");
  const double RESIDUAL = sqrt(GetParameterD_E("RootFinder_Tolerance"));
  //const double W1  = GetParameterD_E("Solving_Field_Update_Weight");
  //const double W2  = 1-W1;
  double theta,phi;
  double *Rnew_NS = 0;/* new R for NS */
  double Max_R_NS = 0;/* maximum radius of NS */
  double guess    = 0;
  double X[3],x[3],N[3];
  char stem[1000],*affix;
  unsigned i,j;
  
  /* populate root finder */
  Root_Finder_T *root = init_root_finder(1);
  root->type      = GetParameterS_E("RootFinder_Method");
  root->tolerance = GetParameterD_E("RootFinder_Tolerance");
  root->MaxIter   = (unsigned)GetParameterI_E("RootFinder_Max_Number_of_Iteration");
  root->x_gss     = &guess;
  root->params    = par;
  root->f[0]      = bbn_NS_surface_enthalpy_eq;
  root->df_dx[0]  = bbn_NS_surface_denthalpy_dr;
  //root->verbose   = 1;
  //root->FD_Right  = 1;
  plan_root_finder(root);
  
  /* parameters for root finder */
  par->root_finder = root;
  
  /* initialize tables */
  init_Legendre_root_function();
  
  Ntheta  = Nphi = 2*lmax+1;
  Rnew_NS = alloc_double(Ntheta*Nphi);
  
  /* for each points of Ylm find the surface of NS */
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      
      Patch_T *h_patch = 0,*patch = 0;
      double y2[3] = {0};
      double h,*dr;
      
      /* find patch and X,Y,Z at NS surface in which theta and phi take place */
      find_XYZ_and_patch_of_theta_phi_NS_CS(X,&patch,theta,phi,grid);
      
      /* find enthalpy at the (X,Y,Z) */
      Interpolation_T *interp_h = init_interpolation();
      interp_h->field = patch->pool[Ind("enthalpy")];
      interp_h->XY_dir_flag  = 1;
      interp_h->X            = X[0];
      interp_h->Y            = X[1];
      interp_h->K            = patch->n[2]-1;
      plan_interpolation(interp_h);
      h = execute_interpolation(interp_h);/* enthalpy */
      
      //if (h <= 0)
        //printf("WARNING: enthalpy = %g\n",h);
      assert(h > 0);
      
      free_interpolation(interp_h);
      
      /* finding x */
      x_of_X(x,X,patch);
      
      /* r^ = sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^ */
      N[0]  = sin(theta)*cos(phi);
      N[1]  = sin(theta)*sin(phi);
      N[2]  = cos(theta);
      y2[0] = x[0]-patch->c[0];
      y2[1] = x[1]-patch->c[1];
      y2[2] = x[2]-patch->c[2];  
      
      if(LSSEQL(h,1))/* if it takes place at NS patch */
      {
        h_patch = patch;
      }
      else/* which means h = 1 occures in neighboring patch */
      {
        /* finding the juxtapose patch of this NS patch */
        affix = regex_find("_[[:alpha:]]{2,5}$",patch->name);
        assert(affix);
        sprintf(stem,"left_NS_surrounding%s",affix);
        free(affix);
        h_patch = GetPatch(stem,grid);
      }
      /* having found h_patch, now find the the position of h = 1 */
      par->x0[0] = x[0];
      par->x0[1] = x[1];
      par->x0[2] = x[2];
      par->patch = h_patch;
      par->N     = N;
      dr = execute_root_finder(root);
      /* if root finder is not OK for some reason */
      if (GRT(root->residual,RESIDUAL))
      {
        printf(". Root finder for NS surface at %s:\n.. ",h_patch->name);
        print_root_finder_exit_status(root);
        printf(".. Residual = %g\n",root->residual);
      }
      
      /*  new coords of R respect to the center of NS */
      y2[0] += N[0]*dr[0];
      y2[1] += N[1]*dr[0];
      y2[2] += N[2]*dr[0];
      Rnew_NS[ij(i,j)] = rms(3,y2,0);
      
      /* find the max NS radius */
      if (Rnew_NS[ij(i,j)] > Max_R_NS)
        Max_R_NS = Rnew_NS[ij(i,j)];
        
      free(dr);
        
    }/* end of for (j = 0; j < Nphi; ++j) */
  }/* end of for (i = 0; i < Ntheta; ++i) */
  
  /* adding maximum radius of NS to grid parameters */
  GridParams->Max_R_NS_l = Max_R_NS;
  
  /* making radius of NS parameter at each patch using Ylm interpolation */
  double *realClm = alloc_ClmYlm(lmax);
  double *imagClm = alloc_ClmYlm(lmax);
  
  /* calculating coeffs */
  get_Ylm_coeffs(realClm,imagClm,Rnew_NS,Ntheta,Nphi,lmax);
  GridParams->NS_R_Ylm->realClm = realClm;
  GridParams->NS_R_Ylm->imagClm = imagClm;
  GridParams->NS_R_Ylm->Lmax    = lmax;
  
  /* if some day you wanna filter Clm's */
  if (0)
  {
    const double e = 0.1;
    unsigned l,m;
    
    for (l = 0; l <= lmax; ++l)
      for (m = 0; m <= l; ++m)
      {
        unsigned lm = lm2n(l,m);
        realClm[lm] /= (1+e*SQR(l)*SQR(l+1));
        imagClm[lm] /= (1+e*SQR(l)*SQR(l+1));
      }
  }
  
  free(Rnew_NS);
  free_root_finder(root);
  
  printf("} Finding the surface of NS, Ylm method ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}
#ifdef ij
#undef ij
#endif

/* given theta, phi and knowing the fact that they are on NS surface, 
// it finds the corresponding patch and X,Y,Z coordinate. */
static void find_XYZ_and_patch_of_theta_phi_NS_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid)
{
  const double tan_phi    = tan(phi);
  const double cos_theta  = cos(theta);
  const double tan_phi2   = SQR(tan_phi);
  const double cos_theta2 = SQR(cos_theta);
  Flag_T found_flg = NO;
  unsigned p;
  
  X[2] = 1;/* since we are on NS surface */
  
  /* check all of NS patches in which (x,y,z) and 
  // (X,Y,Z) and (theta,phi) are consistent */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItNSSurface(patch))
      continue;

    Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
    const double *c = patch->c;
    double a,b;
    double a_sign,b_sign,c_sign;
    double x[3],phi2,theta2,r;
    
    /* we know that theta = 0 or Pi occures only at UP or DOWN patches
    // so don't bother to follow algorithm all the way down.
    // furthermore, this prevent 0 in denominator of unrelated patches. */
    if (EQL(theta,0) || EQL(theta,M_PI))
    {
      if (side == LEFT || side == RIGHT || 
          side == BACK || side == FRONT   )
        continue;
    }
    
    /* first calculate the magnetitude of a and b 
    // which are related to X[0] and X[1] with a sign */
    switch (side)
    {
      case UP:
        a = Sqrt((1 - cos_theta2)/(cos_theta2 + cos_theta2*tan_phi2));
        b = tan_phi*Sqrt((1 - cos_theta2)/(cos_theta2*(1 + tan_phi2)));
      break;
      case DOWN:
        b = Sqrt((1 - cos_theta2)/(cos_theta2 + cos_theta2*tan_phi2));
        a = tan_phi*Sqrt((1 - cos_theta2)/(cos_theta2*(1 + tan_phi2)));
      break;
      case LEFT:
        a = 1/tan_phi;
        b = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/((1. - cos_theta2)*tan_phi2));
      break;
      case RIGHT:
        b = 1/tan_phi;
        a = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/((1. - cos_theta2)*tan_phi2));
      break;
      case BACK:
        a = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/(1 - cos_theta2));
        b = tan_phi;
      break;
      case FRONT:
        b = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/(1 - cos_theta2));
        a = tan_phi;
      break;
      default:
        abortEr(NO_OPTION);
    }
    
    /* having found the magnitude of a and b, we need to find out the sign of them.
    // this is done by paying attention to side, signum(cos_theta) and range of tanphi */
    switch (side)
    {
      case UP:
        arctan_argument_signum(&b_sign,&a_sign,phi);
      break;
      case DOWN:
        arctan_argument_signum(&a_sign,&b_sign,phi);
      break;
      case LEFT:
        arctan_argument_signum(&c_sign,&a_sign,phi);
        if (cos_theta > 0) b_sign = 1;
        else		   b_sign = -1;
      break;
      case RIGHT:
        arctan_argument_signum(&c_sign,&b_sign,phi);
        if (cos_theta > 0) a_sign = 1;
        else		   a_sign = -1;
      break;
      case BACK:
        arctan_argument_signum(&b_sign,&c_sign,phi);
        if (cos_theta > 0) a_sign = 1;
        else		   a_sign = -1;
      break;
      case FRONT:
        arctan_argument_signum(&a_sign,&c_sign,phi);
        if (cos_theta > 0) b_sign = 1;
        else		   b_sign = -1;
      break;
      default:
        abortEr(NO_OPTION);
    }
    
    X[0] = fabs(a)*a_sign;
    X[1] = fabs(b)*b_sign;
    
    /* check if x of X really gives you the correct angles */
    x_of_X(x,X,patch);
    x[0] -= c[0];
    x[1] -= c[1];
    x[2] -= c[2];
    r = rms(3,x,0);
    theta2 = acos(x[2]/r);
    phi2   = arctan(x[1],x[0]);
    if (EQL(theta2,theta) && EQL(phi2,phi))
    {
      found_flg = YES;
      *ppatch = patch;
      break;
    }
  }
  if (found_flg == NO)
    abortEr("(X,Y,Z) or patch could not be found.\n");
}

/* make patches inside the excision region of BH and and extrapolate
// metric fields i.e. beta,eta and psi inside this region */
void bbn_extrapolate_metric_fields_insideBH(Grid_T *const grid)
{
  /* add patches in side the excision region */
  add_patches_insideBH(grid);
  
  /* extrapolate the fields inside the BH */
  extrapolate_insideBH(grid);
} 

/* add patches inside the excision region */
static void add_patches_insideBH(Grid_T *const grid)
{
  if (!strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
    abortEr(NO_OPTION);
    
  const unsigned np1 = grid->np;
  const unsigned np2 = np1+7;/* 6 cubed spherical + 1 box */
  Grid_T *bh_grid = 0;
  unsigned i;
  
  /* allocating */
  grid->patch = realloc(grid->patch,(np2+1)*sizeof(*grid->patch));
  pointerEr(grid->patch);
  grid->patch[np2] = 0;
  
  for (i = np1; i < np2; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    pointerEr(grid->patch[i]);
  }
  grid->np = np2;
  
  /* populate the patches */
  populate_right_BH_central_box(grid,np1);
  populate_right_BH(grid,np1+1);
  
  bh_grid = calloc(1,sizeof(*bh_grid));
  pointerEr(bh_grid);
  
  bh_grid->patch = &grid->patch[np1];
  alloc_nodes(bh_grid);
  make_nodes(bh_grid);
  
  /* test printing coords */
  if (test_print(PRINT_COORDS))
    pr_coords(grid);
    
  free(bh_grid);
}

/* extrapolate the fields B0,B1,eta and psi inside the BH.
// we assume they exponentially decreasing from BH horizon
// such that at center of BH the fields are all zero. 
// formula we use: f(r) = fi*exp(-atten(rh-r)), 
// fi = f(at bh horizon), rh = r at horizon, r = r inside the BH
// a = attenuation factor. */
static void extrapolate_insideBH(Grid_T *const grid)
{
  unsigned p;
  const double atten = 5E1;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if (!IsItInsideBHPatch(patch))
      continue;
     
    /* add fields: */
    ADD_FIELD(B0_U0)
    ADD_FIELD(B0_U1)
    ADD_FIELD(B0_U2)
    ADD_FIELD(B1_U0)
    ADD_FIELD(B1_U1)
    ADD_FIELD(B1_U2)
    ADD_FIELD(psi)
    ADD_FIELD(eta)
    
    GET_FIELD(B0_U0)
    GET_FIELD(B0_U1)
    GET_FIELD(B0_U2)
    GET_FIELD(psi)
    GET_FIELD(eta)
    
    bbn_update_B1_U012(patch);
    
    /* for the centeral box we demand the field values be 0. */
    if (strstr(patch->name,"right_centeral_box"))
      continue;
    
    Patch_T *BHsur_patch;/* corresponding bh surrounding patch */
    const unsigned *n = patch->n;
    double x[3],X[3];
    double rh,r,e;
    unsigned ijk,i,j,k;

    /* find the corresponding BH surrounding patch to be used for extrapolation */
    char stem[1000];
    char *affix = regex_find("_[[:alpha:]]{2,5}$",patch->name);/* finding the side of the patch */
    assert(affix);
    sprintf(stem,"right_BH_surrounding%s",affix);
    free(affix);
    BHsur_patch = GetPatch(stem,grid);
    
    /* prepare interpolation arguments */
    Interpolation_T *interp_eta = init_interpolation();
    interp_eta->field = BHsur_patch->pool[LookUpField_E("eta",BHsur_patch)];
    interp_eta->XYZ_dir_flag = 1;
    plan_interpolation(interp_eta);

    Interpolation_T *interp_psi = init_interpolation();
    interp_psi->field = BHsur_patch->pool[LookUpField_E("psi",BHsur_patch)];
    interp_psi->XYZ_dir_flag = 1;
    plan_interpolation(interp_psi);

    Interpolation_T *interp_B0_U0 = init_interpolation();
    interp_B0_U0->field = BHsur_patch->pool[LookUpField_E("B0_U0",BHsur_patch)];
    interp_B0_U0->XYZ_dir_flag = 1;
    plan_interpolation(interp_B0_U0);

    Interpolation_T *interp_B0_U1 = init_interpolation();
    interp_B0_U1->field = BHsur_patch->pool[LookUpField_E("B0_U1",BHsur_patch)];
    interp_B0_U1->XYZ_dir_flag = 1;
    plan_interpolation(interp_B0_U1);

    Interpolation_T *interp_B0_U2 = init_interpolation();
    interp_B0_U2->field = BHsur_patch->pool[LookUpField_E("B0_U2",BHsur_patch)];
    interp_B0_U2->XYZ_dir_flag = 1;
    plan_interpolation(interp_B0_U2);

    /* smearing the fields exponentially */
    for (i = 0; i < n[0]; ++i)
    {
      for (j = 0; j < n[1]; ++j)
      {
        /* calculate rh at hotizon patch */
        ijk = L(n,i,j,n[2]-1);
        x[0] = patch->node[ijk]->x[0]-patch->c[0];
        x[1] = patch->node[ijk]->x[1]-patch->c[1];
        x[2] = patch->node[ijk]->x[2]-patch->c[2]; 
        rh   = rms(3,x,0);
        
        X_of_x(X,patch->node[ijk]->x,BHsur_patch);

        interp_eta->X  = X[0];
        interp_eta->Y  = X[1];
        interp_eta->Z  = X[2];

        interp_psi->X  = X[0];
        interp_psi->Y  = X[1];
        interp_psi->Z  = X[2];

        interp_B0_U0->X = X[0];
        interp_B0_U0->Y = X[1];
        interp_B0_U0->Z = X[2];

        interp_B0_U1->X = X[0];
        interp_B0_U1->Y = X[1];
        interp_B0_U1->Z = X[2];

        interp_B0_U2->X = X[0];
        interp_B0_U2->Y = X[1];
        interp_B0_U2->Z = X[2];

        double eta_i   = execute_interpolation(interp_eta);
        double psi_i   = execute_interpolation(interp_psi);
        double B0_U0_i = execute_interpolation(interp_B0_U0);
        double B0_U1_i = execute_interpolation(interp_B0_U1);
        double B0_U2_i = execute_interpolation(interp_B0_U2);

        for (k = 0; k < n[2]; ++k)
        {
          ijk = L(n,i,j,k);
          x[0] = patch->node[ijk]->x[0]-patch->c[0];
          x[1] = patch->node[ijk]->x[1]-patch->c[1];
          x[2] = patch->node[ijk]->x[2]-patch->c[2];
          
          r = rms(3,x,0);
          e = exp(-atten*(rh-r));
          eta[ijk]   = eta_i*e;
          psi[ijk]   = psi_i*e;
          B0_U0[ijk] = B0_U0_i*e;
          B0_U1[ijk] = B0_U1_i*e;
          B0_U2[ijk] = B0_U2_i*e;
          
        }/* end of for (k = 0 ; k < n[2]; ++k) */
      }/* end of for (j = 0; j < n[1]; ++j) */
    }/* end of for (i = 0; i < n[0]; ++i) */
    free_interpolation(interp_eta);
    free_interpolation(interp_psi);
    free_interpolation(interp_B0_U0);
    free_interpolation(interp_B0_U1);
    free_interpolation(interp_B0_U2);
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
}

/* extrapolating phi, dphi and W in NS surrounding coords 
// in case they are needed for interpolation to the next grid
// or in calculation of enthalpy at NS surrounding patches.
// for extrapolation we demand:
// the fields spread out in the same fashion as it is changing toward 
// the NS sarface. what we have :
// f(r_out) = (f(r_in) + df)*exp(g(r)), in which df is f(r2)-f(r1), 
// r1 = FACTOR*r2 and r2 is the radius of NS surface, and g(r) is 
// a function of r_out/r2 to control the radial trend of fieldd.
// a = (r_max-r2)/(r2-r1)
// b = (r1*r_max-SQR(r2))/(r1-r2)
// r_out = a*r_in + b, where (r_in) r_out is r (inside)outside NS and
// r_max is the max radius in NS surrounding patch.
// in effect, it means we emulate the trend of the fields inside of NS
// from r1 to r2 in NS surrounding, from r2 to r_max. */
static void extrapolate_fluid_fields_outsideNS_CS(Grid_T *const grid)
{
  const double FACTOR = 0.8;/* r1 = FACTOR*r2 */
  //const double EXP    = 1;/* (r_out/r2)^EXP */
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    /* surrounding patch */
    Patch_T *patch = grid->patch[p];
    if (!IsItNSSurroundingPatch(patch))
      continue;
     
    /* add fields: */
    /* enthalpy */
    ADD_FIELD(enthalpy)
    ADD_FIELD_NoMem(denthalpy_D0)
    ADD_FIELD_NoMem(denthalpy_D1)
    ADD_FIELD_NoMem(denthalpy_D2)
    
    /* irrotational part of fluid */
    ADD_FIELD(phi)
    GET_FIELD(phi)
    ADD_FIELD_NoMem(dphi_D2)
    ADD_FIELD_NoMem(dphi_D1)
    ADD_FIELD_NoMem(dphi_D0)
    DECLARE_AND_EMPTY_FIELD(dphi_D2)
    DECLARE_AND_EMPTY_FIELD(dphi_D1)
    DECLARE_AND_EMPTY_FIELD(dphi_D0)
    
    /* spin part of fluid W^i */
    ADD_FIELD(W_U0)
    ADD_FIELD(W_U1)
    ADD_FIELD(W_U2)
    GET_FIELD(W_U0)
    GET_FIELD(W_U1)
    GET_FIELD(W_U2)
    
    Patch_T *NS_patch;/* corresponding NS patch, to extrapolate out */
    const unsigned *n = patch->n;
    double r1,r2,r_max,r_in,r_out,a,b;
    double phii,W_U0i,W_U1i,W_U2i;/* fields at r1 */
    double phif,W_U0f,W_U1f,W_U2f;/* fields at r2 */
    double dphi,dW_U0,dW_U1,dW_U2;
    //double sphi,sW_U0,sW_U1,sW_U2;/* signs */
    double x[3],X[3];
    double THETA,PHI;
    unsigned ijk,i,j,k;

    /* find the corresponding NS patch to be used for extrapolation */
    char stem[1000];
    char *affix = regex_find("_[[:alpha:]]{2,5}$",patch->name);/* finding the side of the patch */
    assert(affix);
    sprintf(stem,"left_NS%s",affix);
    free(affix);
    NS_patch = GetPatch(stem,grid);
    
    /* prepare interpolation arguments */
    Interpolation_T *interp_phi = init_interpolation();
    interp_phi->field = NS_patch->pool[LookUpField_E("phi",NS_patch)];
    interp_phi->XYZ_dir_flag = 1;
    plan_interpolation(interp_phi);

    Interpolation_T *interp_W_U0 = init_interpolation();
    interp_W_U0->field = NS_patch->pool[LookUpField_E("W_U0",NS_patch)];
    interp_W_U0->XYZ_dir_flag = 1;
    plan_interpolation(interp_W_U0);

    Interpolation_T *interp_W_U1 = init_interpolation();
    interp_W_U1->field = NS_patch->pool[LookUpField_E("W_U1",NS_patch)];
    interp_W_U1->XYZ_dir_flag = 1;
    plan_interpolation(interp_W_U1);

    Interpolation_T *interp_W_U2 = init_interpolation();
    interp_W_U2->field = NS_patch->pool[LookUpField_E("W_U2",NS_patch)];
    interp_W_U2->XYZ_dir_flag = 1;
    plan_interpolation(interp_W_U2);

    /* spreading out the fields value using interpolation */
    for (i = 0; i < n[0]; ++i)
    {
      for (j = 0; j < n[1]; ++j)
      {
        /* calculate r_max at NS surrounding patch */
        ijk = L(n,i,j,n[2]-1);
        x[0] = patch->node[ijk]->x[0]-patch->c[0];
        x[1] = patch->node[ijk]->x[1]-patch->c[1];
        x[2] = patch->node[ijk]->x[2]-patch->c[2]; 
        r_max = rms(3,x,0);

        /* calculate r2 at NS surface */
        ijk = L(n,i,j,0);
        x[0] = patch->node[ijk]->x[0]-patch->c[0];
        x[1] = patch->node[ijk]->x[1]-patch->c[1];
        x[2] = patch->node[ijk]->x[2]-patch->c[2];
        r2   = rms(3,x,0);

        /* calulate r1,a,b */
        r1 = FACTOR*r2;/* small scale to define r1 */
        a  = (r_max-r2)/(r2-r1);
        b  = (r1*r_max-SQR(r2))/(r1-r2);

        /* find the value of phi and W^{i} at r2 */
        THETA = acos(x[2]/r2);
        PHI   = arctan(x[1],x[0]);

        x[0] = r2*sin(THETA)*cos(PHI)+NS_patch->c[0];
        x[1] = r2*sin(THETA)*sin(PHI)+NS_patch->c[1];
        x[2] = r2*cos(THETA)+NS_patch->c[2];

        X_of_x(X,x,NS_patch);

        interp_phi->X  = X[0];
        interp_phi->Y  = X[1];
        interp_phi->Z  = X[2];

        interp_W_U0->X = X[0];
        interp_W_U0->Y = X[1];
        interp_W_U0->Z = X[2];

        interp_W_U1->X = X[0];
        interp_W_U1->Y = X[1];
        interp_W_U1->Z = X[2];

        interp_W_U2->X = X[0];
        interp_W_U2->Y = X[1];
        interp_W_U2->Z = X[2];

        phif  = execute_interpolation(interp_phi);
        W_U0f = execute_interpolation(interp_W_U0);
        W_U1f = execute_interpolation(interp_W_U1);
        W_U2f = execute_interpolation(interp_W_U2);

        /* find the value of phi and W^{i} at r1 */
        x[0] = r1*sin(THETA)*cos(PHI)+NS_patch->c[0];
        x[1] = r1*sin(THETA)*sin(PHI)+NS_patch->c[1];
        x[2] = r1*cos(THETA)+NS_patch->c[2];

        X_of_x(X,x,NS_patch);
        
        interp_phi->X  = X[0];
        interp_phi->Y  = X[1];
        interp_phi->Z  = X[2];

        interp_W_U0->X = X[0];
        interp_W_U0->Y = X[1];
        interp_W_U0->Z = X[2];

        interp_W_U1->X = X[0];
        interp_W_U1->Y = X[1];
        interp_W_U1->Z = X[2];

        interp_W_U2->X = X[0];
        interp_W_U2->Y = X[1];
        interp_W_U2->Z = X[2];

        phii  = execute_interpolation(interp_phi);
        W_U0i = execute_interpolation(interp_W_U0);
        W_U1i = execute_interpolation(interp_W_U1);
        W_U2i = execute_interpolation(interp_W_U2);

        dphi  = phif  - phii;
        dW_U0 = W_U0f - W_U0i;
        dW_U1 = W_U1f - W_U1i;
        dW_U2 = W_U2f - W_U2i;
        //sphi  = dphi  > 0 ?  -1 : 1;
        //sW_U0 = dW_U0 > 0 ?  -1 : 1;
        //sW_U1 = dW_U1 > 0 ?  -1 : 1;
        //sW_U2 = dW_U2 > 0 ?  -1 : 1;
        
        for (k = 0; k < n[2]; ++k)
        {
          ijk = L(n,i,j,k);
          x[0] = patch->node[ijk]->x[0]-patch->c[0];
          x[1] = patch->node[ijk]->x[1]-patch->c[1];
          x[2] = patch->node[ijk]->x[2]-patch->c[2];
          
          r_out = rms(3,x,0);
          r_in  = (r_out-b)/a;

          x[0] = r_in*sin(THETA)*cos(PHI)+NS_patch->c[0];
          x[1] = r_in*sin(THETA)*sin(PHI)+NS_patch->c[1];
          x[2] = r_in*cos(THETA)+NS_patch->c[2];

          /* find the value of phi and W^{i} at NS surface */
          X_of_x(X,x,NS_patch);

          interp_phi->X  = X[0];
          interp_phi->Y  = X[1];
          interp_phi->Z  = X[2];

          interp_W_U0->X = X[0];
          interp_W_U0->Y = X[1];
          interp_W_U0->Z = X[2];

          interp_W_U1->X = X[0];
          interp_W_U1->Y = X[1];
          interp_W_U1->Z = X[2];

          interp_W_U2->X = X[0];
          interp_W_U2->Y = X[1];
          interp_W_U2->Z = X[2];
          
          //double e0 = pow(r_out/r2,EXP);
          
          phi[ijk]  = execute_interpolation(interp_phi)+dphi;
          //phi[ijk] *= exp(sphi*e0)*pow(r_out/r2,sphi);
          W_U0[ijk] = execute_interpolation(interp_W_U0)+dW_U0;
          //W_U0[ijk] *= exp(sW_U0*e0)*pow(r_out/r2,sW_U0);
          W_U1[ijk] = execute_interpolation(interp_W_U1)+dW_U1;
          //W_U1[ijk] *= exp(sW_U1*e0)*pow(r_out/r2,sW_U1);
          W_U2[ijk] = execute_interpolation(interp_W_U2)+dW_U2;
          //W_U2[ijk] *= exp(sW_U2*e0)*pow(r_out/r2,sW_U2);
          
        }/* end of for (k = 0 ; k < n[2]; ++k) */
      }/* end of for (j = 0; j < n[1]; ++j) */
    }/* end of for (i = 0; i < n[0]; ++i) */
    free_interpolation(interp_phi);
    free_interpolation(interp_W_U0);
    free_interpolation(interp_W_U1);
    free_interpolation(interp_W_U2);

    Field_T *phi_field = patch->pool[Ind("phi")];
    dphi_D2->v = Partial_Derivative(phi_field,"z");
    dphi_D1->v = Partial_Derivative(phi_field,"y");
    dphi_D0->v = Partial_Derivative(phi_field,"x");
    
    /* populating enthalpy and its derivatives in NS surroundings */
    Tij_IF_CTS_enthalpy(patch);
    Field_T *enthalpy = patch->pool[Ind("enthalpy")];
    DECLARE_AND_EMPTY_FIELD(denthalpy_D2)
    DECLARE_AND_EMPTY_FIELD(denthalpy_D1)
    DECLARE_AND_EMPTY_FIELD(denthalpy_D0)
    denthalpy_D2->v = Partial_Derivative(enthalpy,"z");
    denthalpy_D1->v = Partial_Derivative(enthalpy,"y");
    denthalpy_D0->v = Partial_Derivative(enthalpy,"x");
    
  }/* end of FOR_ALL_PATCHES(p,grid) */
}

/* use TOV and Kerr-Schil black hole approximation.
// ->return value: resultant grid from this approximation */
static Grid_T *TOV_KerrSchild_approximation(void)
{
  Grid_T *grid = 0;
  struct Grid_Params_S *GridParams = init_GridParams();/* adjust some pars for construction of grid */
  
  /* solve fields for a TOV star located at left side of y axis */
  TOV_T *tov = TOV_init();
  tov->bar_m = GetParameterD_E("NS_baryonic_mass");
  tov->description = "Estimating NS";
  tov = TOV_solution(tov);
  const double ns_R = tov->rbar[tov->N-1];
  
  /* basics of Kerr Schild black hole located at right side of y axis */
  pr_line_custom('=');
  printf("Acquiring Black Hole properties ...\n");
  const double bh_chi  = GetParameterD_E("BH_X_U2");
  const double bh_mass = GetParameterD_E("BH_mass");
  const double bh_R    = bh_mass*(1+sqrt(1-SQR(bh_chi)));
  printf("BH properties:\n");
  printf("--> BH radius (Kerr-Schild Coords.) = %e\n",bh_R);
  printf("--> BH dimensionless spin (z comp.) = %e\n",bh_chi);
  printf("--> BH ADM mass                     = %e\n",bh_mass);
  printf("Acquiring Black Hole properties ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
 
  /* adding some parameters: */
  
  /* center of rotation (approx. Center of Mass) */
  const double D = GetParameterD_E("BH_NS_separation");
  const double C_BH = 0.5*D;/* center of BH it's on +y axis */
  const double C_NS = -C_BH;/* center of NS it's on -y axis*/
  const double ns_mass = tov->ADM_m;/* NS adm mass */
  const double y_CM = (ns_mass*C_NS+bh_mass*C_BH)/(ns_mass+bh_mass);
  add_parameter_double("y_CM",y_CM);
  add_parameter_double("x_CM",0);
  add_parameter_double("y_CM0",y_CM);
  add_parameter_double("x_CM0",0);
  add_parameter_double("NS_mass",ns_mass);
  add_parameter_double("NS_center",C_NS);
  add_parameter_double("r_excision",bh_R);
  
    
  /* -> BH_Omega, the angular frequency of the horizon,
  // is a free vector that determines the spin of BH
  // and it is related to the dimensionless spin by:
  // BH_X = 4*BH_mass*BH_Omega .
  // we only use U2 component, since we assume BH only has spin 
  // in +/- of z direction (PRD 86 084033) */
  AddParameterDoubleF("BH_Omega_U2",bh_chi/(4*bh_mass));

  /* -> the Constant of the integration of Euler equation */
  AddParameterDoubleF("Euler_equation_constant",0);
  
  /* -> centeral rho0 */
  AddParameterDoubleF("rho_center",1E-3);
  
  /* combining these two geometry to create the grid */
  GridParams->Max_R_NS_l = ns_R;
  GridParams->R_BH_r     = bh_R;
  GridParams->a_BH       = bh_chi*bh_mass;
  GridParams->NS_R_type  = "PerfectSphere";
  GridParams->BH_R_type  = "PerfectSphere";
  grid = creat_bbn_grid_CS(GridParams);
  
  /* creating all of the fields needed for construction of Initial Data */
  bbn_allocate_fields(grid);
  
  /* populating the free data part of initial data that we chose ourself */
  bbn_populate_free_data(grid);
  
  /* initialize the fields using TOV and Kerr-Schild solution */
  init_field_TOV_plus_KerrSchild(grid,tov,bh_chi*bh_mass,bh_mass);
  
  /* taking partial derivatives of the fields needed for equations */
  bbn_partial_derivatives_fields(grid);
  
  /* update u0, _J^i, _E and _S */
  Tij_IF_CTS_psi6Sources(grid);
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  bbn_update_Aij(grid);
  
  /* make normal vectorn on BH horizon */
  make_normal_vector_on_BH_horizon(grid);
  
  /* find Euler equation const using enthalpy of TOV star and other fields */
  find_Euler_eq_const_TOV_KerrSchild(grid);
  
  /* add some parameters for momentum and its adjustments */
  Observable_T *obs = init_observable(grid);
  obs->quantity     = "ADM_momentums";
  plan_observable(obs);
  double p_x,p_y,p_z;
  p_x = obs->Px_ADM(obs);
  p_y = obs->Py_ADM(obs);
  p_z = obs->Pz_ADM(obs);
  update_parameter_double_format("P_ADM_x",p_x);
  update_parameter_double_format("P_ADM_y",p_y);
  update_parameter_double_format("P_ADM_z",p_z);
  update_parameter_double_format("v*_boost_x",0);
  update_parameter_double_format("v*_boost_y",0);
  update_parameter_double_format("v*_boost_z",0);
  
  printf("ADM momentums initials:\n");
  printf("-->P_ADM = (%e,%e,%e).\n",p_x,p_y,p_z);
  
  /* freeing */
  free_Grid_Params_S(GridParams);
  TOV_free(tov);
  free_observable(obs);
  
  return grid;
}

/* find Euler equation const using enthalpy of TOV star and 
// other fields that are found in TOV Kerr-Schild approximation.
// this is important, since otherwise when h gets updated, it may
// become 'nan' due to incorrect value of Euler equation constant.
//
// FORMULA I USED: 
// C+h/u0+D_{i}phi*(D^{i}phi+W^{i})/(h*u0)-Beta^{i}*D_{i}phi = 0.
// 
// CPI SCRIPT that I used:
 
Declare = 
{

 # enthalpy
 (obj = Field,name = enthalpy, rank = 0, C_macro);

 # conformal metric inverse
 (obj = Field,name = _gammaI, rank = UU, C_macro);

 # spin part of fluid
 (obj = Field,name = W, rank = U, C_macro);

 # d(phi)/d? for irrotional part of fluid
 (obj = Field,name = dphi, rank = D, C_macro);

 # Beta
 (obj = Field,name = Beta, rank = U, C_macro);

 # conformal factor
 (obj = Field,name = psi, rank = 0, C_macro);

 # u0
 (obj = Field,name = u0, rank = 0, C_macro);
}
# symmetries:
Symm[_gammaI(i,j)  = _gammaI(j,i)];

psim4 = psi**(-4);
dphiP = psim4*_gammaI(i,j)*dphi(-i)*dphi(-j)+dphi(-i)*W(i);
Euler_C = -enthalpy/u0 - dphiP/(enthalpy*u0) + Beta(i)*dphi(-i);

*/
static void find_Euler_eq_const_TOV_KerrSchild(Grid_T *const grid)
{
  unsigned p,ijk;
  double Euler_C;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if(!IsItNSPatch(patch))
      continue;
    
    ijk = 0;/* for an arbitrary point */
    GET_FIELD(enthalpy)
    GET_FIELD(_gammaI_U0U2)
    GET_FIELD(_gammaI_U0U0)
    GET_FIELD(_gammaI_U0U1)
    GET_FIELD(_gammaI_U1U2)
    GET_FIELD(_gammaI_U1U1)
    GET_FIELD(_gammaI_U2U2)
    GET_FIELD(W_U1)
    GET_FIELD(W_U0)
    GET_FIELD(W_U2)
    GET_FIELD(dphi_D2)
    GET_FIELD(dphi_D1)
    GET_FIELD(dphi_D0)
    GET_FIELD(Beta_U1)
    GET_FIELD(Beta_U0)
    GET_FIELD(Beta_U2)
    GET_FIELD(psi)
    GET_FIELD(u0)

    double psim4 = 
pow(psi[ijk], -4);

    double dphiP = 
W_U0[ijk]*dphi_D0[ijk] + W_U1[ijk]*dphi_D1[ijk] + W_U2[ijk]*
dphi_D2[ijk] + psim4*(_gammaI_U0U0[ijk]*pow(dphi_D0[ijk], 2) + 2.0*
_gammaI_U0U1[ijk]*dphi_D0[ijk]*dphi_D1[ijk] + 2.0*_gammaI_U0U2[ijk]*
dphi_D0[ijk]*dphi_D2[ijk] + _gammaI_U1U1[ijk]*pow(dphi_D1[ijk], 2) +
2.0*_gammaI_U1U2[ijk]*dphi_D1[ijk]*dphi_D2[ijk] + _gammaI_U2U2[ijk]*
pow(dphi_D2[ijk], 2));

    Euler_C = 
(-dphiP - pow(enthalpy[ijk], 2) + enthalpy[ijk]*u0[ijk]*(Beta_U0[ijk]*
dphi_D0[ijk] + Beta_U1[ijk]*dphi_D1[ijk] + Beta_U2[ijk]*dphi_D2[ijk]))/
(enthalpy[ijk]*u0[ijk]);
    
    break;/* only for 1 patch we find Euler constant */
  }
  update_parameter_double_format("Euler_equation_constant",Euler_C);
}

/* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
// _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
static void bbn_update_Aij(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Updating _A^{ij}, _dA^{ij} and _A^{ij}*A_{ij} ...\n");
  unsigned p;

  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    bbn_update_psi10A_UiUj(patch);
  }
  
  printf("Updating _A^{ij}, _dA^{ij} and _A^{ij}*A_{ij} ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* make normal vectorn on BH horizon */
static void make_normal_vector_on_BH_horizon(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Making normal vector on BH horizon ...\n");
  
  unsigned p,nn,ijk;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItHorizonPatch(patch))
      continue;
      
    nn = patch->nn;
    
    GET_FIELD(_gamma_D2D2)
    GET_FIELD(_gamma_D0D2)
    GET_FIELD(_gamma_D0D0)
    GET_FIELD(_gamma_D0D1)
    GET_FIELD(_gamma_D1D2)
    GET_FIELD(_gamma_D1D1)
    
    /* normal vector on horizon */
    PREP_FIELD(_HS_U0);
    PREP_FIELD(_HS_U1);
    PREP_FIELD(_HS_U2);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      /* minus sign to point outside the black hole */
      _HS_U0[ijk] = dq2_dq1(patch,_c_,_x_,ijk);
      _HS_U1[ijk] = dq2_dq1(patch,_c_,_y_,ijk);
      _HS_U2[ijk] = dq2_dq1(patch,_c_,_z_,ijk);
      
      double N2 = 
pow(_HS_U0[ijk], 2)*_gamma_D0D0[ijk] + 2.0*_HS_U0[ijk]*_HS_U1[ijk]*
_gamma_D0D1[ijk] + 2.0*_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D2[ijk] +
pow(_HS_U1[ijk], 2)*_gamma_D1D1[ijk] + 2.0*_HS_U1[ijk]*_HS_U2[ijk]*
_gamma_D1D2[ijk] + pow(_HS_U2[ijk], 2)*_gamma_D2D2[ijk];
        
      double N = sqrt(N2);
      
      /* normalizing */
      _HS_U0[ijk] /= N;
      _HS_U1[ijk] /= N;
      _HS_U2[ijk] /= N;
      
    }
     
  }
  
  printf("Making normal vector on BH horizon ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
}

/* initialize the fields using TOV and Kerr-Schild solution.
// the idea is to superimpose two fields of each solutions. */
static void init_field_TOV_plus_KerrSchild(Grid_T *const grid,const TOV_T *const tov, const double a_BH, const double M_BH)
{
  pr_line_custom('=');
  printf("Initializing the fields using TOV and Kerr-Schild solution ...\n");
  
  Transformation_T *t = initialize_transformation();
  const double M_NS = tov->ADM_m;/* NS adm mass */
  const double D = GetParameterD_E("BH_NS_separation");
  const double C_BH = 0.5*D;/* center of BH it's on +y axis */
  const double C_NS = -C_BH;/* center of NS it's on -y axis*/
  const double R_Schwar = tov->r[tov->N-1];/* NS's Schwarzchild radius */
  const double a2_BH = SQR(a_BH);/* spin vector of BH */
  const double y_CM = GetParameterD_E("y_CM");
  const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity");
  const double Omega_NS_x = GetParameterD_E("NS_Omega_U0");
  const double Omega_NS_y = GetParameterD_E("NS_Omega_U1");
  const double Omega_NS_z = GetParameterD_E("NS_Omega_U2");
  double Bx,By,Bz;/* B = v/c */
  unsigned p;
  
  Bx = -Omega_BHNS*(C_BH-y_CM);
  By = 0;
  Bz = 0;
  t->boost->Bx = Bx;
  t->boost->By = By;
  t->boost->Bz = Bz;
  t->boost->B2 = SQR(Bx)+SQR(By)+SQR(Bz);

  /* black hole parts */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    PREP_FIELD(Beta_U0)
    PREP_FIELD(Beta_U1)
    PREP_FIELD(Beta_U2)
    PREP_FIELD(_gammaI_U0U2)
    PREP_FIELD(_gammaI_U0U0)
    PREP_FIELD(_gammaI_U0U1)
    PREP_FIELD(_gammaI_U1U2)
    PREP_FIELD(_gammaI_U1U1)
    PREP_FIELD(_gammaI_U2U2)

    ADD_FIELD(KSbeta_D0)
    ADD_FIELD(KSbeta_D1)
    ADD_FIELD(KSbeta_D2)
    GET_FIELD(KSbeta_D0)
    GET_FIELD(KSbeta_D1)
    GET_FIELD(KSbeta_D2)
    
    ADD_FIELD(KSalpha)
    GET_FIELD(KSalpha)
    
    /* beta and alpha needed */
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x = patch->node[ijk]->x[0];
      double y = patch->node[ijk]->x[1]-C_BH;
      double z = patch->node[ijk]->x[2];
      double x_mu[4] = {0/* time component */,x,y,z};/* x^mu in boost coords */
      double Lm1_x_mu[4];/* Lorentz^-1 x^mu, inverse boost */
      t->boost->inverse = 1;
      Lorentz_boost(t,x_mu,Lm1_x_mu);
      double _x    = Lm1_x_mu[1];
      double _y    = Lm1_x_mu[2];
      double _z    = Lm1_x_mu[3];
      double rbar  = bbn_KerrShcild_r(_x,_y,_z,a_BH);
      double rbar2 = SQR(rbar);
      double _k0 = (rbar*_x+a_BH*_y)/(rbar2+a2_BH);
      double _k1 = (rbar*_y-a_BH*_x)/(rbar2+a2_BH);
      double _k2 = _z/rbar;
      double _kt = 1;
      double _k_mu[4] = {_kt,_k0,_k1,_k2};
      double L_k_mu[4];/* Lorentz *k^mu */
      t->boost->inverse = 0;
      Lorentz_boost(t,_k_mu,L_k_mu);
      double kt = L_k_mu[0];
      double k0 = L_k_mu[1];
      double k1 = L_k_mu[2];
      double k2 = L_k_mu[3];
      double H  = bbn_KerrSchild_H(M_BH,rbar,a_BH,_z);
      double C = 2.*H;
      
      KSalpha[ijk] = 1/sqrt(1+C*kt*kt);
      KSbeta_D0[ijk]  = C*k0*kt;
      KSbeta_D1[ijk]  = C*k1*kt;
      KSbeta_D2[ijk]  = C*k2*kt;
      
      /* note the followings are multiplied by _gammaI, 
      // they need also multiplication by (psi)^-4 to make gammaI 
      // which will be done after psi is made */
      double shift_U0 = 
KSbeta_D0[ijk]*_gammaI_U0U0[ijk] + KSbeta_D1[ijk]*_gammaI_U0U1[ijk] + 
KSbeta_D2[ijk]*_gammaI_U0U2[ijk];

      double shift_U1 = 
KSbeta_D0[ijk]*_gammaI_U0U1[ijk] + KSbeta_D1[ijk]*_gammaI_U1U1[ijk] + 
KSbeta_D2[ijk]*_gammaI_U1U2[ijk];

      double shift_U2 = 
KSbeta_D0[ijk]*_gammaI_U0U2[ijk] + KSbeta_D1[ijk]*_gammaI_U1U2[ijk] + 
KSbeta_D2[ijk]*_gammaI_U2U2[ijk];


      /* populating: */
      Beta_U1[ijk] = shift_U1;
      Beta_U0[ijk] = shift_U0;
      Beta_U2[ijk] = shift_U2;
      
    }
    
  }/* end of black hole part */
  free_transformation(t);
  
  /* initialization psi, eta and matter fields: */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    PREP_FIELD(psi)
    PREP_FIELD(eta)
    PREP_FIELD(KSalpha)
    
    if (IsItNSPatch(patch))
    {
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

      EoS_T *eos = initialize_EoS();
      
      PREP_FIELD(enthalpy)
      PREP_FIELD(rho0)
      PREP_FIELD(phi)
      PREP_FIELD(W_U0)
      PREP_FIELD(W_U1)
      PREP_FIELD(W_U2)
      
      for (ijk = 0; ijk < nn; ++ijk)
      {
        /* note that we naturally using isotropic coords. 
        // for our coordiante, so bar in rbar is dropped */
        double x = patch->node[ijk]->x[0];
        double y = patch->node[ijk]->x[1]-C_NS;
        double z = patch->node[ijk]->x[2];
        double r = sqrt(SQR(x)+SQR(y)+SQR(z));
        double alpha;
        double enthalpy_h;
        
        interp_psi->N_cubic_spline_1d->h = r;
        interp_h->N_cubic_spline_1d->h = r;
        
        /* psi */
        psi[ijk] = execute_interpolation(interp_psi);
        /* + 1 for KS but we won't add and so we won't subtrac 1 either */
        
        /* eta */
        enthalpy_h = execute_interpolation(interp_h);
        alpha = sqrt(1-2*M_NS/R_Schwar)/enthalpy_h/* NS part */ + 
                KSalpha[ijk]/* BH part */ - 1./* supper-position */;
        eta[ijk] = psi[ijk]*alpha;
        
        /* enthalpy */
        enthalpy[ijk] = enthalpy_h;
        
        /* rho0 */
        eos->h = enthalpy_h;
        rho0[ijk] = eos->rest_mass_density(eos);
        
        /* phi Newtonian approximation */
        phi[ijk] = -Omega_BHNS*(C_NS-y_CM)*x;
        
        /* spin part */
        W_U0[ijk] = Omega_NS_y*z-Omega_NS_z*y;
        W_U1[ijk] = Omega_NS_z*x-Omega_NS_x*z;
        W_U2[ijk] = Omega_NS_x*y-Omega_NS_y*x;
      }
      free_interpolation(interp_psi);
      free_interpolation(interp_h);
      free_EoS(eos);
    }
    else/* outside NS */
    {
      for (ijk = 0; ijk < nn; ++ijk)
      {
        /* note that we naturally using isotropic for our coordiante, 
        // so bar in rbar is dropped */
        double x    = patch->node[ijk]->x[0];
        double y    = patch->node[ijk]->x[1]-C_NS;
        double z    = patch->node[ijk]->x[2];
        double r = sqrt(SQR(x)+SQR(y)+SQR(z));
        double alpha;
        
        /* psi */
        psi[ijk] = 1+0.5*M_NS/r;
        /* + 1 for KS but we won't add and so we won't subtrac 1 either */
        
        /* eta */
        alpha = (1-0.5*M_NS/r)/(1+0.5*M_NS/r)+KSalpha[ijk]-1;
        eta[ijk] = psi[ijk]*alpha;
        
      }
    }
    
  }/* end of initialization psi, eta and matter fields */
  
  /* initializing Beta and B's */
  FOR_ALL_PATCHES(p,grid)
  {
     Patch_T *patch = grid->patch[p];
     unsigned nn = patch->nn;
     unsigned ijk;
     double psim4;/* psi^-4 */
     
     bbn_update_B1_U012(patch);
      
     PREP_FIELD(B0_U0)
     PREP_FIELD(B0_U1)
     PREP_FIELD(B0_U2)
     
     PREP_FIELD(Beta_U0)
     PREP_FIELD(Beta_U1)
     PREP_FIELD(Beta_U2)
     
     GET_FIELD(B1_U0)
     GET_FIELD(B1_U1)
     GET_FIELD(B1_U2)
     GET_FIELD(psi)
     
     /* for outermost patches the better approximation is B0 = 0 */
     if (IsItOutermostPatch(patch))
     {
       for (ijk = 0; ijk < nn; ++ijk)
       {
         double x = patch->node[ijk]->x[0];
         double y = patch->node[ijk]->x[1];
         double z = patch->node[ijk]->x[2];
         double r = sqrt(SQR(x)+SQR(y)+SQR(z));
         psim4    = pow(psi[ijk],-4);
         
         /* Beta */
         Beta_U0[ijk] *= psim4;
         Beta_U1[ijk] *= psim4;
         Beta_U2[ijk] *= psim4;
         
         /* B0 */
         B0_U0[ijk] = Beta_U0[ijk]-B1_U0[ijk];
         B0_U1[ijk] = Beta_U1[ijk]-B1_U1[ijk];
         B0_U2[ijk] = Beta_U2[ijk]-B1_U2[ijk];
         B0_U0[ijk] /= r;
         B0_U1[ijk] /= r;
         B0_U2[ijk] /= r;
       }
     }
     else
     {
       for (ijk = 0; ijk < nn; ++ijk)
       {
         psim4 = pow(psi[ijk],-4);
         
         /* Beta */
         Beta_U0[ijk] *= psim4;
         Beta_U1[ijk] *= psim4;
         Beta_U2[ijk] *= psim4;
         
         /* B0 */
         B0_U0[ijk] = Beta_U0[ijk]-B1_U0[ijk];
         B0_U1[ijk] = Beta_U1[ijk]-B1_U1[ijk];
         B0_U2[ijk] = Beta_U2[ijk]-B1_U2[ijk];
         
       }
     }
  }/* end of * initializing Beta and B */
  
  /* freeing */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    DECLARE_FIELD(KSbeta_D0)
    DECLARE_FIELD(KSbeta_D1)
    DECLARE_FIELD(KSbeta_D2)
    DECLARE_FIELD(KSalpha)

    REMOVE_FIELD(KSbeta_D0)
    REMOVE_FIELD(KSbeta_D1)
    REMOVE_FIELD(KSbeta_D2)
    REMOVE_FIELD(KSalpha)
  }
  
  printf("Initializing the fields using TOV and Kerr-Schild solution ==> Done.\n");
  pr_clock();
  pr_line_custom('=');

}

/* given the max of radius of NS, radius of BH and their separation,
// and BH spin as the parameters, create a grid with these properties.
// NOTE: WE assume we are using cubed spherical grid.
// ->return value: grid of NS and BH in which inside of the BH excised. */
static Grid_T *creat_bbn_grid_CS(struct Grid_Params_S *const GridParams)
{
  Grid_T *grid = alloc_grid();/* adding a new grid */
  /* calculate the characteristics of this grid */
  const double Max_R_NS_l = GridParams->Max_R_NS_l;/* maximum radius of NS */
  const double R_BH_r     = GridParams->R_BH_r;
  const unsigned gn = grid->gn;
  const double C    = GetParameterD_E("BH_NS_separation");
  const unsigned N_Outermost_Split = (unsigned)GetParameterI_E("Number_of_Outermost_Split"); 
  double *R_outermost = alloc_double(N_Outermost_Split);
  double box_size_l,box_size_r;
  unsigned nlb[3]/*left box*/,n;
  char var[100] = {'\0'};
  char par[100] = {'\0'};
  char val[100] = {'\0'};
  const char *kind;
  unsigned i;
  
  /* finding the kind of grid */
  kind = GetParameterS_E("grid_kind");
  if (!strcmp_i(kind,"BBN_CubedSpherical_grid"))
    abortEr("This function only works with cubed spherical grid.\n");
    
  grid->kind = dup_s(kind);
  
  assert(GRT(C,0));
  assert(GRT(Max_R_NS_l,0));
  assert(GRT(R_BH_r,0));
  assert(LSS(2*Max_R_NS_l,C));
  assert(LSS(2*R_BH_r,C));
  
  /* making NS and BH surfaces function */
  NS_BH_surface_CubedSpherical_grid(grid,GridParams);
  
  box_size_l = GetParameterD_E("left_central_box_length_ratio")*Max_R_NS_l;
  box_size_r = GetParameterD_E("left_central_box_length_ratio")*R_BH_r;/* use same ratio as NS */
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R_outermost[i] = GetParameterD_E(var);
    
    if (LSS(R_outermost[i],2*C))
      abortEr("the radius of outermost patches must be greater than twice of BBN distance.");
    
    if (i > 0)
      if (LSSEQL(R_outermost[i],R_outermost[i-1]))
        abortEr("The radius of outermost must be increasing.");
    
  }
  
  /* filling n */
  
  /* left box */
  nlb[0] = (unsigned)GetParameterI("n_a");
  nlb[1] = (unsigned)GetParameterI("n_b");
  nlb[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("left_NS_n_a");
  if (n != INT_MAX)   nlb[0] = n;
  n = (unsigned)GetParameterI("left_NS_n_b");
  if (n != INT_MAX)   nlb[1] = n;
  n = (unsigned)GetParameterI("left_NS_n_c");
  if (n != INT_MAX)   nlb[2] = n;
    
  if(nlb[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(nlb[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(nlb[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* adding the results to the parameter data base */
  
  /* n_a, n_b, n_c */
  /* left box */
  sprintf(par,"grid%u_left_centeral_box_n_a",gn);
  sprintf(val,"%u",nlb[0]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_left_centeral_box_n_b",gn);
  sprintf(val,"%u",nlb[1]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_left_centeral_box_n_c",gn);
  sprintf(val,"%u",nlb[2]);
  add_parameter_string(par,val);
  
  /* size a,b,c */
  sprintf(par,"grid%u_left_centeral_box_size_a",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_b",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_c",gn);
  add_parameter_double(par,box_size_l);
  
  /* surrounding box length */
  sprintf(par,"grid%u_surrounding_box_length",gn);
  add_parameter_double(par,C);
  
  /* right box. NOTE: this is needed when we fill the excision region */
  nlb[0] = (unsigned)GetParameterI("n_a");
  nlb[1] = (unsigned)GetParameterI("n_b");
  nlb[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("right_BH_n_a");
  if (n != INT_MAX)   nlb[0] = n;
  n = (unsigned)GetParameterI("right_BH_n_b");
  if (n != INT_MAX)   nlb[1] = n;
  n = (unsigned)GetParameterI("right_BH_n_c");
  if (n != INT_MAX)   nlb[2] = n;
    
  if(nlb[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(nlb[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(nlb[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* adding the results to the parameter data base */
  
  /* n_a, n_b, n_c */
  /* right box */
  sprintf(par,"grid%u_right_centeral_box_n_a",gn);
  sprintf(val,"%u",nlb[0]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_centeral_box_n_b",gn);
  sprintf(val,"%u",nlb[1]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_centeral_box_n_c",gn);
  sprintf(val,"%u",nlb[2]);
  add_parameter_string(par,val);
  
  /* size a,b,c we take it the same as left box. no biggie! */
  sprintf(par,"grid%u_right_centeral_box_size_a",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_right_centeral_box_size_b",gn);
  add_parameter_double(par,box_size_r);
  
  sprintf(par,"grid%u_right_centeral_box_size_c",gn);
  add_parameter_double(par,box_size_r);
  
  /* R1 and R2 outermost */
  sprintf(par,"grid%u_outermost%u_R2",gn,0);
  add_parameter_double(par,R_outermost[0]);
    
  for (i = 1; i < N_Outermost_Split; i++)
  {
    /* R1: */
    sprintf(par,"grid%u_outermost%u_R1",gn,i);
    add_parameter_double(par,R_outermost[i-1]);
    
    /* R2: */
    sprintf(par,"grid%u_outermost%u_R2",gn,i);
    add_parameter_double(par,R_outermost[i]);
    
  }
  
  /* assuming the center of left NS at (0,-C/2,0) */
  sprintf(par,"grid%u_left_NS_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_left_NS_center_b",gn);
  add_parameter_double(par,-C/2);
  
  sprintf(par,"grid%u_left_NS_center_c",gn);
  add_parameter_double(par,0.0);
  
  /* assuming the center of right BH at (0,C/2,0) */
  sprintf(par,"grid%u_right_BH_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_right_BH_center_b",gn);
  add_parameter_double(par,C/2);
  
  sprintf(par,"grid%u_right_BH_center_c",gn);
  add_parameter_double(par,0.0);
  
  free(R_outermost);

  make_patches(grid);/* making patch(es) to cover the grid */
  realize_geometry(grid);/* realizing the geometry of whole grid
                     // including the way patches have been sewed,
                     // normal to the boundary, 
                     // outer-boundary, inner boundary and etc. */
  
  return grid;
  
}

/* making  NS and BH surfaces function */
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  const double Max_R_NS_l = GridParams->Max_R_NS_l;/* maximum radius of NS */
  const double R_BH_r     = GridParams->R_BH_r;
  const double a_BH       = GridParams->a_BH;
  const double y_CM       = GetParameterD_E("y_CM");
  const double C_BH       = 0.5*GetParameterD_E("BH_NS_separation");/* center of BH it's on +y axis */
  const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity");
  const double g2         = 1-SQR(-Omega_BHNS*(C_BH-y_CM));/* inverse square of Lorentz factor  */
  double *R;
  char par[1000] = {'\0'};
  unsigned N[3],n,i,j,k,N_total;
  Patch_T patch[1] = {0};
  struct Collocation_s coll_s[2] = {0};
  double X[3],r;
  
  /* left NS */
  
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("left_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("left_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("left_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface */
  R = alloc_double(N_total);
  /* if NS is perfect sphere like TOV*/
  if (strcmp_i(GridParams->NS_R_type,"PerfectSphere"))
  {
    for (i = 0; i < N[0]; ++i)
      for (j = 0; j < N[1]; ++j)
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = Max_R_NS_l;
    
    sprintf(par,"grid%u_left_NS_surface_function_up",grid->gn);
    add_parameter_array(par,R,N_total);
    sprintf(par,"grid%u_left_NS_surface_function_down",grid->gn);
    add_parameter_array(par,R,N_total);
    sprintf(par,"grid%u_left_NS_surface_function_back",grid->gn);
    add_parameter_array(par,R,N_total);
    sprintf(par,"grid%u_left_NS_surface_function_front",grid->gn);
    add_parameter_array(par,R,N_total);
    sprintf(par,"grid%u_left_NS_surface_function_left",grid->gn);
    add_parameter_array(par,R,N_total);
    sprintf(par,"grid%u_left_NS_surface_function_right",grid->gn);
    add_parameter_array(par,R,N_total);
  }
  /* if NS radius is varied and we know its expansion in Ylm bases */
  else if (strcmp_i(GridParams->NS_R_type,"SphericalHarmonic"))
  {
    /* we need interpolation */
    const double *const realClm = GridParams->NS_R_Ylm->realClm;
    const double *const imagClm = GridParams->NS_R_Ylm->imagClm;
    const unsigned Lmax   = GridParams->NS_R_Ylm->Lmax;
    double theta,phi;
    
    X[3] = 1;/* since we are on the NS surface */
    
    /* filling min */
    patch->min[0] = -1;
    patch->min[1] = -1;
    patch->min[2] = 0;

    /* filling max */
    patch->max[0] = 1;
    patch->max[1] = 1;
    patch->max[2] = 1;
    
    /* collocation */
    patch->collocation[0] = Chebyshev_Extrema;
    patch->collocation[1] = Chebyshev_Extrema;
    patch->collocation[2] = Chebyshev_Extrema;

    /* basis */
    patch->basis[0] = Chebyshev_Tn_BASIS;
    patch->basis[1] = Chebyshev_Tn_BASIS;
    patch->basis[2] = Chebyshev_Tn_BASIS;
      
    patch->n[0] = N[0];
    patch->n[1] = N[1];
    patch->n[2] = N[2];
  
    initialize_collocation_struct(patch,&coll_s[0],0);
    initialize_collocation_struct(patch,&coll_s[1],1);
    
    /* surface up */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,UP);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_up",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface down */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,DOWN);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_down",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface back */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,BACK);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_back",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface front */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,FRONT);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_front",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface left */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,LEFT);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_left",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface right */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,RIGHT);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_right",grid->gn);
    add_parameter_array(par,R,N_total);
  }
  else
    abortEr(NO_OPTION);
  
  free(R);
  
  /* right BH: */
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = 0;

  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;

  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  
  /* check for override */
  n = (unsigned)GetParameterI("right_BH_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("right_BH_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("right_BH_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
   
  patch->n[0] = N[0];
  patch->n[1] = N[1];
  patch->n[2] = N[2];
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  
  N_total = N[0]*N[1]*N[2];
  
  R = alloc_double(N_total);
  
  /* Boosted_Kerr-Schild radius */
  if (strcmp_i(GridParams->BH_R_type,"Boosted_Kerr-Schild"))
  {
    /* surface up */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        r = sqrt(
                 (1+SQR(X[0])+SQR(X[1]))/
                 ((g2*SQR(X[0])+SQR(X[1]))/(SQR(R_BH_r)+SQR(a_BH)) + 1/SQR(R_BH_r))
                );
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_up",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface down */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        r = sqrt(
                 (1+SQR(X[0])+SQR(X[1]))/
                 ((SQR(X[0])+g2*SQR(X[1]))/(SQR(R_BH_r)+SQR(a_BH)) + 1/SQR(R_BH_r))
                );
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_down",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface back */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);/* a = z/x */
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);/* b = y/x */
        r = sqrt(
                 (1+SQR(X[0])+SQR(X[1]))/
                 (((g2+SQR(X[1])))/(SQR(R_BH_r)+SQR(a_BH)) + SQR(X[0])/SQR(R_BH_r))
                );
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_back",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface front */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);/* a = y/x */
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);/* b = z/x */
        r = sqrt(
                 (1+SQR(X[0])+SQR(X[1]))/
                 (((g2+SQR(X[0])))/(SQR(R_BH_r)+SQR(a_BH)) + SQR(X[1])/SQR(R_BH_r))
                );
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_front",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface left */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);/* a = x/y */
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);/* b = z/y */
        r = sqrt(
                 (1+SQR(X[0])+SQR(X[1]))/
                 (((1+g2*SQR(X[0])))/(SQR(R_BH_r)+SQR(a_BH)) + SQR(X[1])/SQR(R_BH_r))
                );
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_left",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface right */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);/* a = z/y */
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);/* b = x/y */
        r = sqrt(
                 (1+SQR(X[0])+SQR(X[1]))/
                 (((1+g2*SQR(X[1])))/(SQR(R_BH_r)+SQR(a_BH)) + SQR(X[0])/SQR(R_BH_r))
                );
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_right",grid->gn);
    add_parameter_array(par,R,N_total);
  }
  else if (strcmp_i(GridParams->BH_R_type,"PerfectSphere"))
  {
    for (i = 0; i < N[0]; ++i)
      for (j = 0; j < N[1]; ++j)
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = R_BH_r;
    
    sprintf(par,"grid%u_right_BH_surface_function_up",grid->gn);
    add_parameter_array(par,R,N_total);
    sprintf(par,"grid%u_right_BH_surface_function_down",grid->gn);
    add_parameter_array(par,R,N_total);
    sprintf(par,"grid%u_right_BH_surface_function_back",grid->gn);
    add_parameter_array(par,R,N_total);
    sprintf(par,"grid%u_right_BH_surface_function_front",grid->gn);
    add_parameter_array(par,R,N_total);
    sprintf(par,"grid%u_right_BH_surface_function_left",grid->gn);
    add_parameter_array(par,R,N_total);
    sprintf(par,"grid%u_right_BH_surface_function_right",grid->gn);
    add_parameter_array(par,R,N_total);
  }
  else
    abortEr(NO_OPTION);
  
  free(R);
}

/* given (X,Y,Z) in the specified slice of NS in cubed spherical coords
// it finds the associated polar and azimuthal angels on the surface of NS */
static void find_theta_phi_of_XYZ_NS_CS(double *const theta,double *const phi,const double *const X,const Flag_T side)
{
  const double a = X[0];
  const double b = X[1];
  const double d = sqrt(1+SQR(a)+SQR(b));
  
  switch (side)
  {
    case UP:
      *phi   = arctan(b,a);
      *theta = acos(1/d);
    break;
    case DOWN:
      *phi   = arctan(a,b);
      *theta = acos(-1/d);
    break;
    case LEFT:
      *phi   = arctan(-1,a);
      *theta = acos(b/d);
    break;
    case RIGHT:
      *phi   = arctan(1,b);
      *theta = acos(a/d);
    break;
    case BACK:
      *phi   = arctan(b,-1);
      *theta = acos(a/d);
    break;
    case FRONT:
      *phi   = arctan(a,1);
      *theta = acos(b/d);
    break;
    default:
      abortEr(NO_OPTION);
  }
  
}

/* initialize Grid_Params struct */
static struct Grid_Params_S *init_GridParams(void)
{
  struct Grid_Params_S *par = calloc(1,sizeof(*par));
  return par;
}

/* free Grid_Params struct */
static void free_Grid_Params_S(struct Grid_Params_S *par)
{
  _free(par->NS_R_Ylm->realClm);
  _free(par->NS_R_Ylm->imagClm);
  _free(par);
}

/* find r such that f(h(r)) = h(r)-1 = 0.
// the root finder moving along the r and it seeks for h = 1. */
static double bbn_NS_surface_enthalpy_eq(void *params,const double *const x)
{
  const struct NS_surface_RootFinder_S *const pars = params;
  Patch_T *const patch0  = pars->patch;/* supposed to find in this patch */
  const double dx        = x[0];
  const double *const x0 = pars->x0;
  const double *const N  = pars->N;
  const double y[3]      = {x0[0]+dx*N[0],x0[1]+dx*N[1],x0[2]+dx*N[2]};
  Patch_T *patch = 0;
  double X[3],h;
  int h_ind;
  char hint[1000],*root_name;
  
  root_name = strstr(patch0->name,"_");/* the patch->name convention is grid\d?_root */
  assert(root_name);
  root_name++;
  sprintf(hint,"%s",root_name);
    
  find_X_and_patch(y,hint,patch0->grid,X,&patch);
  
  /* find enthalpy at the (X,Y,Z) */
  h_ind = _Ind("enthalpy");
  if (h_ind < 0)/* if there is no enthalpy defined in the patch */
    return -1;
    
  Interpolation_T *interp_h = init_interpolation();
  interp_h->field = patch->pool[h_ind];
  interp_h->XYZ_dir_flag  = 1;
  interp_h->X            = X[0];
  interp_h->Y            = X[1];
  interp_h->Z            = X[2];
  plan_interpolation(interp_h);
  h = execute_interpolation(interp_h);/* enthalpy */
  free_interpolation(interp_h);
  
  return h-1;
}

/* denthalpy(r)/dr for NS surface root finder */
static double bbn_NS_surface_denthalpy_dr(void *params,const double *const x,const unsigned dir)
{
  assert(dir == 0);
  const struct NS_surface_RootFinder_S *const pars = params;
  Patch_T *const patch0  = pars->patch;/* supposed to find in this patch */
  const double dx        = x[0];
  const double *const x0 = pars->x0;
  const double *const N  = pars->N;
  const double y[3]      = {x0[0]+dx*N[0],x0[1]+dx*N[1],x0[2]+dx*N[2]};
  Patch_T *patch = 0;
  double X[3],dh_dx,dh_dy,dh_dz;
  int dh_dx_ind,dh_dy_ind,dh_dz_ind;
  Interpolation_T *interp_dh_dx = 0,
                  *interp_dh_dy = 0,
                  *interp_dh_dz = 0;
  char hint[1000],*root_name;
  
  root_name = strstr(patch0->name,"_");/* the patch->name convention is grid\d?_root */
  assert(root_name);
  root_name++;
  sprintf(hint,"%s",root_name);
    
  find_X_and_patch(y,hint,patch0->grid,X,&patch);
  
  /* find denthalpy/dr at the (X,Y,Z): */
  
  dh_dx_ind = _Ind("denthalpy_D0");
  if (dh_dx_ind < 0)/* if there is no enthalpy defined in the patch */
    return 1;
  interp_dh_dx = init_interpolation();
  interp_dh_dx->field = patch->pool[dh_dx_ind];
  interp_dh_dx->X = X[0];
  interp_dh_dx->Y = X[1];
  interp_dh_dx->Z = X[2];
  interp_dh_dx->XYZ_dir_flag = 1;
  plan_interpolation(interp_dh_dx);
  dh_dx = execute_interpolation(interp_dh_dx);
  free_interpolation(interp_dh_dx);
  
  dh_dy_ind = _Ind("denthalpy_D1");
  if (dh_dy_ind < 0)/* if there is no enthalpy defined in the patch */
    return 1;
  interp_dh_dy = init_interpolation();
  interp_dh_dy->field = patch->pool[dh_dy_ind];
  interp_dh_dy->X = X[0];
  interp_dh_dy->Y = X[1];
  interp_dh_dy->Z = X[2];
  interp_dh_dy->XYZ_dir_flag = 1;
  plan_interpolation(interp_dh_dy);
  dh_dy = execute_interpolation(interp_dh_dy);
  free_interpolation(interp_dh_dy);
  
  dh_dz_ind = _Ind("denthalpy_D2");
  if (dh_dz_ind < 0)/* if there is no enthalpy defined in the patch */
    return 1;
  interp_dh_dz = init_interpolation();
  interp_dh_dz->field = patch->pool[dh_dz_ind];
  interp_dh_dz->X = X[0];
  interp_dh_dz->Y = X[1];
  interp_dh_dz->Z = X[2];
  interp_dh_dz->XYZ_dir_flag = 1;
  plan_interpolation(interp_dh_dz);
  dh_dz = execute_interpolation(interp_dh_dz);
  free_interpolation(interp_dh_dz);
  
  /* Grad h . r^ = dh/dr */
  return N[0]*dh_dx+N[1]*dh_dy+N[2]*dh_dz;
}

/* find NS surface using h = 1 */
static void find_NS_surface(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  GridParams->NS_R_type = GetParameterS_E("NS_surface_finder_method");
  
  /* find NS surface using cubed spherical points */
  //if (strstr_i(GridParams->NS_R_type,"CubedSpherical"))
    //find_NS_surface_CS_method_CS(grid,GridParams);
    
  /* find NS surface using spherical harmonic points */
  if (strstr_i(GridParams->NS_R_type,"SphericalHarmonic"))
    find_NS_surface_Ylm_method_CS(grid,GridParams);
  else
    abortEr(NO_OPTION);
}

/* extrapolate fluid fields outside NS */
static void extrapolate_fluid_fields_outsideNS(Grid_T *const grid)
{
  if (strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
  {
    /* extrapolate fluid fields outside of NS in case their value needed. */
    extrapolate_fluid_fields_outsideNS_CS(grid);
  }
  else
    abortEr(NO_OPTION);
  
}
