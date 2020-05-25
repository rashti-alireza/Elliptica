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
  
  if (!grid_prev)/* if grid is empty come up with an initialization */
  {
    /* if we wanna use checkpoint file */
    if (Pcmps("BH_NS_initialization","checkpoint_file"))
      grid_next = load_checkpoint_file();
    
    /* can we resume from a useful checkpoint file */
    else if (IsThereAnyUsefulCheckpointFile())
      grid_next = load_checkpoint_file();
      
    /* if we use TOV and Kerr-Schil black hole approximation */
    else if (Pcmps("BH_NS_initialization","TOV_KerrSchild"))
      grid_next = TOV_KerrSchild_approximation();
      
    else
      Error0(NO_OPTION);
  }
  else/* use previous grid to make the next one with new adjustments */
  {
    grid_next = make_next_grid_using_previous_grid(grid_prev);
  }
  
  return grid_next;   
}

/* loading the grid and parameters from checkpoint file.
// there is no check at this stage, we assume every thing is consistence */
static Grid_T *load_checkpoint_file(void)
{
  /* print some descriptions */
  pr_line_custom('=');
  printf("{ Initializing from checkpoint file ...\n");
  
  Grid_T *grid = 0;
  FILE *file   = 0;
  const char *const checkpoint_file_path = Pgets("checkpoint_file_path");
  
  if (access(checkpoint_file_path,F_OK))/* if file does not exist */
    Error1("Checkpoint file does not exist at\n%s\n",checkpoint_file_path);
    
  file = Fopen(checkpoint_file_path,"r");
  IsNull(file);
  
  grid = bbn_init_from_checkpoint(file);
  
  fclose(file);
  
  printf("} Initializing from checkpoint file ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
  return grid;
}

/* check if there is a consistence checkpoint file to be used 
// for initialization.
// -> return value : 1 if exists, 0 otherwise. */
static int IsThereAnyUsefulCheckpointFile(void)
{
  int ret = 0;
  const int FOLDER_TYPE   = 4;
  const unsigned MAX_LIST_NUM = 10;
  DIR *prev_dir;
  struct dirent *ent;
  struct stat st = {0};/* status of files */
  FILE *checkpoint_file = 0;
  Parameter_T *par;
  const char *const cur_out_dir = Pgets("output_directory_path");
  const char *const cur_folder_name  = strrchr(cur_out_dir,'/')+1;
  const char *const cur_folder_affix = strrchr(cur_folder_name,'_')+1;
  int cur_folder_index = atoi(cur_folder_affix);
  char prev_out_dir[MAX_ARRx4];
  char prev_folder_name[MAX_ARRx2];
  char folder_stem[MAX_ARR];
  char out_dir_stem[MAX_ARRx2];
  char *prev_data_folder_list[MAX_LIST_NUM];
  char prev_data_file_path[MAX_ARRx5];
  char *aux,str[MAX_ARRx5];
  long latest_mtime = 0;
  unsigned count,i;
  
  /* if there is no previous folder */
  if (!cur_folder_index)
    return 0;
  
  /* if there is a previous folder check the contents
  // see if there is any useful checkpoint file and finally 
  // find the latest one. */
  
  sprintf(out_dir_stem,"%s",cur_out_dir);
  aux    = strrchr(out_dir_stem,'/');
  aux[0] = '\0';
  
  sprintf(folder_stem,"%s",cur_folder_name);
  aux    = strrchr(folder_stem,'_');
  aux[0] = '\0';
  
  sprintf(prev_folder_name,"%s"FOLDER_AFFIX,
          folder_stem,cur_folder_index-1);
  
  sprintf(prev_out_dir,"%s/%s",out_dir_stem,prev_folder_name);
  
  /* open previous directory */        
  prev_dir = opendir(prev_out_dir);
  if (!prev_dir)/* if cannot be opened */
    return 0;
  
  count = 0;  
  ent = readdir(prev_dir);
  while (ent)
  {
    if (ent->d_type == FOLDER_TYPE)/* if this a folder */
    {
      if (!strcmp(ent->d_name,".") || !strcmp(ent->d_name,".."))
      {
        ent = readdir(prev_dir);
        continue;
      }
      assert(count < MAX_LIST_NUM);/* if it grater increase MAX_LIST_NUM */
      prev_data_folder_list[count] = dup_s(ent->d_name);
      count++;
    }
    ent = readdir(prev_dir);
  }
  closedir(prev_dir);  
  
  /* having found the data folder, find the latest checkpoint file */
  latest_mtime = 0;
  for (i = 0; i < count; ++i)
  {
    sprintf(str,"%s/%s/%s",prev_out_dir,
       prev_data_folder_list[i],CHECKPOINT_FILE_NAME);
    
    /* if the file exists */
    if(!stat(str, &st))
    {
      /* if it is later */
      if (st.st_mtime > latest_mtime)
      {
        /* if the checkpoint file written completely */
        if (bbn_IsCheckpointFileCompleted(str))
        {
          sprintf(prev_data_file_path,"%s",str);
          latest_mtime = st.st_mtime;
          ret = 1;/* ---> where we set ret = 1 */
        }
      }
    }/* end of if(!stat(str, &st)) */
  }
  
  /* free  */
  for (i  = 0; i < count; ++i)
    free(prev_data_folder_list[i]);
  
  /* some checks to make sure the checkpoint file is valid */
  /* NS & BH masses, NS & BH spins, EoS, separation, BH_KerrSchild_RollOff */
  if (ret && 0)
  {
    checkpoint_file = Fopen(prev_data_file_path,"r");
    IsNull(checkpoint_file);
    const char *check_quantities[] = {"NS_baryonic_mass",
                                      "NS_Omega_U0",
                                      "NS_Omega_U1",
                                      "NS_Omega_U2",
                                      "EoS_K",
                                      "BH_irreducible_mass",
                                      "BH_chi_U2",
                                      "BH_NS_separation",
                                      "BH_KerrSchild_RollOff",
                                      0};
    i = 0;
    while (check_quantities[i])
    {
      par = bbn_parameter_query_from_checkpoint_file(check_quantities[i],checkpoint_file);
      if (!Pcmps(check_quantities[i],par->rv)) 
      {
        printf("~> checkpoint file is not matched for %s\n",check_quantities[i]);
        fflush(stdout);
        ret = 0;/* --> set ret = 0 */
        free_given_parameter(par);
        break;
      }
      free_given_parameter(par);
      i++;
    }
    fclose(checkpoint_file);
  }
  
  /* set the path of checkpoint file and some preparation */
  if (ret)
  {
    Psets("checkpoint_file_path",prev_data_file_path);
    
    printf("~> checkpoint file found at:\n%s\n",prev_data_file_path);
    
    /* remove the current directory */
    sprintf(str,"rm -rf %s",cur_out_dir);
    shell_command(str);
  }
  
  return ret;
}

/* finding different quantities and then make the next grid using previous grid
// first find the values of the following parameters and some adjustment:
// the followings are not in order!
//
// . Euler equation constant.
// . orbital angular velocity
// . center of rotation (center of mass)
// . NS center
// . drag NS to the center
// . find NS surface
// . BH_radius
// . Omega_BH
//
// ->return value: the next grid called 'grid_next' */
static Grid_T *make_next_grid_using_previous_grid(Grid_T *const grid_prev)
{
  pr_line_custom('=');
  printf("{ Initializing the Fields and Grid Using Previous Solutions ...\n");
  
  Grid_T *grid_next = 0;
  struct Grid_Params_S *GridParams = init_GridParams();/* adjust some pars for construction of next grid */
  GridParams->grid_prev = grid_prev;
  
  /* NOTE: the order of function calls are important */

  /* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S */
  bbn_update_stress_energy_tensor(grid_prev,0);
  
  /* find Euler equation constant to meet NS baryonic mass */
  find_Euler_eq_const(grid_prev);
  
  /* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S */
  bbn_update_stress_energy_tensor(grid_prev,0);
  
  /* adjust the apparent horizon radius to acquire the desired BH mass */
  adjust_AH_radius(grid_prev,GridParams);
  
  /* P_ADM control */
  P_ADM_control(grid_prev);
  
  /* adjust the Omega_BH to acquire the desired BH spin
  // NOTE: this function should be after adjust_AH_radius 
  // and P_ADM_control since it needs some parameters such as
  // BH_irreducible_mass_current and BH_Px_ADM etc. */
  adjust_BH_Omega(grid_prev,GridParams);
 
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
  
  /* make new grid with new parameters */
  grid_next = creat_bbn_grid_CS(GridParams);
  
  /* fields: */
  /* adding all of the fields needed for construction of Initial Data */
  bbn_add_fields(grid_next);
  
  /* populating the free data part of initial data that we chose ourself */
  bbn_populate_free_data(grid_next);

  /* use previous grid to interpolate values of the fields for 
  // the next grid  and initialzing some other fields */
  interpolate_and_initialize_to_next_grid(grid_next,grid_prev);
  
  /* make normal vectorn on BH horizon 
  // note: this MUST be before "bbn_partial_derivatives_fields" */
  bbn_make_normal_vector_on_BH_horizon(grid_next);
  
  /* taking partial derivatives of the fields needed for equations */
  bbn_partial_derivatives_fields(grid_next);
  
  /* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S */
  bbn_update_stress_energy_tensor(grid_next,1);
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  bbn_update_Aij(grid_next);
  
  /* freeing */
  free_Grid_Params_S(GridParams);
  
  printf("} Initializing the Fields and Grid Using Previous Solutions ==> Done.\n");
  pr_line_custom('=');
  
  return grid_next;
}

/* keep NS center fixed at (0,C,0) by adjusting enthalpy */
static void keep_NS_center_fixed(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Adjusting NS center ...\n");
  
  const double C    = -0.5*Pgetd("BH_NS_separation");
  struct NC_Center_RootFinder_S par[1] = {0};
  Root_Finder_T root_finder[1] = {0};
  const double x_center[3] = {0,C,0};
  double dhx0,dhz0;
  double dh1[3] = {0},dh2[3] = {0};
  
  par->patch = GetPatch("left_central_box",grid);
  par->root_finder = root_finder;
  
  dh1[0] = dh_dx0_root_finder_eq(par,x_center);
  dh1[1] = dh_dx1_root_finder_eq(par,x_center);
  dh1[2] = dh_dx2_root_finder_eq(par,x_center);

  dhx0 = dh1[0];
  dhz0 = dh1[2];;
  
  if (Pcmps("NS_adjust_center_method","draw_enthalpy"))
  {
    adjust_NS_center_draw_enthalpy(grid);
  }
  else if (Pcmps("NS_adjust_center_method","tune_enthalpy"))
  {
    adjust_NS_center_tune_enthalpy(grid,dhx0,dhz0);
  }
  else
    Error0(NO_OPTION);
  
  dh2[0] = dh_dx0_root_finder_eq(par,x_center);
  dh2[1] = dh_dx1_root_finder_eq(par,x_center);
  dh2[2] = dh_dx2_root_finder_eq(par,x_center);
  
  /* print initial values after adjustments */
  printf("Enthalpy derivatives at NS center after NS center adjustment:\n");
  
  printf("|--> dh(%g,%g,%g)/dx = %+g\n",
    x_center[0],x_center[1],x_center[2],dh2[0]);
  printf("|--> dh(%g,%g,%g)/dy = %+g\n",
    x_center[0],x_center[1],x_center[2],dh2[1]);
  printf("|--> dh(%g,%g,%g)/dz = %+g\n",
    x_center[0],x_center[1],x_center[2],dh2[2]);
    
  printf("Changes in enthalpy derivatives after NS center adjustment:\n");
  printf("|--> dh2/dx-dh1/dx = %+g\n",dh2[0]-dh1[0]);
  printf("|--> dh2/dy-dh1/dy = %+g\n",dh2[1]-dh1[1]);
  printf("|--> dh2/dz-dh1/dz = %+g\n",dh2[2]-dh1[2]);
  
  printf("} Adjusting NS center  ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* controlling P_ADM */
static void P_ADM_control(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Adjusting ADM momentums ...\n");
  
  char *adjust[3];
  const char *const par = Pgets("P_ADM_control_method");
  
  printf("|--> %s\n",par);
  
  parse_adjust_parameter(par,adjust);
  
  void (*P_ADM_control_0)(Grid_T *const grid) =
                              get_func_P_ADM_adjustment(adjust[0]);
  
  void (*P_ADM_control_1)(Grid_T *const grid) =
                              get_func_P_ADM_adjustment(adjust[1]);
                              
  void (*P_ADM_control_2)(Grid_T *const grid) =
                              get_func_P_ADM_adjustment(adjust[2]);
  
  /* update P_ADM and J_ADM momentum parameters */
  Observable_T *obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  double p1[3] = {0};
  double p2[3] = {0};
  double j_adm[3] = {0};
  
  obs->quantity = "ADM(P,J)|BBN";
  plan_observable(obs);
  
  /* get previous P_ADMs */
  p1[0] = Pgetd("Px_ADM");
  p1[1] = Pgetd("Py_ADM");
  p1[2] = Pgetd("Pz_ADM");
  
  /* get the current P_ADMs */
  p2[0] = obs->Px(obs);
  p2[1] = obs->Py(obs);
  p2[2] = obs->Pz(obs);
  
  /* get the current J_ADMs  */
  j_adm[0] = obs->Jx(obs);
  j_adm[1] = obs->Jy(obs);
  j_adm[2] = obs->Jz(obs);
  
  printf("|--> Current P_ADM = (%e,%e,%e)\n",p2[0],p2[1],p2[2]);
  printf("|--> Current J_ADM = (%e,%e,%e)\n",j_adm[0],j_adm[1],j_adm[2]);
  
  Psetd("Px_ADM",p2[0]);
  Psetd("Py_ADM",p2[1]);
  Psetd("Pz_ADM",p2[2]);
  
  Psetd("Px_ADM_prev",p1[0]);
  Psetd("Py_ADM_prev",p1[1]);
  Psetd("Pz_ADM_prev",p1[2]);
  
  Psetd("Jx_ADM",j_adm[0]);
  Psetd("Jy_ADM",j_adm[1]);
  Psetd("Jz_ADM",j_adm[2]);
  
  free_observable(obs);
  
  /* NS adms */
  obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "ADM(P,J)|NS";
  plan_observable(obs);
  Psetd("NS_Px_ADM",obs->Px(obs));
  Psetd("NS_Py_ADM",obs->Py(obs));
  Psetd("NS_Pz_ADM",obs->Pz(obs));
  Psetd("NS_Jx_ADM",obs->Jx(obs));
  Psetd("NS_Jy_ADM",obs->Jy(obs));
  Psetd("NS_Jz_ADM",obs->Jz(obs));
  free_observable(obs);
  
  /* BH adms */
  obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "ADM(P,J)|BH";
  plan_observable(obs);
  Psetd("BH_Px_ADM",obs->Px(obs));
  Psetd("BH_Py_ADM",obs->Py(obs));
  Psetd("BH_Pz_ADM",obs->Pz(obs));
  Psetd("BH_Jx_ADM",obs->Jx(obs));
  Psetd("BH_Jy_ADM",obs->Jy(obs));
  Psetd("BH_Jz_ADM",obs->Jz(obs));
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
  
  printf("} Adjusting ADM momentums ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* getting adjustment str, returns the relevant function. */
fAdjustment_t *get_func_P_ADM_adjustment(const char *const adjust)
{
  fAdjustment_t *f = 0;
  
  if (!adjust)
  {
    f = 0;
  }
  else if (strcmp_i(adjust,"BH_Vz"))
  {
    f = Pz_ADM_is0_by_BH_Vz;
  }
  else if (strcmp_i(adjust,"x_CM&y_CM") || strcmp_i(adjust,"y_CM&x_CM"))
  {
    f = Pxy_ADM_is0_by_xy_CM_roots;
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
  else if (strcmp_i(adjust,"Px_BH_center"))
  {
    f = Px_ADM_is0_by_BH_center_y;
  }
  else if (strcmp_i(adjust,"Py_BH_center"))
  {
    f = Py_ADM_is0_by_BH_center_x;
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
    Error0(NO_OPTION);
  
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

/* adjust various quantities according to parameter:
// "force_balance_equation" */
static void force_balance_eq(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Applying force balance equation ...\n");

  struct NC_Center_RootFinder_S dh_par[1] = {0};
  Root_Finder_T root_finder[1]            = {0};
  const double D            = Pgetd("BH_NS_separation");
  const double NS_center[3] = {0,-D/2,0};/* since we keep the NS center always here */
  char *adjust[3];
  const char *const par = Pgets("force_balance_equation");
  double dh1[3] = {0},dh2[3] = {0};
  dh_par->patch = GetPatch("left_central_box",grid);
  dh_par->root_finder = root_finder;
  
  /* initial values before adjustments */
  dh1[0] = dh_dx0_root_finder_eq(dh_par,NS_center);
  dh1[1] = dh_dx1_root_finder_eq(dh_par,NS_center);
  dh1[2] = dh_dx2_root_finder_eq(dh_par,NS_center);
  
  printf("%s\n",par);
  
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
  
  /* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S */
  bbn_update_stress_energy_tensor(grid,0);
  
  dh2[0] = dh_dx0_root_finder_eq(dh_par,NS_center);
  dh2[1] = dh_dx1_root_finder_eq(dh_par,NS_center);
  dh2[2] = dh_dx2_root_finder_eq(dh_par,NS_center);
  
  /* print initial values after adjustments */
  printf("Enthalpy derivatives at NS center after force balance eq.:\n");
  
  printf("|--> dh(%g,%g,%g)/dx = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh2[0]);
  printf("|--> dh(%g,%g,%g)/dy = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh2[1]);
  printf("|--> dh(%g,%g,%g)/dz = %+g\n",
    NS_center[0],NS_center[1],NS_center[2],dh2[2]);
    
  printf("Changes in enthalpy derivatives after force balance eq.:\n");
  printf("|--> dh2/dx-dh1/dx = %+g\n",dh2[0]-dh1[0]);
  printf("|--> dh2/dy-dh1/dy = %+g\n",dh2[1]-dh1[1]);
  printf("|--> dh2/dz-dh1/dz = %+g\n",dh2[2]-dh1[2]);
  
  printf("} Applying force balance equation ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* adjust Px ADM by changing the center of BH. */
static void Px_ADM_is0_by_BH_center_y(Grid_T *const grid)
{
  const double W1 = Pgetd("Solving_Field_Update_Weight");
  const double W2 = 1-W1;
  const double dP = Pgetd("P_ADM_control_tolerance");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double BH_center_y = Pgetd("BH_center_y");
  double px0,px,BH_center_y_new,dy;
  
  /* get P_ADM */
  px = Pgetd("Px_ADM");
  px0  = Pgetd("Px_ADM_prev");
  
  /* changing center of mass */
  dy              = px/Omega_BHNS;
  BH_center_y_new = W2*BH_center_y+dy*W1;
  
  const double dPx_Px = fabs(px0-px)/fabs(px);
  /* having found new x_CM now update */
  if (GRT(dPx_Px,dP))
  {
    printf("|--> |Px_ADM2 - Px_ADM1|/|Px_ADM2| = %g > %g\n",dPx_Px,dP);
    Psetd("BH_center_y",BH_center_y_new);
  }
  else
  {
    printf("|--> |Px_ADM2 - Px_ADM1|/|Px_ADM2| = %g <= %g\n"
           "     |--> no BH center y-axis update.\n",dPx_Px,dP);
  }
  
  UNUSED(grid);
}

/* adjust Py ADM by changing the center of BH. */
static void Py_ADM_is0_by_BH_center_x(Grid_T *const grid)
{
  const double W1 = Pgetd("Solving_Field_Update_Weight");
  const double W2 = 1-W1;
  const double dP = Pgetd("P_ADM_control_tolerance");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double BH_center_x = Pgetd("BH_center_x");
  double py0,py,BH_center_x_new,dx;
  
  /* get P_ADM */
  py = Pgetd("Py_ADM");
  py0  = Pgetd("Py_ADM_prev");
  
  /* changing center of mass */
  dx              = -py/Omega_BHNS;
  BH_center_x_new = W2*BH_center_x+dx*W1;
  
  const double dPy_Py = fabs(py0-py)/fabs(py);
  /* having found new x_CM now update */
  if (GRT(dPy_Py,dP))
  {
    printf("|--> |Py_ADM2 - Py_ADM1|/|Py_ADM2| = %g > %g\n",dPy_Py,dP);
    Psetd("BH_center_x",BH_center_x_new);
  }
  else
  {
    printf("|--> |Py_ADM2 - Py_ADM1|/|Py_ADM2| = %g <= %g\n"
           "     |--> no BH center x-axis update.\n",dPy_Py,dP);
  }
  
  UNUSED(grid);
}

/* adjust the boost velocity at the outer boundary to diminish P_ADM
// it only makes changes in the specified direction x. */
static void Px_ADM_is0_by_x_boost(Grid_T *const grid)
{
  const double SMALL_FAC = 1E-2;
  const double dP   = Pgetd("P_ADM_control_tolerance");
  const double W    = Pgetd("Solving_Field_Update_Weight");
  double p1[3] = {0};
  double p2[3] = {0};
  double v1[3] = {0};
  double v2[3] = {0};
  double v0[3] = {0};
  double  v[3] = {0};
  static unsigned iter = 0;
  
  /* get previous P_ADMs */
  p1[0] = Pgetd("Px_ADM_prev");
  
  /* get the current P_ADMs */
  p2[0] = Pgetd("Px_ADM");
  
  if (iter == 0)
  {
    Psetd("v1_boost_x",p2[0]*SMALL_FAC);
  }
  else if (iter == 1)
  {
    Psetd("v2_boost_x",p2[0]*SMALL_FAC);
  }
  else if (iter > 1)
  {
    /* get the penultimate boost velocity */
    v1[0] = Pgetd("v1_boost_x");
    
    /* get the ultimate boost velocity */
    v2[0] = Pgetd("v2_boost_x");
    
    /* get the boost velocity */
    v0[0] = Pgetd("v*_boost_x");
    
    /* calculate the new boost velocity */
    v[0] = (v2[0]*p1[0]-v1[0]*p2[0])/(p1[0]-p2[0]);
    
    /* take cure of 0 denominator */
    if (EQL(p1[0],p2[0])) v[0] = 0;
    
    /* change the boost velocity relaxed */
    v[0] = W*v[0]+(1-W)*v0[0];
    
    /* update parameters */
    Psetd("v1_boost_x",v2[0]);
    Psetd("v2_boost_x",v[0]);
    
    const double dPx_Px = fabs(p2[0]-p1[0])/fabs(p2[0]);
    /* if change in momentum is big */
    if (GRT(dPx_Px,dP))
    {
      printf("|--> |Px_ADM2 - Px_ADM1|/|Px_ADM2| = %g > %g\n",dPx_Px,dP);
      Psetd("v*_boost_x",v[0]);
    }
    else
    {
      printf("|--> |Px_ADM2 - Px_ADM1|/|Px_ADM2| = %g <= %g\n"
             "     |--> no v*_boost_x update.\n",dPx_Px,dP);
    }
      
  }
  
  iter++;
  UNUSED(grid);
} 

/* adjust the boost velocity at the outer boundary to diminish P_ADM
// it only makes changes in the specified direction y. */
static void Py_ADM_is0_by_y_boost(Grid_T *const grid)
{
  const double SMALL_FAC = 1E-2;
  const double dP   = Pgetd("P_ADM_control_tolerance");
  const double W    = Pgetd("Solving_Field_Update_Weight");
  double p1[3] = {0};
  double p2[3] = {0};
  double v1[3] = {0};
  double v2[3] = {0};
  double v0[3] = {0};
  double  v[3] = {0};
  static unsigned iter = 0;
  
  /* get previous P_ADMs */
  p1[1] = Pgetd("Py_ADM_prev");
  
  /* get the current P_ADMs */
  p2[1] = Pgetd("Py_ADM");
  
  if (iter == 0)
  {
    Psetd("v1_boost_y",p2[1]*SMALL_FAC);
  }
  else if (iter == 1)
  {
    Psetd("v2_boost_y",p2[1]*SMALL_FAC);
  }
  else if (iter > 1)
  {
    /* get the penultimate boost velocity */
    v1[1] = Pgetd("v1_boost_y");
    
    /* get the ultimate boost velocity */
    v2[1] = Pgetd("v2_boost_y");
    
    /* get the boost velocity */
    v0[1] = Pgetd("v*_boost_y");
    
    /* calculate the new boost velocity */
    v[1] = (v2[1]*p1[1]-v1[1]*p2[1])/(p1[1]-p2[1]);
    
    /* take cure of 0 denominator */
    if (EQL(p1[1],p2[1])) v[1] = 0;
    
    /* change the boost velocity relaxed */
    v[1] = W*v[1]+(1-W)*v0[1];
    
    /* update parameters */
    Psetd("v1_boost_y",v2[1]);
    Psetd("v2_boost_y",v[1]);
    
    const double dPy_Py = fabs(p2[1]-p1[1])/fabs(p2[1]);
    if (GRT(dPy_Py,dP))
    {
      printf("|--> |Py_ADM2 - Py_ADM1|/|Py_ADM2| = %g > %g\n",dPy_Py,dP);
      Psetd("v*_boost_y",v[1]);
    }
    else
    {
      printf("|--> |Py_ADM2 - Py_ADM1|/|Py_ADM2| = %g <= %g\n"
           "     |--> no v*_boost_y update.\n",dPy_Py,dP);
    }
      
  }
  
  iter++;
  UNUSED(grid);
} 

/* adjust the boost velocity at the outer boundary to diminish P_ADM
// it only makes changes in the specified direction z. */
static void Pz_ADM_is0_by_z_boost(Grid_T *const grid)
{
  const double SMALL_FAC = 1E-2;
  const double dP   = Pgetd("P_ADM_control_tolerance");
  const double W    = Pgetd("Solving_Field_Update_Weight");
  double p1[3] = {0};
  double p2[3] = {0};
  double v1[3] = {0};
  double v2[3] = {0};
  double v0[3] = {0};
  double  v[3] = {0};
  static unsigned iter = 0;
  
  /* get previous P_ADMs */
  p1[2] = Pgetd("Pz_ADM_prev");
  
  /* get the current P_ADMs */
  p2[2] = Pgetd("Pz_ADM");
  
  if (iter == 0)
  {
    Psetd("v1_boost_z",p2[2]*SMALL_FAC);
  }
  else if (iter == 1)
  {
    Psetd("v2_boost_z",p2[2]*SMALL_FAC);
  }
  else if (iter > 1)
  {
    /* get the penultimate boost velocity */
    v1[2] = Pgetd("v1_boost_z");
    
    /* get the ultimate boost velocity */
    v2[2] = Pgetd("v2_boost_z");
    
    /* get the boost velocity */
    v0[2] = Pgetd("v*_boost_z");
    
    /* calculate the new boost velocity */
    v[2] = (v2[2]*p1[2]-v1[2]*p2[2])/(p1[2]-p2[2]);
    
    /* take cure of 0 denominator */
    if (EQL(p1[2],p2[2])) v[2] = 0;
    
    /* change the boost velocity relaxed */
    v[2] = W*v[2]+(1-W)*v0[2];
    
    /* update parameters */
    Psetd("v1_boost_z",v2[2]);
    Psetd("v2_boost_z",v[2]);
    
    const double dPz_Pz = fabs(p2[2]-p1[2])/fabs(p2[2]);
    if (GRT(dPz_Pz,dP))
    {
      printf("|--> |Pz_ADM2 - Pz_ADM1|/|Pz_ADM2| = %g > %g\n",dPz_Pz,dP);
      Psetd("v*_boost_z",v[2]);
    }
    else
    {
      printf("|--> |Pz_ADM2 - Pz_ADM1|/|Pz_ADM2| = %g <= %g\n"
           "     |--> no v*_boost_z update.\n",dPz_Pz,dP);
    }
  }
  
  iter++;
  UNUSED(grid);
} 

/* adjust the center of NS at the designated point, in case it moved. 
// in this method, we tune enthalpy values such that the derivative of 
// the enthalpy be 0 at (0,NS_C,0), using Taylor expansion, assuming dh/dy|center = 0. */
static void adjust_NS_center_tune_enthalpy(Grid_T *const grid,const double dhx0,const double dhz0)
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
      WRITE_v(enthalpy)
      //READ_v(denthalpy_D2)
      //READ_v(denthalpy_D0)
    
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double x = patch->node[ijk]->x[0];
        double z = patch->node[ijk]->x[2];
        
        enthalpy[ijk] = enthalpy[ijk]-dhx0*x-dhz0*z;
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
  const double C    = -0.5*Pgetd("BH_NS_separation");
  sprintf(par_name,"grid%u_NS_center_xyz",grid->gn);
  NS_center = Pgetdd(par_name);
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
    REALLOC_v_WRITE_v(shifted_enthalpy);
    
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
      if (!IsItNSPatch(patchp) && !IsItNSSurroundingPatch(patchp))
        shifted_enthalpy[ijk] = 1;
      else
        shifted_enthalpy[ijk] = interpolate_from_patch_prim("enthalpy",Xp,patchp);
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
  force_balance_eq_root_finders(grid,0,"BH_NS_angular_velocity");
}

/* force_balance_equation : adjust Omega at direction d/dy */
static void force_balance_ddy_Omega(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,1,"BH_NS_angular_velocity");
}

/* force_balance_equation : adjust Omega at direction d/dz */
static void force_balance_ddz_Omega(Grid_T *const grid)
{
  force_balance_eq_root_finders(grid,2,"BH_NS_angular_velocity");
}

/* find parameter par using force balance equation in direction dir */
static void force_balance_eq_root_finders(Grid_T *const grid,const int dir, const char *const par)
{
  const double D            = Pgetd("BH_NS_separation");
  const double Vr           = Pgetd("BH_NS_infall_velocity");
  const double NS_center[3] = {0,-D/2,0};/* since we keep the NS center always here */
  const double RESIDUAL     = sqrt(Pgetd("RootFinder_Tolerance"));
  const double Omega_BHNS   = Pgetd("BH_NS_angular_velocity");
  const double y_CM         = Pgetd("y_CM");
  const double x_CM         = Pgetd("x_CM");
  const double W1  = Pgetd("Solving_Field_Update_Weight");
  const double W2  = 1-W1;
  double *new_par,old_par = 0;
  double guess[1],X[3];
  struct Force_Balance_RootFinder_S params[1];
  Patch_T *patch;
  char desc[1000] = {'\0'};
  
  patch = GetPatch("left_central_box",grid);
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
  if (strcmp_i("BH_NS_angular_velocity",par))
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
    Error0(NO_OPTION);
  
  sprintf(desc,"Solving Force Balance Eq. for '%s' at direction 'x^%d'",par,dir);
  
  Root_Finder_T *root = init_root_finder(1);
  root->description   = desc;
  root->verbose       = 1;
  root->type          = Pgets("RootFinder_Method");
  root->tolerance     = Pgetd("RootFinder_Tolerance");
  root->MaxIter       = (unsigned)Pgeti("RootFinder_Max_Number_of_Iteration");
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
  Psetd(par,new_par[0]);
  
  /* since B1 has been changed let's update the pertinent fields */
  update_B1_dB1_Beta_dBete_Aij_dAij(grid);
  
  free_root_finder(root);
  free(new_par);
}

/* find the NS center using d(enthalpy)/dx^i = 0 */
static void find_NS_center(Grid_T *const grid)
{
  printf("{ Finding NS center ...\n");
  
  double *NS_center;
  Root_Finder_T *root = init_root_finder(3);
  struct NC_Center_RootFinder_S params[1];
  const double RESIDUAL = sqrt(Pgetd("RootFinder_Tolerance"));
  Flag_T success_f = NONE;
  double guess[3];/* initial guess for root finder */
  char par_name[1000];
  unsigned p;
  
  guess[0] = guess[2] = 0;
  guess[1] = -0.5*Pgetd("BH_NS_separation");
  params->root_finder = root;
  root->type        = Pgets("RootFinder_Method");
  root->tolerance   = Pgetd("RootFinder_Tolerance");
  root->MaxIter     = (unsigned)Pgeti("RootFinder_Max_Number_of_Iteration");
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
      sprintf(par_name,"grid%u_NS_center_xyz",grid->gn);
      add_parameter_array(par_name,NS_center,3);
      /* save the patch stem where the NS center takes place */
      sprintf(par_name,"grid%u_NS_center_patch",grid->gn);
      char *stem = strstr(patch->name,"_");
      assert(stem);
      stem++;
      add_parameter(par_name,stem);
      
      printf("|--> Current NS center found at (%g,%g,%g)\n",NS_center[0],NS_center[1],NS_center[2]);
      printf("|--> Change in x direction = %+g\n",NS_center[0]-guess[0]);
      printf("|--> Change in y direction = %+g\n",NS_center[1]-guess[1]);
      printf("|--> Change in z direction = %+g\n",NS_center[2]-guess[2]);
      
      free(NS_center);
      break;
    }
    free(NS_center);
  }
  if (success_f != YES)
  {
    print_root_finder_exit_status(root);
    Error0("NS center could not be found.\n");
  }
  free_root_finder(root);
  
  printf("} Finding NS center ==> Done.\n");  
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
  pr_line_custom('=');
  printf("{ Finding Euler equation constant using NS baryonic mass ...\n");
  
  Root_Finder_T *root = init_root_finder(1);
  const double W1  = Pgetd("Solving_Field_Update_Weight");
  const double W2  = 1-W1;
  double *Euler_const = 0;
  double guess[1] = {Pgetd("Euler_equation_constant")};
  const double RESIDUAL = sqrt(Pgetd("RootFinder_Tolerance"));
  struct Euler_eq_const_RootFinder_S params[1] = {0};
  Observable_T *obs = 0;
  double bar_mass,adm_mass,kommar_mass;
  
  bar_mass = bbn_NS_baryonic_mass(grid,guess[0]);
  obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "ADM(M)|NS";
  plan_observable(obs);
  adm_mass = obs->M(obs);
  free_observable(obs);
  
  obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "Kommar(M)|NS";
  plan_observable(obs);
  kommar_mass = obs->M(obs);
  free_observable(obs);

  printf("|--> current NS baryonic mass = %e\n",bar_mass);
  printf("|--> current NS ADM mass      = %e\n",adm_mass);
  printf("|--> current NS Kommar mass   = %e\n",kommar_mass);
  
  Psetd("NS_baryonic_mass_current",bar_mass);
  Psetd("NS_ADM_mass",adm_mass);
  Psetd("NS_Kommar_mass",kommar_mass);
  
  params->grid = grid;
  params->NS_baryonic_mass = Pgetd("NS_baryonic_mass");
  
  root->type        = Pgets("RootFinder_Method");
  root->tolerance   = Pgetd("RootFinder_Tolerance");
  root->MaxIter     = (unsigned)Pgeti("RootFinder_Max_Number_of_Iteration");
  root->x_gss       = guess;
  root->params      = params;
  root->f[0]        = Euler_eq_const_rootfinder_eq;
  root->verbose     = 1;
  plan_root_finder(root);
  
  Euler_const       = execute_root_finder(root);
  /* if root finder is not OK for some reason */
  if (GRT(root->residual,RESIDUAL))
    Euler_const[0] = guess[0];/* don't update */
  else
    Euler_const[0] = W1*Euler_const[0]+W2*guess[0];
  
  Psetd("Euler_equation_constant",Euler_const[0]);
  free(Euler_const);
  free_root_finder(root);
  
  printf("} Finding Euler equation constant using NS baryonic mass ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}
/* find y_CM and x_CM by demanding Px_ADM = 0 and Py_ADM = 0
// using root finder. */
static void Pxy_ADM_is0_by_xy_CM_roots(Grid_T *const grid)
{
  printf("|--> Solving {Px(y_CM) = 0 && Py(x_CM) = 0} ...\n");
  const double W      = Pgetd("P_ADM_control_update_weight");
  const double dP     = Pgetd("P_ADM_control_tolerance");
  const unsigned Iter = (unsigned)Pgeti("P_ADM_control_iteration");
  Root_Finder_T *root = 0;
  struct PxPy_RootFinder_S params[1] = {0};
  const double x0[2] = {Pgetd("x_CM"),Pgetd("y_CM")};/* NOTE: index 0 is for x and 1 for y */
  const double guess[2] = {0,0};
  double x_new[2] = {0},*dx = 0;
  Grid_T *freedata_grid = 0;/* don't update for inside BH patches */
  Patch_T **freedata_patch = 0;/* all but inside BH patches */
  double p_adm[3] = {0},j_adm[3] = {0};
  unsigned i,p;
  
  /* populate Aij grid */
  freedata_grid = calloc(1,sizeof(*freedata_grid));
  IsNull(freedata_grid);
  
  i = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    if (!IsItInsideBHPatch(grid->patch[p]))
    {
      freedata_patch = realloc(freedata_patch,
                  (i+1)*sizeof(*freedata_patch));
      IsNull(freedata_patch);
      freedata_patch[i++] = grid->patch[p];
    }
  }
  freedata_grid->kind  = grid->kind;
  freedata_grid->patch = freedata_patch;
  freedata_grid->gn    = grid->gn;
  freedata_grid->np    = i;
  freedata_grid->nn    = UINT_MAX;
  
  params->grid  = grid;
  params->freedata_grid = freedata_grid;
  params->x_CM0 = x0[0];
  params->y_CM0 = x0[1];
  
  root = init_root_finder(2);
  root->verbose       = 1;
  root->type          = Pgets("RootFinder_Method");
  root->tolerance     = dP;
  root->MaxIter       = Iter;
  root->x_gss         = guess;
  root->params        = params;
  root->f[0]          = x_CM_root_of_Py;
  root->f[1]          = y_CM_root_of_Px;
  plan_root_finder(root);
  dx = execute_root_finder(root);
  free_root_finder(root);
  
  /* updating */
  x_new[0] = -W*dx[0]+x0[0];
  x_new[1] = -W*dx[1]+x0[1];
  Psetd("x_CM",x_new[0]);
  Psetd("y_CM",x_new[1]);
  
  /* updating the free data */
  bbn_populate_free_data(freedata_grid);
  update_B1_dB1_Beta_dBete_Aij_dAij(grid);

  free(dx);
  free(freedata_grid);
  free(freedata_patch);
  
  Observable_T *obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "ADM(P,J)|BBN";
  plan_observable(obs);
  
  /* get the current P_ADMs */
  p_adm[0] = obs->Px(obs);
  p_adm[1] = obs->Py(obs);
  p_adm[2] = obs->Pz(obs);
  
  /* get the current J_ADMs  */
  j_adm[0] = obs->Jx(obs);
  j_adm[1] = obs->Jy(obs);
  j_adm[2] = obs->Jz(obs);
  
  printf("|--> After CM update P_ADM = (%e,%e,%e)\n",p_adm[0],p_adm[1],p_adm[2]);
  printf("|--> After CM update J_ADM = (%e,%e,%e)\n",j_adm[0],j_adm[1],j_adm[2]);
  free_observable(obs);
}

/* solving Py = 0 by finding x_CM */
static double x_CM_root_of_Py(void *params,const double *const x)
{
  struct PxPy_RootFinder_S *const par = params;
  Grid_T *const grid    = par->grid;
  Grid_T *const freedata_grid = par->freedata_grid;
  Observable_T *obs = 0;
  const double x_cm = par->x_CM0+x[0];/* index 0 is for x_cm */
  double residual;

  /* updating free data and B's and related */
  Psetd_0print("x_CM",x_cm);
  bbn_free_data_gammas(freedata_grid);
  bbn_free_data_Gamma(freedata_grid);
  //bbn_free_data_dGamma(freedata_grid);
  //bbn_free_data_Ricci(freedata_grid);
  bbn_free_data_tr_KSKij(freedata_grid);
  update_B1_dB1_Beta_dBete_Aij_dAij(grid);

  obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "ADM(P,J)|BBN";
  plan_observable(obs);
  residual = obs->Py(obs);
  free_observable(obs);

  return residual;
}

/* solving Px = 0 by finding y_CM */
static double y_CM_root_of_Px(void *params,const double *const x)
{
  struct PxPy_RootFinder_S *const par = params;
  Grid_T *const grid    = par->grid;
  Grid_T *const freedata_grid = par->freedata_grid;
  Observable_T *obs = 0;
  const double y_cm = par->y_CM0+x[1];/* index 1 is for y_cm */
  double residual;

  /* updating free data and B's and related */
  Psetd_0print("y_CM",y_cm);
  bbn_free_data_gammas(freedata_grid);
  bbn_free_data_Gamma(freedata_grid);
  //bbn_free_data_dGamma(freedata_grid);
  //bbn_free_data_Ricci(freedata_grid);
  bbn_free_data_tr_KSKij(freedata_grid);
  update_B1_dB1_Beta_dBete_Aij_dAij(grid);

  obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "ADM(P,J)|BBN";
  plan_observable(obs);
  residual = obs->Px(obs);
  free_observable(obs);

  return residual;
}

/* find y_CM by demanding Px_ADM = 0 */
static void Px_ADM_is0_by_y_CM(Grid_T *const grid)
{
  double dy_CM = 0,px0,y_CM_new,px;
  const double W    = Pgetd("P_ADM_control_update_weight");
  const double dP   = Pgetd("P_ADM_control_tolerance");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double y_CM0 = Pgetd("y_CM");
  const double M_NS  = Pgetd("NS_baryonic_mass");
  const double M_BH  = Pgetd("BH_irreducible_mass");
  Grid_T *freedata_grid = 0;/* don't update for inside BH patches */
  Patch_T **freedata_patch = 0;/* all but inside BH patches */
  Observable_T *obs = 0;
  double p_adm[3] = {0},j_adm[3] = {0};
  unsigned i,p;
  
  printf("|--> adjusting Px_ADM by y_CM ...\n");
  
  /* populate Aij grid */
  freedata_grid = calloc(1,sizeof(*freedata_grid));
  IsNull(freedata_grid);
  
  i = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    if (!IsItInsideBHPatch(grid->patch[p]))
    {
      freedata_patch = realloc(freedata_patch,
                  (i+1)*sizeof(*freedata_patch));
      IsNull(freedata_patch);
      freedata_patch[i++] = grid->patch[p];
    }
  }
  freedata_grid->kind  = grid->kind;
  freedata_grid->patch = freedata_patch;
  freedata_grid->gn    = grid->gn;
  freedata_grid->np    = i;
  freedata_grid->nn    = UINT_MAX;
  
  /* get P_ADM */
  px  = Pgetd("Px_ADM");
  px0 = Pgetd("Px_ADM_prev");
  
  /* changing center of mass */
  dy_CM    = -px/(Omega_BHNS*(M_NS+M_BH));
  y_CM_new = y_CM0+dy_CM*W;
  
  const double dPx_Px = fabs(px0-px)/fabs(px);
  
  /* having found new x_CM now update */
  if (GRT(dPx_Px,dP))
  {
    printf("|--> |Px_ADM2 - Px_ADM1|/|Px_ADM2| = %g > %g\n",dPx_Px,dP);
    Psetd("y_CM",y_CM_new);
    bbn_populate_free_data(freedata_grid);
    update_B1_dB1_Beta_dBete_Aij_dAij(grid);
    obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
    obs->quantity = "ADM(P,J)|BBN";
    plan_observable(obs);
    
    /* get the current P_ADMs */
    p_adm[0] = obs->Px(obs);
    p_adm[1] = obs->Py(obs);
    p_adm[2] = obs->Pz(obs);
    
    /* get the current J_ADMs  */
    j_adm[0] = obs->Jx(obs);
    j_adm[1] = obs->Jy(obs);
    j_adm[2] = obs->Jz(obs);
    
    printf("|--> After CM update P_ADM = (%e,%e,%e)\n",p_adm[0],p_adm[1],p_adm[2]);
    printf("|--> After CM update J_ADM = (%e,%e,%e)\n",j_adm[0],j_adm[1],j_adm[2]);
    
    free_observable(obs);
  }
  else
  {
    printf("|--> |Px_ADM2 - Px_ADM1|/|Px_ADM2| = %g <= %g\n"
           "     |--> no y_CM update.\n",dPx_Px,dP);
  }
  
  free(freedata_grid);
  free(freedata_patch);
}

/* find BH_Vz by demanding Pz_ADM = 0.
// in effect, boosting BH in z direction to get Pz_ADM = 0. */
static void Pz_ADM_is0_by_BH_Vz(Grid_T *const grid)
{
  double  pz0,pz,BH_Vz_new,dBH_Vz;
  const double W     = Pgetd("P_ADM_control_update_weight");
  const double dP    = Pgetd("P_ADM_control_tolerance");
  const double BH_Vz = Pgetd("BH_Vz");
  const double M_BH  = Pgetd("BH_irreducible_mass");
  Grid_T *freedata_grid = 0;/* don't update for inside BH patches */
  Patch_T **freedata_patch = 0;/* all but inside BH patches */
  Observable_T *obs = 0;
  double p_adm[3] = {0},j_adm[3] = {0};
  unsigned i,p;
  
  printf("|--> adjusting Pz_ADM by BH_Vz ...\n");
  
  /* populate Aij grid */
  freedata_grid = calloc(1,sizeof(*freedata_grid));
  IsNull(freedata_grid);
  
  i = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    if (!IsItInsideBHPatch(grid->patch[p]))
    {
      freedata_patch = realloc(freedata_patch,
                  (i+1)*sizeof(*freedata_patch));
      IsNull(freedata_patch);
      freedata_patch[i++] = grid->patch[p];
    }
  }
  freedata_grid->kind  = grid->kind;
  freedata_grid->patch = freedata_patch;
  freedata_grid->gn    = grid->gn;
  freedata_grid->np    = i;
  freedata_grid->nn    = UINT_MAX;
  
  /* get P_ADM */
  pz  = Pgetd("Pz_ADM");
  pz0 = Pgetd("Pz_ADM_prev");
  /* adjust BH_Vz */
  //dBH_Vz    = -pz/M_BH;
  dBH_Vz    = pz/M_BH;
  BH_Vz_new = W*dBH_Vz+BH_Vz;
  
  const double dPz_Pz = fabs(pz0-pz)/fabs(pz);
  /* having found new BH_Vz now update */
  if (GRT(dPz_Pz,dP))
  {
    printf("|--> |Pz_ADM2 - Pz_ADM1|/|Pz_ADM2| = %g > %g\n",dPz_Pz,dP);
    Psetd("BH_Vz",BH_Vz_new);
    bbn_populate_free_data(freedata_grid);
    update_B1_dB1_Beta_dBete_Aij_dAij(grid);
    obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
    obs->quantity = "ADM(P,J)|BBN";
    plan_observable(obs);
    
    /* get the current P_ADMs */
    p_adm[0] = obs->Px(obs);
    p_adm[1] = obs->Py(obs);
    p_adm[2] = obs->Pz(obs);
    
    /* get the current J_ADMs  */
    j_adm[0] = obs->Jx(obs);
    j_adm[1] = obs->Jy(obs);
    j_adm[2] = obs->Jz(obs);
    
    printf("|--> After BH_Vz update P_ADM = (%e,%e,%e)\n",p_adm[0],p_adm[1],p_adm[2]);
    printf("|--> After BH_Vz update J_ADM = (%e,%e,%e)\n",j_adm[0],j_adm[1],j_adm[2]);
    
    free_observable(obs);
  }
  else
  {
    printf("|--> |Pz_ADM2 - Pz_ADM1|/|Pz_ADM2| = %g <= %g\n"
           "     |--> no BH_Vz update.\n",dPz_Pz,dP);
  }
  
  free(freedata_grid);
  free(freedata_patch);
}

/* find x_CM by demanding Py_ADM = 0 */
static void Py_ADM_is0_by_x_CM(Grid_T *const grid)
{
  double  dx_CM = 0,py0,x_CM_new,py;
  const double W    = Pgetd("P_ADM_control_update_weight");
  const double dP   = Pgetd("P_ADM_control_tolerance");
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double x_CM0 = Pgetd("x_CM");
  const double M_NS  = Pgetd("NS_baryonic_mass");
  const double M_BH  = Pgetd("BH_irreducible_mass");
  Grid_T *freedata_grid = 0;/* don't update for inside BH patches */
  Patch_T **freedata_patch = 0;/* all but inside BH patches */
  Observable_T *obs = 0;
  double p_adm[3] = {0},j_adm[3] = {0};
  unsigned i,p;
  
  printf("|--> adjusting Py_ADM by x_CM ...\n");
  
  /* populate Aij grid */
  freedata_grid = calloc(1,sizeof(*freedata_grid));
  IsNull(freedata_grid);
  
  i = 0;
  FOR_ALL_PATCHES(p,grid)
  {
    if (!IsItInsideBHPatch(grid->patch[p]))
    {
      freedata_patch = realloc(freedata_patch,
                  (i+1)*sizeof(*freedata_patch));
      IsNull(freedata_patch);
      freedata_patch[i++] = grid->patch[p];
    }
  }
  freedata_grid->kind  = grid->kind;
  freedata_grid->patch = freedata_patch;
  freedata_grid->gn    = grid->gn;
  freedata_grid->np    = i;
  freedata_grid->nn    = UINT_MAX;
  
  /* get P_ADM */
  py  = Pgetd("Py_ADM");
  py0 = Pgetd("Py_ADM_prev");
  
  /* changing center of mass */
  dx_CM    = py/(Omega_BHNS*(M_NS+M_BH));
  x_CM_new = x_CM0+dx_CM*W;
  
  const double dPy_Py = fabs(py0-py)/fabs(py);
  /* having found new x_CM now update */
  if (GRT(dPy_Py,dP))
  {
    printf("|--> |Py_ADM2 - Py_ADM1|/|Py_ADM2| = %g > %g\n",dPy_Py,dP);
    Psetd("x_CM",x_CM_new);
    bbn_populate_free_data(freedata_grid);
    update_B1_dB1_Beta_dBete_Aij_dAij(grid);
    obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
    obs->quantity = "ADM(P,J)|BBN";
    plan_observable(obs);
    
    /* get the current P_ADMs */
    p_adm[0] = obs->Px(obs);
    p_adm[1] = obs->Py(obs);
    p_adm[2] = obs->Pz(obs);
    
    /* get the current J_ADMs  */
    j_adm[0] = obs->Jx(obs);
    j_adm[1] = obs->Jy(obs);
    j_adm[2] = obs->Jz(obs);
    
    printf("|--> After CM update P_ADM = (%e,%e,%e)\n",p_adm[0],p_adm[1],p_adm[2]);
    printf("|--> After CM update J_ADM = (%e,%e,%e)\n",j_adm[0],j_adm[1],j_adm[2]);
    
    free_observable(obs);
  }
  else
  {
    printf("|--> |Py_ADM2 - Py_ADM1|/|Py_ADM2| = %g <= %g\n"
           "     |--> no x_CM update.\n",dPy_Py,dP);
  }
  
  free(freedata_grid);
  free(freedata_patch);
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

/* adjust the apparent horizon radius to acquire the desired BH irreducible mass */
static void adjust_AH_radius(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  pr_line_custom('=');
  printf("{ Adjusting apparent horizon radius to meet BH mass ...\n");
  
  const double target_bh_mass  = Pgetd("BH_irreducible_mass");
  const double current_r_excision = Pgetd("r_excision");
  const double W  = Pgetd("BH_r_excision_update_weight");
  const double dM_tolerance = Pgetd("BH_mass_tolerance");
  const double irr_mass     = bbn_BH_irreducible_mass(grid);
  double kommar_mass, adm_mass;
  Observable_T *obs = 0;
  double dr, r_excision,dM;
  
  obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "ADM(M)|BH";
  plan_observable(obs);
  adm_mass = obs->M(obs);
  free_observable(obs);
  
  obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity = "Kommar(M)|BH";
  plan_observable(obs);
  kommar_mass = obs->M(obs);
  free_observable(obs);
  
  printf("|--> current BH irreducible mass = %e\n",irr_mass);
  printf("|--> current BH ADM mass         = %e\n",adm_mass);
  printf("|--> current BH Kommar mass      = %e\n",kommar_mass);
  
  Psetd("BH_irreducible_mass_current",irr_mass);
  Psetd("BH_ADM_mass",adm_mass);
  Psetd("BH_Kommar_mass",kommar_mass);
  
  dM = fabs(irr_mass/target_bh_mass-1);
  dr = -current_r_excision*(irr_mass/target_bh_mass-1);
  if (EQL(W,0))
  {
    dr = 0;
    printf("|--> updating weight factor is zero.\n");
  }
  if (LSSEQL(dM,dM_tolerance)) 
  {
    dr = 0;
    printf("|--> |dM/M| = %g < Tol. = %g\n",dM,dM_tolerance);
  }

  r_excision = current_r_excision + W*dr;
  
  GridParams->R_BH_r    = r_excision;
  GridParams->BH_R_type = "PerfectSphere";
  
  Psetd("r_excision",r_excision);
  
  if (EQL(dr,0))/* => no change in AH surface */
    Pseti("did_AH_surface_change?",0);
  else          /* => change in AH surface */
    Pseti("did_AH_surface_change?",1);
  
  printf("} Adjusting apparent horizon radius to meet BH mass ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* adjust the Omega_BH to acquire the desired BH spin */
static void adjust_BH_Omega(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  pr_line_custom('=');
  printf("{ Adjusting BH Omega ...\n");
  
  const double irr_massc = Pgetd("BH_irreducible_mass_current");
  const double irr_mass  = Pgetd("BH_irreducible_mass");
  const double irr_mc2   = Pow2(irr_massc);
  const double W         = Pgetd("BH_spin_update_weight");
  const double chi_xt    = Pgetd("BH_chi_U0");
  const double chi_yt    = Pgetd("BH_chi_U1");
  const double chi_zt    = Pgetd("BH_chi_U2");
  const double Omega_x   = Pgetd("BH_Omega_U0");
  const double Omega_y   = Pgetd("BH_Omega_U1");
  const double Omega_z   = Pgetd("BH_Omega_U2");
  double chr_mass = 0;
  double S_BH[3] = {0},s_BH2;
  double dOmega_x = 0,dOmega_y = 0,dOmega_z = 0;
  double chi_xc,chi_yc,chi_zc;
  
  /* find current BH spin */
  bbn_define_spin_JRP(S_BH,grid,"BH");
  
  s_BH2    = Pow2(S_BH[0])+Pow2(S_BH[1])+Pow2(S_BH[2]);
  chr_mass = sqrt(irr_mc2+s_BH2/(4*irr_mc2));
  chi_xc   = S_BH[0]/Pow2(chr_mass);
  chi_yc   = S_BH[1]/Pow2(chr_mass);
  chi_zc   = S_BH[2]/Pow2(chr_mass);
  
  /* adjust Omegas */
  if (!EQL(chi_xt,0))
    dOmega_x = -(chi_xc-chi_xt)/(4*chr_mass) + 
            (irr_massc-irr_mass)/(4*Pow2(irr_massc))*chi_xc;
  if (!EQL(chi_yt,0))
    dOmega_y = -(chi_yc-chi_yt)/(4*chr_mass) + 
            (irr_massc-irr_mass)/(4*Pow2(irr_massc))*chi_yc;
  if (!EQL(chi_zt,0))
    dOmega_z = -(chi_zc-chi_zt)/(4*chr_mass) + 
            (irr_massc-irr_mass)/(4*Pow2(irr_massc))*chi_zc;
 
  if (EQL(W,0))
  {
    printf("|--> updating weight factor is zero.\n");
  }
  else
  {
    Psetd("BH_omega_U0",Omega_x+W*dOmega_x);
    Psetd("BH_omega_U1",Omega_y+W*dOmega_y);
    Psetd("BH_omega_U2",Omega_z+W*dOmega_z);
    Psetd("BH_Christodoulou_mass",chr_mass);
  }
  
  GridParams->a_BH   = DBL_MAX;/* to catch error */
  printf("} Adjusting BH Omega ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* root finder eqution for Euler equation constant */
static double Euler_eq_const_rootfinder_eq(void *params,const double *const x)
{
  struct Euler_eq_const_RootFinder_S *const par = params;
  
  return bbn_NS_baryonic_mass(par->grid,x[0]) - par->NS_baryonic_mass;
}

/* use previous grid to interpolate values of the fields that will be solved for the next grid */
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev)
{
  pr_line_custom('=');
  printf("{ Interpolating & initializing to the next grid ...\n");
  
  /* if it is ready */
  if (Pgeti("use_previous_data"))
  {
    printf("~> Using the fields of the previous grid.\n");
    printf("} Interpolating & initializing to the next grid ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
    return;
  }
  
  const unsigned np = grid_next->np;
  const int change_res_flg = Pgeti("did_resolution_change?");
  const int change_NS_flg  = Pgeti("did_NS_surface_change?");
  const int change_AH_flg  = Pgeti("did_AH_surface_change?");
 
  unsigned p;
 
  /* the following fields are interpolated: */
  /* B0_U[0-2],psi,eta,phi,enthalpy */
  
  /* to avoid race condition between threads write all coeffs */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid_prev->np; ++p)/* note: grid_prev has more patches! */
  {
    Patch_T *patch = grid_prev->patch[p];
    Field_T *R1_f  = 0;
    Field_T *R2_f  = 0;
    
    /* surface fields also are used for the interpolation in X_of_x function */
    if (patch->coordsys == CubedSpherical)
    {
      R1_f = patch->CoordSysInfo->CubedSphericalCoord->R1_f;
      R2_f = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
      if (R1_f)
        make_coeffs_2d(R1_f,0,1);/* X and Y direction */
      if (R2_f)
        make_coeffs_2d(R2_f,0,1);/* X and Y direction */
    }
    else if (patch->coordsys == Cartesian)
    {
      R1_f  = R2_f = 0;
    }
    else
      Error0(NO_OPTION);
      
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
    
    if (IsItNSPatch(patch) || IsItNSSurroundingPatch(patch))
    {
      DECLARE_FIELD(phi)
      DECLARE_FIELD(enthalpy)
      make_coeffs_3d(phi);
      make_coeffs_3d(enthalpy);
    }
  }
  
  /* ORDER IS IMPORTANT */
  if (change_res_flg)/* make from scratch */
  {
    printf("~> resolution is changed:\n");
    
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
      
      REALLOC_v_WRITE_v(B0_U0)
      REALLOC_v_WRITE_v(B0_U1)
      REALLOC_v_WRITE_v(B0_U2)
      REALLOC_v_WRITE_v(psi)
      REALLOC_v_WRITE_v(eta)
      
      if (IsItNSPatch(patch))
      {
        REALLOC_v_WRITE_v(phi)
        REALLOC_v_WRITE_v(enthalpy)
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
        
        printf("|--> %s:\n"
               "     |--> interpolating B0_U0    ~> Done.\n"
               "     |--> interpolating B0_U1    ~> Done.\n"
               "     |--> interpolating B0_U2    ~> Done.\n"
               "     |--> interpolating psi      ~> Done.\n"
               "     |--> interpolating eta      ~> Done.\n"
               "     |--> interpolating phi      ~> Done.\n"
               "     |--> interpolating enthalpy ~> Done.\n"
               ,patch->name);
        fflush(stdout);
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
        
        printf("|--> %s:\n"
               "     |--> interpolating B0_U0    ~> Done.\n"
               "     |--> interpolating B0_U1    ~> Done.\n"
               "     |--> interpolating B0_U2    ~> Done.\n"
               "     |--> interpolating psi      ~> Done.\n"
               "     |--> interpolating eta      ~> Done.\n"
               ,patch->name);
        fflush(stdout);
      }
    }/* end of for (p = 0; p < np; ++p) */
  }/* end of else if (change_res_flg) */
  else if (!change_NS_flg)/* only NS, filling box and outermost are reusable */
  {
    printf("~> NS, filling box and outermost patches are the same:\n");
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < np; ++p)
    {
      Patch_T *patch = grid_next->patch[p];
      unsigned nn    = patch->nn;
      char *root_name;
      
      root_name = strstr(patch->name,"_");/* the patch->name convention is grid\d?_root */
      assert(root_name);
      root_name++;
      
      if (IsItNSPatch(patch))
      {
        Patch_T *patchp = GetPatch(root_name,grid_prev);
        unsigned ijk;
        
        prep_and_call(phi)
        prep_and_call(enthalpy)
        prep_and_call(B0_U0)
        prep_and_call(B0_U1)
        prep_and_call(B0_U2)
        prep_and_call(psi)
        prep_and_call(eta)
        
        for (ijk = 0; ijk < nn; ++ijk)
        {
          copy_values(phi)
          copy_values(enthalpy)
          copy_values(B0_U0)
          copy_values(B0_U1)
          copy_values(B0_U2)
          copy_values(psi)
          copy_values(eta)
        }
        
        printf("|--> %s:\n"
               "     |--> copying B0_U0    ~> Done.\n"
               "     |--> copying B0_U1    ~> Done.\n"
               "     |--> copying B0_U2    ~> Done.\n"
               "     |--> copying psi      ~> Done.\n"
               "     |--> copying eta      ~> Done.\n"
               "     |--> copying phi      ~> Done.\n"
               "     |--> copying enthalpy ~> Done.\n"
               ,patch->name);
        fflush(stdout);
        
      }
      else if (
               IsItNSSurroundingPatch(patch) ||
               IsItOutermostPatch(patch)     ||
               IsItFillingBoxPatch(patch)
               )
      {
        Patch_T *patchp = GetPatch(root_name,grid_prev);
        unsigned ijk;
        
        prep_and_call(B0_U0)
        prep_and_call(B0_U1)
        prep_and_call(B0_U2)
        prep_and_call(psi)
        prep_and_call(eta)
        
        for (ijk = 0; ijk < nn; ++ijk)
        {
          copy_values(B0_U0)
          copy_values(B0_U1)
          copy_values(B0_U2)
          copy_values(psi)
          copy_values(eta)
        }
        printf("|--> %s:\n"
               "     |--> copying B0_U0    ~> Done.\n"
               "     |--> copying B0_U1    ~> Done.\n"
               "     |--> copying B0_U2    ~> Done.\n"
               "     |--> copying psi      ~> Done.\n"
               "     |--> copying eta      ~> Done.\n"
               ,patch->name);
        fflush(stdout);
         
      }
      else
      {
        continue;
      }
    }/* end of for (p = 0; p < np; ++p) */
    
    Patch_T *BH_patches[6] = {0};/* 6 cubed spherical */
    unsigned count = 0;
    for (p = 0; p < np; ++p)
    {
      Patch_T *patch = grid_next->patch[p];
      if (IsItHorizonPatch(patch))
      {
        BH_patches[count] = patch;
        count++;
      }
    }
    assert(count == 6);
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < 6; ++p)
    {
      Patch_T *patch = BH_patches[p];
      unsigned nn = patch->nn;
      char hint[100],*root_name;
      unsigned ijk;
      
      root_name = strstr(patch->name,"_");/* the patch->name convention is grid\d?_root */
      assert(root_name);
      root_name++;
      sprintf(hint,"%s",root_name);
      
      REALLOC_v_WRITE_v(B0_U0)
      REALLOC_v_WRITE_v(B0_U1)
      REALLOC_v_WRITE_v(B0_U2)
      REALLOC_v_WRITE_v(psi)
      REALLOC_v_WRITE_v(eta)
      
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
      
      printf("|--> %s:\n"
             "     |--> interpolating B0_U0    ~> Done.\n"
             "     |--> interpolating B0_U1    ~> Done.\n"
             "     |--> interpolating B0_U2    ~> Done.\n"
             "     |--> interpolating psi      ~> Done.\n"
             "     |--> interpolating eta      ~> Done.\n"
             ,patch->name);
      fflush(stdout);
    
    }/* end of for (p = 0; p < 6; ++p) */
  }/* end of else if (!change_NS_flg) */
  else if (!change_AH_flg)/* only BH, filling box and outermost are reusable */
  {
    printf("~> BH, filling box and outermost patches are the same:\n");
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < np; ++p)
    {
      Patch_T *patch = grid_next->patch[p];
      unsigned nn    = patch->nn;
      char *root_name;
      
      root_name = strstr(patch->name,"_");/* the patch->name convention is grid\d?_root */
      assert(root_name);
      root_name++;
      
      if (
          IsItHorizonPatch(patch)    ||
          IsItOutermostPatch(patch)  ||
          IsItFillingBoxPatch(patch)
         )
      {
        Patch_T *patchp = GetPatch(root_name,grid_prev);
        unsigned ijk;
        
        prep_and_call(B0_U0)
        prep_and_call(B0_U1)
        prep_and_call(B0_U2)
        prep_and_call(psi)
        prep_and_call(eta)
        
        for (ijk = 0; ijk < nn; ++ijk)
        {
          copy_values(B0_U0)
          copy_values(B0_U1)
          copy_values(B0_U2)
          copy_values(psi)
          copy_values(eta)
        }
        printf("|--> %s:\n"
               "     |--> copying B0_U0    ~> Done.\n"
               "     |--> copying B0_U1    ~> Done.\n"
               "     |--> copying B0_U2    ~> Done.\n"
               "     |--> copying psi      ~> Done.\n"
               "     |--> copying eta      ~> Done.\n"
               ,patch->name);
        fflush(stdout);
         
      }
      else
      {
        continue;
      }
    }/* end of for (p = 0; p < np; ++p) */
    
    /* interpolating from the NS and surrounding patches: */
    Patch_T *NS_patches[13] = {0};/* 2*6 cubed spherical + 1 box*/
    unsigned count = 0;
    for (p = 0; p < np; ++p)
    {
      Patch_T *patch = grid_next->patch[p];
      if (
          IsItNSPatch(patch)           ||
          IsItNSSurroundingPatch(patch)
          )
      {
        NS_patches[count] = patch;
        count++;
      }
    }
    assert(count == 13);
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < 13; ++p)
    {
      Patch_T *patch = NS_patches[p];
      unsigned nn = patch->nn;
      char hint[100],*root_name;
      unsigned ijk;
      
      root_name = strstr(patch->name,"_");/* the patch->name convention is grid\d?_root */
      assert(root_name);
      root_name++;
      sprintf(hint,"%s",root_name);
      
      REALLOC_v_WRITE_v(B0_U0)
      REALLOC_v_WRITE_v(B0_U1)
      REALLOC_v_WRITE_v(B0_U2)
      REALLOC_v_WRITE_v(psi)
      REALLOC_v_WRITE_v(eta)
      
      if (IsItNSPatch(patch))
      {
        REALLOC_v_WRITE_v(phi)
        REALLOC_v_WRITE_v(enthalpy)
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
        
        printf("|--> %s:\n"
               "     |--> interpolating B0_U0    ~> Done.\n"
               "     |--> interpolating B0_U1    ~> Done.\n"
               "     |--> interpolating B0_U2    ~> Done.\n"
               "     |--> interpolating psi      ~> Done.\n"
               "     |--> interpolating eta      ~> Done.\n"
               "     |--> interpolating phi      ~> Done.\n"
               "     |--> interpolating enthalpy ~> Done.\n"
               ,patch->name);
        fflush(stdout);
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
        
        printf("|--> %s:\n"
               "     |--> interpolating B0_U0    ~> Done.\n"
               "     |--> interpolating B0_U1    ~> Done.\n"
               "     |--> interpolating B0_U2    ~> Done.\n"
               "     |--> interpolating psi      ~> Done.\n"
               "     |--> interpolating eta      ~> Done.\n"
               ,patch->name);
        fflush(stdout);
      }
    }/* end of for (p = 0; p < 6; ++p) */
  }/* end of else if (!change_AH_flg) */
  else/* if both NS and AH are changed but the resolution */
  {
    printf("~> Filling box and outermost patches are the same:\n");
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < np; ++p)
    {
      Patch_T *patch = grid_next->patch[p];
      unsigned nn    = patch->nn;
      char *root_name;
      
      root_name = strstr(patch->name,"_");/* the patch->name convention is grid\d?_root */
      assert(root_name);
      root_name++;
      
      if (
          IsItOutermostPatch(patch)  ||
          IsItFillingBoxPatch(patch)
         )
      {
        Patch_T *patchp = GetPatch(root_name,grid_prev);
        unsigned ijk;
      
        prep_and_call(B0_U0)
        prep_and_call(B0_U1)
        prep_and_call(B0_U2)
        prep_and_call(psi)
        prep_and_call(eta)
        
        for (ijk = 0; ijk < nn; ++ijk)
        {
          copy_values(B0_U0)
          copy_values(B0_U1)
          copy_values(B0_U2)
          copy_values(psi)
          copy_values(eta)
        }
        printf("|--> %s:\n"
               "     |--> copying B0_U0    ~> Done.\n"
               "     |--> copying B0_U1    ~> Done.\n"
               "     |--> copying B0_U2    ~> Done.\n"
               "     |--> copying psi      ~> Done.\n"
               "     |--> copying eta      ~> Done.\n"
               ,patch->name);
        fflush(stdout);
         
      }
      else
      {
        continue;
      }
    }/* end of for (p = 0; p < np; ++p) */
    
    /* interpolating from the NS, BH surrounding patches: */
    Patch_T *NS_BH_patches[19] = {0};/* 3*6 cubed spherical + 1 box*/
    unsigned count = 0;
    for (p = 0; p < np; ++p)
    {
      Patch_T *patch = grid_next->patch[p];
      if (
          IsItNSPatch(patch)           ||
          IsItHorizonPatch(patch)      ||
          IsItNSSurroundingPatch(patch)
          
          )
      {
        NS_BH_patches[count] = patch;
        count++;
      }
    }
    assert(count == 19);
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < 19; ++p)
    {
      Patch_T *patch = NS_BH_patches[p];
      unsigned nn = patch->nn;
      char hint[100],*root_name;
      unsigned ijk;
      
      root_name = strstr(patch->name,"_");/* the patch->name convention is grid\d?_root */
      assert(root_name);
      root_name++;
      sprintf(hint,"%s",root_name);
      
      REALLOC_v_WRITE_v(B0_U0)
      REALLOC_v_WRITE_v(B0_U1)
      REALLOC_v_WRITE_v(B0_U2)
      REALLOC_v_WRITE_v(psi)
      REALLOC_v_WRITE_v(eta)
      
      if (IsItNSPatch(patch))
      {
        REALLOC_v_WRITE_v(phi)
        REALLOC_v_WRITE_v(enthalpy)
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
        
        printf("|--> %s:\n"
               "     |--> interpolating B0_U0    ~> Done.\n"
               "     |--> interpolating B0_U1    ~> Done.\n"
               "     |--> interpolating B0_U2    ~> Done.\n"
               "     |--> interpolating psi      ~> Done.\n"
               "     |--> interpolating eta      ~> Done.\n"
               "     |--> interpolating phi      ~> Done.\n"
               "     |--> interpolating enthalpy ~> Done.\n"
               ,patch->name);
        fflush(stdout);
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
        
        printf("|--> %s:\n"
               "     |--> interpolating B0_U0    ~> Done.\n"
               "     |--> interpolating B0_U1    ~> Done.\n"
               "     |--> interpolating B0_U2    ~> Done.\n"
               "     |--> interpolating psi      ~> Done.\n"
               "     |--> interpolating eta      ~> Done.\n"
               ,patch->name);
        fflush(stdout);
      }
    }/* end of for (p = 0; p < 19; ++p) */
  }/* end of else */
  
  /* initializing some other fields: */
  /* W_U[0-2],Beta_U[0-2],B1_U[0-2] */
  const double Omega_NS_x = Pgetd("NS_Omega_U0");
  const double Omega_NS_y = Pgetd("NS_Omega_U1");
  const double Omega_NS_z = Pgetd("NS_Omega_U2");
  const double C_NSx = Pgetd("NS_Center_x");
  const double C_NSy = Pgetd("NS_Center_y");
  const double C_NSz = Pgetd("NS_Center_z");
  
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
      REALLOC_v_WRITE_v(W_U0)
      REALLOC_v_WRITE_v(W_U1)
      REALLOC_v_WRITE_v(W_U2)
      
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double x = patch->node[ijk]->x[0]-C_NSx;
        double y = patch->node[ijk]->x[1]-C_NSy;
        double z = patch->node[ijk]->x[2]-C_NSz;
        
        /* spin part */
        W_U0[ijk] = Omega_NS_y*z-Omega_NS_z*y;
        W_U1[ijk] = Omega_NS_z*x-Omega_NS_x*z;
        W_U2[ijk] = Omega_NS_x*y-Omega_NS_y*x;
      }
    }/* end of if (IsItNSPatch(patch)) */
    
  }/* end of for (p = 0; p < np; ++p) */
  
  printf(
         "|--> initializing W_U0    ~> Done.\n"
         "|--> initializing W_U1    ~> Done.\n"
         "|--> initializing W_U2    ~> Done.\n"
         "|--> initializing Beta_U0 ~> Done.\n"
         "|--> initializing Beta_U1 ~> Done.\n"
         "|--> initializing Beta_U2 ~> Done.\n"
         "|--> initializing B1_U2   ~> Done.\n"
         "|--> initializing B1_U2   ~> Done.\n"
         "|--> initializing B1_U2   ~> Done.\n");
  fflush(stdout);
    
  printf("} Interpolating & initializing to the next grid ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* given field name, X and patch, finds the value of the field in X  
// using interpolation.
// ->return value: f(X) */
static double interpolate_from_patch_prim(const char *const field,const double *const X,Patch_T *const patch)
{
  double interp;
  Interpolation_T *interp_s = init_interpolation();
  Field_T *const F_field    = patch->pool[Ind(field)];
  
  interp_s->field = F_field;
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
  Needle_T *needle = alloc_needle();
  const double LOW_RES_ERR = 1E-9;
  unsigned *found;
  unsigned p;
  
  needle->grid = grid;
  needle->x    = x;
  
  /* find this point everywhere */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    needle_in(needle,patch);
      
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
          Error0("Point not found!\n");
        }
      }
      else
      {
        fprintf(stderr,"The point (%g,%g,%g) could not be found!\n",x[0],x[1],x[2]);
        Error0("Point not found!\n");
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
  const unsigned lmax = (unsigned)Pgeti("NS_surface_Ylm_expansion_max_l");
  const double RESIDUAL = sqrt(Pgetd("RootFinder_Tolerance"));
  const double max_h_L2_res = Pgetd("NS_enthalpy_allowed_residual");
  double h_L2_res = 0;
  //const double W1  = Pgetd("Solving_Field_Update_Weight");
  //const double W2  = 1-W1;
  double theta,phi;
  double *Rnew_NS = 0;/* new R for NS */
  double Max_R_NS = 0;/* maximum radius of NS */
  double Min_R_NS = DBL_MAX;/* minimum radius of NS */
  double guess    = 0;
  double *h_res   = 0;/* residual of h */
  double X[3],x[3],N[3];
  char stem[1000],*affix;
  int NS_surface_finder_work_flg = 1;/* whether surface finder worked or not */
  unsigned i,j;
  unsigned l,m;
  
  /* populate root finder */
  Root_Finder_T *root = init_root_finder(1);
  root->type      = Pgets("RootFinder_Method");
  root->tolerance = Pgetd("RootFinder_Tolerance");
  root->MaxIter   = (unsigned)Pgeti("RootFinder_Max_Number_of_Iteration");
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
  h_res   = alloc_double(Ntheta*Nphi);
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
      h_res[ij(i,j)] = root->residual;
      /* if root finder is not OK for some reason */
      if (GRT(root->residual,RESIDUAL))
      {
        printf(". Root finder for NS surface at %s:\n.. ",h_patch->name);
        print_root_finder_exit_status(root);
        printf(".. Residual = %g\n",root->residual);
        //NS_surface_finder_work_flg = 0;/* since the residual is large don't update */
      }
      
      /*  new coords of R respect to the center of NS */
      y2[0] += N[0]*dr[0];
      y2[1] += N[1]*dr[0];
      y2[2] += N[2]*dr[0];
      Rnew_NS[ij(i,j)] = root_square(3,y2,0);
      free(dr);
      
      /* find the max NS radius */
      if (Rnew_NS[ij(i,j)] > Max_R_NS)
        Max_R_NS = Rnew_NS[ij(i,j)];
      /* find the min NS radius */
      if (Rnew_NS[ij(i,j)] < Min_R_NS)
        Min_R_NS = Rnew_NS[ij(i,j)];
    }/* end of for (j = 0; j < Nphi; ++j) */
  }/* end of for (i = 0; i < Ntheta; ++i) */
  
  h_L2_res = L2_norm(Ntheta*Nphi,h_res,0);
  if (h_L2_res > max_h_L2_res)
    NS_surface_finder_work_flg = 0;/* since the residual is large don't update */
  else
    NS_surface_finder_work_flg = 1;
    
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
  
  /* printing */
  printf("|--> Max NS radius       = %e\n",Max_R_NS);
  printf("|--> Min NS radius       = %e\n",Min_R_NS);
  printf("|--> L2 norm of enthalpy = %e\n",h_L2_res);
  l = lmax;
  for (m = 0; m <= l; ++m)
  {
    unsigned lm = lm2n(l,m);
    printf("|--> Truncation error [Real(C[%u][%u])] = %e\n",l,m,realClm[lm]);
    printf("|--> Truncation error [Imag(C[%u][%u])] = %e\n",l,m,imagClm[lm]);
  }
  
  /* if some day you wanna filter Clm's */
  if (0)
  {
    const double e = 0.1;
    for (l = 0; l <= lmax; ++l)
      for (m = 0; m <= l; ++m)
      {
        unsigned lm = lm2n(l,m);
        realClm[lm] /= (1+e*Pow2(l)*Pow2(l+1));
        imagClm[lm] /= (1+e*Pow2(l)*Pow2(l+1));
      }
  }
  
  free(Rnew_NS);
  free(h_res);
  free_root_finder(root);
  
  Pseti("did_NS_surface_finder_work?",NS_surface_finder_work_flg);
  
  Psetd("NS_max_radius",Max_R_NS);
  Psetd("NS_min_radius",Min_R_NS);
  
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
  const double tan_phi2   = Pow2(tan_phi);
  const double cos_theta2 = Pow2(cos_theta);
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
    double a = 0, b = 0;
    double a_sign = 0,b_sign = 0,c_sign = 0;
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
        Error0(NO_OPTION);
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
        Error0(NO_OPTION);
    }
    
    X[0] = fabs(a)*a_sign;
    X[1] = fabs(b)*b_sign;
    
    /* check if x of X really gives you the correct angles */
    x_of_X(x,X,patch);
    x[0] -= c[0];
    x[1] -= c[1];
    x[2] -= c[2];
    r = root_square(3,x,0);
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
    Error0("(X,Y,Z) or patch could not be found.\n");
}


/* given theta, phi and knowing the fact that they are on BH surface, 
// it finds the corresponding patch and X,Y,Z coordinate. */
static void find_XYZ_and_patch_of_theta_phi_BH_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid)
{
  const double tan_phi    = tan(phi);
  const double cos_theta  = cos(theta);
  const double tan_phi2   = Pow2(tan_phi);
  const double cos_theta2 = Pow2(cos_theta);
  Flag_T found_flg = NO;
  unsigned p;
  
  X[2] = 0;/* since we are on BH surface from BH surrounding side */
  
  /* check all of BH patches in which (x,y,z) and 
  // (X,Y,Z) and (theta,phi) are consistent */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItHorizonPatch(patch))
      continue;

    Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
    const double *c = patch->c;
    double a = 0, b = 0;
    double a_sign = 0,b_sign = 0,c_sign = 0;
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
        Error0(NO_OPTION);
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
        Error0(NO_OPTION);
    }
    
    X[0] = fabs(a)*a_sign;
    X[1] = fabs(b)*b_sign;
    
    /* check if x of X really gives you the correct angles */
    x_of_X(x,X,patch);
    x[0] -= c[0];
    x[1] -= c[1];
    x[2] -= c[2];
    r = root_square(3,x,0);
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
    Error0("(X,Y,Z) or patch could not be found.\n");
}

/* make patches inside the excision region of BH and and extrapolate
// metric fields i.e. beta,eta and psi inside this region */
void bbn_extrapolate_metric_fields_insideBH(Grid_T *const grid)
{
  if (Pgeti("STOP"))
    return;
    
  pr_line_custom('=');
  printf("{ Extrapolating metric fields inside the BH ...\n");
  
  unsigned p;
  
  /* add patches in side the excision region */
  add_patches_insideBH(grid);
  
  if (Pcmps("extrapolate_inside_BH_method","Ylm"))
  {
    /* add and update necessary fields inside the BH */
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      if (!IsItInsideBHPatch(patch))
        continue;
      
      bbn_add_fields_in_patch(patch);
      
      /* making B1 */
      //bbn_update_B1_U012(patch);
    }
    
    /* extrapolate the fields inside the BH */
    extrapolate_insideBH_CS_C0_Ylm(grid,"B0_U0");
    extrapolate_insideBH_CS_C0_Ylm(grid,"B0_U1");
    extrapolate_insideBH_CS_C0_Ylm(grid,"B0_U2");
    extrapolate_insideBH_CS_C0_Ylm(grid,"psi");
    extrapolate_insideBH_CS_C0_Ylm(grid,"eta");
  }
  else if (Pcmps("extrapolate_inside_BH_method","linear"))
  {
    extrapolate_insideBH_CS_linear(grid);
  }
  else if (Pcmps("extrapolate_inside_BH_method","WTGR"))
  {
    extrapolate_insideBH_CS_WTGR(grid);
  }
  else
    Error0(NO_OPTION);
  
  printf("} Extrapolating metric fields inside the BH ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
} 

/* add patches inside the excision region */
static void add_patches_insideBH(Grid_T *const grid)
{
  if (!strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
    Error0(NO_OPTION);
    
  const unsigned np1 = grid->np;
  const unsigned np2 = np1+7;/* 6 cubed spherical + 1 box */
  Grid_T *bh_grid = 0;
  unsigned i;
  
  /* allocating */
  grid->patch = realloc(grid->patch,(np2+1)*sizeof(*grid->patch));
  IsNull(grid->patch);
  grid->patch[np2] = 0;
  
  for (i = np1; i < np2; i++)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    IsNull(grid->patch[i]);
  }
  grid->np = np2;
  
  /* populate the patches */
  populate_right_BH_central_box(grid,np1);
  populate_right_BH(grid,np1+1);
  
  bh_grid = calloc(1,sizeof(*bh_grid));
  IsNull(bh_grid);
  
  bh_grid->patch = &grid->patch[np1];
  alloc_nodes(bh_grid);
  make_nodes(bh_grid);
  make_JacobianT(bh_grid);
  /* test printing coords */
  if (test_print(PRINT_COORDS))
    pr_coords(grid);
    
  free(bh_grid);
}

/* extrapolate the fields B0,B1,eta and psi inside the BH.
// we demand the fields have C^1 continuity across the apparent horizon.
// we assume: f = a*(r-rh)+b, r is coordinate distance to the center of BH. 
// note: at the central box, we put B0 to 0 and psi and eta to 1 and
// B1 is calculated fromm its formula "Omega cross r". */
static void extrapolate_insideBH_CS_linear(Grid_T *const grid)
{
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if (!IsItInsideBHPatch(patch))
      continue;
      
    Patch_T *BHsur_patch;/* corresponding bh surrounding patch */
    const unsigned *n = patch->n;
    double x[3],X[3],N[3];
    double rh,r,R;
    double a_BU0,b_BU0,
           a_BU1,b_BU1,
           a_BU2,b_BU2,
           a_eta,b_eta,
           a_psi,b_psi;
    unsigned ijk,i,j,k;
 
    /* add fields: */
    bbn_add_fields_in_patch(patch);
    
    REALLOC_v_WRITE_v(B0_U0)
    REALLOC_v_WRITE_v(B0_U1)
    REALLOC_v_WRITE_v(B0_U2)
    REALLOC_v_WRITE_v(psi)
    REALLOC_v_WRITE_v(eta)
    
    /* making B1 */
    bbn_update_B1_U012(patch);
    
    /* for the central box we have: */
    if (strstr(patch->name,"right_central_box"))
    {
      unsigned nn = patch->nn;
      for (ijk = 0; ijk < nn; ++ijk)
      {
        psi[ijk] = 1;
        eta[ijk] = 1;
        B0_U0[ijk] = B0_U1[ijk] = B0_U2[ijk] = 0;
      }
      continue;
    }
    
    /* find the corresponding BH surrounding patch to be used for extrapolation */
    char stem[1000];
    char *affix = regex_find("_[[:alpha:]]{2,5}$",patch->name);/* finding the side of the patch */
    assert(affix);
    sprintf(stem,"right_BH_surrounding%s",affix);
    free(affix);
    BHsur_patch = GetPatch(stem,grid);
    
    /* populate B0, eta and psi */
    for (i = 0; i < n[0]; ++i)
    {
      for (j = 0; j < n[1]; ++j)
      {
        /* calculate rh at the horizon patch */
        ijk  = L(n,i,j,n[2]-1);
        x[0] = patch->node[ijk]->x[0]-patch->c[0];
        x[1] = patch->node[ijk]->x[1]-patch->c[1];
        x[2] = patch->node[ijk]->x[2]-patch->c[2]; 
        rh   = root_square(3,x,0);
        
        X_of_x(X,patch->node[ijk]->x,BHsur_patch);
        
        N[0] = interpolate_from_patch_prim("_HS_U0",X,BHsur_patch);
        N[1] = interpolate_from_patch_prim("_HS_U1",X,BHsur_patch);
        N[2] = interpolate_from_patch_prim("_HS_U2",X,BHsur_patch);
        
        double B0_U0_i    = interpolate_from_patch_prim("B0_U0",X,BHsur_patch);
        double dB0_U0D0_i = interpolate_from_patch_prim("dB0_U0D0",X,BHsur_patch);
        double dB0_U0D1_i = interpolate_from_patch_prim("dB0_U0D1",X,BHsur_patch);
        double dB0_U0D2_i = interpolate_from_patch_prim("dB0_U0D2",X,BHsur_patch);
        
        double B0_U1_i    = interpolate_from_patch_prim("B0_U1",X,BHsur_patch);
        double dB0_U1D0_i = interpolate_from_patch_prim("dB0_U1D0",X,BHsur_patch);
        double dB0_U1D1_i = interpolate_from_patch_prim("dB0_U1D1",X,BHsur_patch);
        double dB0_U1D2_i = interpolate_from_patch_prim("dB0_U1D2",X,BHsur_patch);
        
        double B0_U2_i    = interpolate_from_patch_prim("B0_U2",X,BHsur_patch);
        double dB0_U2D0_i = interpolate_from_patch_prim("dB0_U2D0",X,BHsur_patch);
        double dB0_U2D1_i = interpolate_from_patch_prim("dB0_U2D1",X,BHsur_patch);
        double dB0_U2D2_i = interpolate_from_patch_prim("dB0_U2D2",X,BHsur_patch);
        
        double eta_i     = interpolate_from_patch_prim("eta",X,BHsur_patch);
        double deta_D0_i = interpolate_from_patch_prim("deta_D0",X,BHsur_patch);
        double deta_D1_i = interpolate_from_patch_prim("deta_D1",X,BHsur_patch);
        double deta_D2_i = interpolate_from_patch_prim("deta_D2",X,BHsur_patch);
        
        double psi_i     = interpolate_from_patch_prim("psi",X,BHsur_patch);
        double dpsi_D0_i = interpolate_from_patch_prim("dpsi_D0",X,BHsur_patch);
        double dpsi_D1_i = interpolate_from_patch_prim("dpsi_D1",X,BHsur_patch);
        double dpsi_D2_i = interpolate_from_patch_prim("dpsi_D2",X,BHsur_patch);
        
        a_BU0 = (N[0]*dB0_U0D0_i+N[1]*dB0_U0D1_i+N[2]*dB0_U0D2_i);
        b_BU0 = B0_U0_i;
        
        a_BU1 = (N[0]*dB0_U1D0_i+N[1]*dB0_U1D1_i+N[2]*dB0_U1D2_i);
        b_BU1 = B0_U1_i;
        
        a_BU2 = (N[0]*dB0_U2D0_i+N[1]*dB0_U2D1_i+N[2]*dB0_U2D2_i);
        b_BU2 = B0_U2_i;
        
        a_eta = (N[0]*deta_D0_i+N[1]*deta_D1_i+N[2]*deta_D2_i);
        b_eta = eta_i;
        
        a_psi = (N[0]*dpsi_D0_i+N[1]*dpsi_D1_i+N[2]*dpsi_D2_i);
        b_psi = psi_i;
        
        for (k = 0; k < n[2]; ++k)
        {
          ijk = L(n,i,j,k);
          x[0] = patch->node[ijk]->x[0]-patch->c[0];
          x[1] = patch->node[ijk]->x[1]-patch->c[1];
          x[2] = patch->node[ijk]->x[2]-patch->c[2];
          
          r = root_square(3,x,0);
          R = r-rh;
          eta[ijk]   = a_eta*R+b_eta;
          psi[ijk]   = b_psi*exp(R/rh);/* "a_psi*R+b_psi" not good! */
          UNUSED(a_psi);
          B0_U0[ijk] = a_BU0*R+b_BU0;
          B0_U1[ijk] = a_BU1*R+b_BU1;
          B0_U2[ijk] = a_BU2*R+b_BU2;
          
        }/* end of for (k = 0 ; k < n[2]; ++k) */
      }/* end of for (j = 0; j < n[1]; ++j) */
    }/* end of for (i = 0; i < n[0]; ++i) */
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
}

/* for BAM initial data reader,
// extrapolate the fields Beta,eta, psi, _gamma's
// inside the BH, using the method developed by Wolfgang and Geroge,
// more info "http://fau.digital.flvc.org/islandora/object/fau%3A4224". */
static void extrapolate_insideBH_CS_WTGR(Grid_T *const grid)
{
  printf("|--> BH-filler method = WTGR.\n");
  fflush(stdout);
  
  const double EPS            = 1E-12;/* to avoid division by zero */
  const double EPS2           = 1E-6;/* to increase r_fill radius a bit */
  const double r_fill         = Pgetd("BH_R_size")*(1+EPS2);
  const double Ma             = Pgetd("BH_irreducible_mass");
  const double u0_Beta_U0     = 0;
  const double u0_Beta_U1     = 0;
  const double u0_Beta_U2     = 0;
  const double u0__gamma_D0D0 = 1;
  const double u0__gamma_D0D1 = 1;
  const double u0__gamma_D0D2 = 1;
  const double u0__gamma_D1D1 = 1;
  const double u0__gamma_D1D2 = 1;
  const double u0__gamma_D2D2 = 1;
  const double u0_K           = 0;
  Needle_T *patch_numbers = 0;
  const unsigned npi = 7;/* number of patches inside BH */
  const unsigned npo = 6;/* number of patches outside BH */
  unsigned p;
  
  /* check if it is perfect sphere */
  if (!Pcmps("BH_R_type","PerfectSphere"))
    Error0("This function is used when the BH surface is a perfect sphere!");
  
  /* update coeffs to avoid race condition */
  patch_numbers       = alloc_needle();
  patch_numbers->grid = grid;
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_up",grid));
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_down",grid));
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_left",grid));
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_right",grid));
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_back",grid));
  needle_in(patch_numbers,GetPatch("right_BH_surrounding_front",grid));
  assert(patch_numbers->Nin == npo);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npo; p++)
  {
    Patch_T *patch = grid->patch[patch_numbers->in[p]];
    unsigned f;
    
    bbn_preparing_conformal_metric_derivatives(patch);
    
    Field_T *R1_f  = patch->CoordSysInfo->CubedSphericalCoord->R1_f;
    Field_T *R2_f  = patch->CoordSysInfo->CubedSphericalCoord->R2_f;
    if (R1_f)
      make_coeffs_2d(R1_f,0,1);/* X and Y direction */
    if (R2_f)
      make_coeffs_2d(R2_f,0,1);/* X and Y direction */
    
    /* make coeffs for all fields inside this patch */
    for (f = 0; f < patch->nfld; ++f)
    {
      if (patch->pool[f]->v      &&
          patch->pool[f] != R1_f && 
          patch->pool[f] != R2_f    )
        make_coeffs_3d(patch->pool[f]);
    }
    
  }
  free_needle(patch_numbers);
  
  /* extrapolate */
  patch_numbers       = alloc_needle();
  patch_numbers->grid = grid;
  needle_in(patch_numbers,GetPatch("right_BH_up",grid));
  needle_in(patch_numbers,GetPatch("right_BH_down",grid));
  needle_in(patch_numbers,GetPatch("right_BH_left",grid));
  needle_in(patch_numbers,GetPatch("right_BH_right",grid));
  needle_in(patch_numbers,GetPatch("right_BH_back",grid));
  needle_in(patch_numbers,GetPatch("right_BH_front",grid));
  needle_in(patch_numbers,GetPatch("right_central_box",grid));
  assert(patch_numbers->Nin == npi);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < npi; p++)
  {
    Patch_T *patch = grid->patch[patch_numbers->in[p]];
    unsigned nn = patch->nn;
    double Y;
    double theta = 0,phi = 0,dr = 0;
    double x_on_BHsurf[3]={0},
           X_on_BHsurf[3]={0},
           N[3] = {0};
    unsigned ijk;
    
    /* fill needle */
    Needle_T *needle = alloc_needle();
    needle->grid = grid;
    needle_in(needle,GetPatch("right_BH_surrounding_up",grid));
    needle_in(needle,GetPatch("right_BH_surrounding_down",grid));
    needle_in(needle,GetPatch("right_BH_surrounding_left",grid));
    needle_in(needle,GetPatch("right_BH_surrounding_right",grid));
    needle_in(needle,GetPatch("right_BH_surrounding_back",grid));
    needle_in(needle,GetPatch("right_BH_surrounding_front",grid));
    
    bbn_add_fields_in_patch(patch);
    REALLOC_v_WRITE_v(Beta_U0)
    REALLOC_v_WRITE_v(Beta_U1)
    REALLOC_v_WRITE_v(Beta_U2)
    REALLOC_v_WRITE_v(B0_U0)
    REALLOC_v_WRITE_v(B0_U1)
    REALLOC_v_WRITE_v(B0_U2)
    REALLOC_v_WRITE_v(psi)
    REALLOC_v_WRITE_v(eta)
    REALLOC_v_WRITE_v(K)
    REALLOC_v_WRITE_v(_gamma_D2D2)
    REALLOC_v_WRITE_v(_gamma_D0D2)
    REALLOC_v_WRITE_v(_gamma_D0D0)
    REALLOC_v_WRITE_v(_gamma_D0D1)
    REALLOC_v_WRITE_v(_gamma_D1D2)
    REALLOC_v_WRITE_v(_gamma_D1D1)
    
    REALLOC_v_WRITE_v(_gammaI_U2U2)
    REALLOC_v_WRITE_v(_gammaI_U0U2)
    REALLOC_v_WRITE_v(_gammaI_U0U0)
    REALLOC_v_WRITE_v(_gammaI_U0U1)
    REALLOC_v_WRITE_v(_gammaI_U1U2)
    REALLOC_v_WRITE_v(_gammaI_U1U1)
    
    /* making B1 */
    bbn_update_B1_U012(patch);
    READ_v(B1_U0)
    READ_v(B1_U1)
    READ_v(B1_U2)
    for (ijk = 0; ijk < nn; ++ijk)
    {
      DEF_RELATIVE_x
      DEF_RELATIVE_y
      DEF_RELATIVE_z
      DEF_RELATIVE_r
      Patch_T *BHsurf_patch = 0;
      if (!EQL(r,0))
      {
        N[0]  = x/r;
        N[1]  = y/r;
        N[2]  = z/r;
        theta = acos(z/r);
      }
      else
      {
        N[0]  = 0;
        N[1]  = 0;
        N[2]  = 0;
        r     = EPS;
        theta = 0;
      }
      dr    = r - r_fill;
      phi   = arctan(y,x);
      Y     = 0.5*(1+tanh(48./125.*(r_fill/(r_fill-r)-3./2.*(r_fill/r))));
      assert(isfinite(Y));
      
      x_on_BHsurf[0] = r_fill*sin(theta)*cos(phi)+patch->c[0];
      x_on_BHsurf[1] = r_fill*sin(theta)*sin(phi)+patch->c[1];
      x_on_BHsurf[2] = r_fill*cos(theta)         +patch->c[2];
      
      /* find the patch and X which has this point */
      needle->x = x_on_BHsurf;
      point_finder(needle);
      if (!needle->Nans)
        Error0("Could not find the given point!\n");
      BHsurf_patch = grid->patch[needle->ans[0]];
      assert(X_of_x(X_on_BHsurf,x_on_BHsurf,BHsurf_patch));
      _free(needle->ans);
      needle->ans  = 0;
      needle->Nans = 0;
      
      /* extrapolate */
      double u0_psi = 2+Ma/(2*r);
      double u0_eta = 0.1*u0_psi;
      
      WTGR_EXTRAPOLATE_scalar(psi)
      WTGR_EXTRAPOLATE_scalar(eta)
      WTGR_EXTRAPOLATE_scalar(K)

      WTGR_EXTRAPOLATE_Beta(Beta_U0)
      WTGR_EXTRAPOLATE_Beta(Beta_U1)
      WTGR_EXTRAPOLATE_Beta(Beta_U2)
      
      B0_U0[ijk] = Beta_U0[ijk]-B1_U0[ijk];
      B0_U1[ijk] = Beta_U1[ijk]-B1_U1[ijk];
      B0_U2[ijk] = Beta_U2[ijk]-B1_U2[ijk];
      
      WTGR_EXTRAPOLATE_gammabar(gamma_D2D2)
      WTGR_EXTRAPOLATE_gammabar(gamma_D0D2)
      WTGR_EXTRAPOLATE_gammabar(gamma_D0D0)
      WTGR_EXTRAPOLATE_gammabar(gamma_D0D1)
      WTGR_EXTRAPOLATE_gammabar(gamma_D1D2)
      WTGR_EXTRAPOLATE_gammabar(gamma_D1D1)
      
      /* _gammaI =  _gamma inverse */
      COMPUTE_gammaI(_gamma_D0D0[ijk],_gamma_D0D1[ijk],_gamma_D0D2[ijk],
                     _gamma_D0D1[ijk],_gamma_D1D1[ijk],_gamma_D1D2[ijk],
                     _gamma_D0D2[ijk],_gamma_D1D2[ijk],_gamma_D2D2[ijk])
                     
      /* quick test check _gamma * _gammaI = delta */
      if (0)
      {
          double delta_U0D0 = 
        _gammaI_U0U0[ijk]*_gamma_D0D0[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D0D2[ijk];

          double delta_U0D1 = 
        _gammaI_U0U0[ijk]*_gamma_D0D1[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D1D2[ijk];

          double delta_U0D2 = 
        _gammaI_U0U0[ijk]*_gamma_D0D2[ijk] + _gammaI_U0U1[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U0U2[ijk]*_gamma_D2D2[ijk];

          double delta_U1D2 = 
        _gammaI_U0U1[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U1U2[ijk]*_gamma_D2D2[ijk];

          double delta_U1D0 = 
        _gammaI_U0U1[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D0D2[ijk];

         double delta_U1D1 = 
        _gammaI_U0U1[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U1[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D1D2[ijk];

          double delta_U2D2 = 
        _gammaI_U0U2[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D1D2[ijk] + _gammaI_U2U2[ijk]*_gamma_D2D2[ijk];

          double delta_U2D0 = 
        _gammaI_U0U2[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D0D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D0D2[ijk];

          double delta_U2D1 = 
        _gammaI_U0U2[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*
        _gamma_D1D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D1D2[ijk];

        if(!EQL(delta_U1D1,1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D1,0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D2,0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U1D2,0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U0D0,1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D1,0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D2,1))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U2D0,0))  Error0("_gammaI is not correct!\n");
        if(!EQL(delta_U1D0,0))  Error0("_gammaI is not correct!\n");
      }
    }
    /* free */
    free_needle(needle);
    
    /* update derivatives */
    bbn_update_derivative_Beta_U0(patch);
    bbn_update_derivative_Beta_U1(patch);
    bbn_update_derivative_Beta_U2(patch);
    bbn_update_derivative_psi(patch);
    bbn_update_derivative_eta(patch);
    
    /* compute _Aij inside BH patches */
    bbn_preparing_conformal_metric_derivatives(patch);
    bbn_free_data_Gamma_patch(patch);
    bbn_free_conformal_metric_derivatives(patch);
    bbn_update_psi10A_UiUj(patch);
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  /* free */
  free_needle(patch_numbers);
  bbn_free_conformal_metric_derivatives(GetPatch("right_BH_surrounding_up",grid));
  bbn_free_conformal_metric_derivatives(GetPatch("right_BH_surrounding_down",grid));
  bbn_free_conformal_metric_derivatives(GetPatch("right_BH_surrounding_left",grid));
  bbn_free_conformal_metric_derivatives(GetPatch("right_BH_surrounding_right",grid));
  bbn_free_conformal_metric_derivatives(GetPatch("right_BH_surrounding_back",grid));
  bbn_free_conformal_metric_derivatives(GetPatch("right_BH_surrounding_front",grid));
}

#define ij(i,j) ((j)+Nphi*(i))
/* extrapolating the given field inside of BH using Ylm method.
// field(r,theta,phi) = Sum{C_lm * r^-(l+1) * Ylm(theta,phi)} 
//
// in this method we:
// first : pick a R.
// second: finding the coeffs Clm using:
//         field(R_min,theta,phi) = Sum{Clm * r^-(l+1) * Ylm(theta,phi)}
// third : we interpolate using the Ylm expansion. */
static void extrapolate_insideBH_CS_C0_Ylm(Grid_T *const grid,const char *const field_name)
{
  const double FRACTION = 1.;/* if you wanna use R_min = FRACTION*R_min */
  const unsigned lmax = (unsigned)Pgeti("NS_surface_Ylm_expansion_max_l");
  unsigned Ntheta,Nphi;/* total number of theta and phi points */
  double *field_R_min = 0;/* field(R_min,theta,phi) */
  double *realClm,*imagClm;/* Clm coeffs */
  double R_min = Pgetd("BH_R_size");
  unsigned ijk,i,j,l,m,lm,p;
  
  /* check if it is perfect sphere */
  if (!Pcmps("BH_R_type","PerfectSphere"))
    Error0("This function is used when the BH surface is a perfect sphere!");
    
  /* initialize tables */
  init_Legendre_root_function();
  
  Ntheta  = Nphi = 2*lmax+1;
  field_R_min = alloc_double(Ntheta*Nphi);
  
  R_min = FRACTION*R_min;
  
  /* populate field(R_min,theta,phi) */
  for (i = 0; i < Ntheta; ++i)
  {
    double theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      double phi = j*2*M_PI/Nphi;
      Patch_T *patch = 0;
      double X[3],x[3];
      
      /* find patch for the given theta and phi */
      find_XYZ_and_patch_of_theta_phi_BH_CS(X,&patch,theta,phi,grid);
      
      /* r = R_min(sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^) */
      x[0]  = R_min*sin(theta)*cos(phi);
      x[1]  = R_min*sin(theta)*sin(phi);
      x[2]  = R_min*cos(theta);
      x[0] += patch->c[0];
      x[1] += patch->c[1];
      x[2] += patch->c[2];
      
      assert(X_of_x(X,x,patch));
        
      /* find field at the (X,Y,Z) */
      Interpolation_T *interp = init_interpolation();
      interp->field = patch->pool[Ind(field_name)];
      interp->XYZ_dir_flag = 1;
      interp->X            = X[0];
      interp->Y            = X[1];
      interp->Z            = X[2];
      plan_interpolation(interp);
      field_R_min[ij(i,j)] = execute_interpolation(interp);
      
      free_interpolation(interp);
    }
  }/* end of for (i = 0; i < Ntheta; ++i) */
  
  /* finding Clm in the expansion:
  // field(r,theta,phi) = Sum{Clm * r^-(l+1) * Ylm(theta,phi)} 
  // note: the coeffs need to be multiplied by R_min^(l+1) */
  realClm = alloc_ClmYlm(lmax);
  imagClm = alloc_ClmYlm(lmax);
  get_Ylm_coeffs(realClm,imagClm,field_R_min,Ntheta,Nphi,lmax);
  
  /* multiplying coeff by R_min^(l+1) */
  for (l = 0; l <= lmax; ++l)
  {
    double mult = pow(R_min,l+1);
    for (m = 0; m <= l; ++m)
    {
      lm = lm2n(l,m);
      realClm[lm] *= mult;
      imagClm[lm] *= mult;
    }
  }
  
  /* having found Clm's we using Ylm interpolation
  // field(r,theta,phi) = Sum{C_lm * r^-(l+1) * Ylm(theta,phi)} 
  // to extrapolate the field outside the BH. */
  FOR_ALL_PATCHES(p,grid)
  {
    /* surrounding patch */
    Patch_T *patch = grid->patch[p];
    const unsigned nn = patch->nn;
    Field_T *field;
    double *v;
    double x[3],r,theta,phi;
    
    if (!IsItInsideBHPatch(patch))
      continue;
      
    field = patch->pool[Ind(field_name)];
    empty_field(field);
    field->v = alloc_double(patch->nn);
    v = field->v;
    
    /* interpolate */
    for (ijk = 0; ijk < nn; ++ijk)
    {
      x[0]   = patch->node[ijk]->x[0]-patch->c[0];
      x[1]   = patch->node[ijk]->x[1]-patch->c[1];
      x[2]   = patch->node[ijk]->x[2]-patch->c[2];
      r      = root_square(3,x,0);
      theta  = acos(x[2]/r);
      phi    = arctan(x[1],x[0]);
      v[ijk] = interpolate_Clm_r_Ylm_3d(realClm,imagClm,lmax,r,theta,phi);
    }
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  /* free */
  _free(field_R_min);
  _free(realClm);
  _free(imagClm);
}
#ifdef ij
#undef ij
#endif

/* extrapolating phi, dphi and W in NS surrounding coords 
// in case they are needed for interpolation to the next grid
// or in calculation of enthalpy at NS surrounding patches.
// for extrapolation we demand:
// the phi field and its normal derivative at the NS surface be
// continues and it decreases exponentially, i.e. phi = a*exp(-att*(r-r0))+b.
// W fields keep their forms and enthalpy is made using the general formula. */
static void extrapolate_outsideNS_CS_exp_continuity_method(Grid_T *const grid)
{
  const double Omega_NS_x = Pgetd("NS_Omega_U0");
  const double Omega_NS_y = Pgetd("NS_Omega_U1");
  const double Omega_NS_z = Pgetd("NS_Omega_U2");
  const double att = 0.1;/* exp(-att*(r2-r1)) */
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    /* surrounding patch */
    Patch_T *patch = grid->patch[p];
    const unsigned *n = patch->n;
    const unsigned nn = patch->nn;
    unsigned ijk,i,j,k;
    
    if (!IsItNSSurroundingPatch(patch))
      continue;
     
    /* irrotational part of fluid */
    REALLOC_v_WRITE_v(phi)
    DECLARE_AND_EMPTY_FIELD(dphi_D2)
    DECLARE_AND_EMPTY_FIELD(dphi_D1)
    DECLARE_AND_EMPTY_FIELD(dphi_D0)
    
    /* spin part of fluid W^i */
    REALLOC_v_WRITE_v(W_U0)
    REALLOC_v_WRITE_v(W_U1)
    REALLOC_v_WRITE_v(W_U2)
    
    /* populate the spin part */
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x = patch->node[ijk]->x[0]-patch->c[0];
      double y = patch->node[ijk]->x[1]-patch->c[1];
      double z = patch->node[ijk]->x[2]-patch->c[2];
      
      /* spin part */
      W_U0[ijk] = Omega_NS_y*z-Omega_NS_z*y;
      W_U1[ijk] = Omega_NS_z*x-Omega_NS_x*z;
      W_U2[ijk] = Omega_NS_x*y-Omega_NS_y*x;
    }

    /* populat the irrotational part */
    Patch_T *NS_patch;/* corresponding NS patch, to extrapolate out */
    double r0,r,a,b;
    double phii,d0phii,d1phii,d2phii;/* fields at r0 */
    double x[3],X[3],N[3],Norm;
    
    /* find the corresponding NS patch to be used for extrapolation */
    char stem[1000];
    char *affix = regex_find("_[[:alpha:]]{2,5}$",patch->name);/* finding the side of the patch */
    assert(affix);
    sprintf(stem,"left_NS%s",affix);
    free(affix);
    NS_patch = GetPatch(stem,grid);
    
    /* prepare interpolation arguments for phi and dphi */
    Interpolation_T *interp_phi = init_interpolation();
    interp_phi->field = NS_patch->pool[LookUpField_E("phi",NS_patch)];
    interp_phi->XYZ_dir_flag = 1;
    plan_interpolation(interp_phi);
    
    Interpolation_T *interp_d0phi = init_interpolation();
    Interpolation_T *interp_d1phi = init_interpolation();
    Interpolation_T *interp_d2phi = init_interpolation();
    interp_d0phi->field = NS_patch->pool[LookUpField_E("dphi_D0",NS_patch)];
    interp_d1phi->field = NS_patch->pool[LookUpField_E("dphi_D1",NS_patch)];
    interp_d2phi->field = NS_patch->pool[LookUpField_E("dphi_D2",NS_patch)];
    interp_d0phi->XYZ_dir_flag = 1;
    interp_d1phi->XYZ_dir_flag = 1;
    interp_d2phi->XYZ_dir_flag = 1;
    
    plan_interpolation(interp_d0phi);
    plan_interpolation(interp_d1phi);
    plan_interpolation(interp_d2phi);

    /* spreading out the fields value using interpolation */
    for (i = 0; i < n[0]; ++i)
    {
      for (j = 0; j < n[1]; ++j)
      {
        /* calculate r0 at NS surface */
        ijk = L(n,i,j,0);
        x[0] = patch->node[ijk]->x[0]-patch->c[0];
        x[1] = patch->node[ijk]->x[1]-patch->c[1];
        x[2] = patch->node[ijk]->x[2]-patch->c[2];
        r0   = root_square(3,x,0);

        /* unit normal vector on NS surface pointing outward */
        N[0] = dq2_dq1(patch,_c_,_x_,ijk);
        N[1] = dq2_dq1(patch,_c_,_y_,ijk);
        N[2] = dq2_dq1(patch,_c_,_z_,ijk);        
        Norm = root_square(3,N,0);
        N[0] /= Norm;
        N[1] /= Norm;
        N[2] /= Norm;
        
        X_of_x(X,patch->node[ijk]->x,NS_patch);

        interp_phi->X    = X[0];
        interp_phi->Y    = X[1];
        interp_phi->Z    = X[2];

        interp_d0phi->X  = X[0];
        interp_d0phi->Y  = X[1];
        interp_d0phi->Z  = X[2];
        
        interp_d1phi->X  = X[0];
        interp_d1phi->Y  = X[1];
        interp_d1phi->Z  = X[2];
        
        interp_d2phi->X  = X[0];
        interp_d2phi->Y  = X[1];
        interp_d2phi->Z  = X[2];

        phii    = execute_interpolation(interp_phi);
        d0phii  = execute_interpolation(interp_d0phi);
        d1phii  = execute_interpolation(interp_d1phi);
        d2phii  = execute_interpolation(interp_d2phi);
        
        a = -(N[0]*d0phii+N[1]*d1phii+N[2]*d2phii)/att;
        b = phii-a;
        
        for (k = 0; k < n[2]; ++k)
        {
          ijk = L(n,i,j,k);
          x[0] = patch->node[ijk]->x[0]-patch->c[0];
          x[1] = patch->node[ijk]->x[1]-patch->c[1];
          x[2] = patch->node[ijk]->x[2]-patch->c[2];
          
          r = root_square(3,x,0);
          
          phi[ijk]  = a*exp(-att*(r-r0))+b;
          
        }/* end of for (k = 0 ; k < n[2]; ++k) */
      }/* end of for (j = 0; j < n[1]; ++j) */
    }/* end of for (i = 0; i < n[0]; ++i) */
    free_interpolation(interp_phi);
    free_interpolation(interp_d0phi);
    free_interpolation(interp_d1phi);
    free_interpolation(interp_d2phi);

    Field_T *phi_field = patch->pool[Ind("phi")];
    dphi_D2->v = Partial_Derivative(phi_field,"z");
    dphi_D1->v = Partial_Derivative(phi_field,"y");
    dphi_D0->v = Partial_Derivative(phi_field,"x");
    
    /*
    Tij_IF_CTS_enthalpy(patch);
      
    Field_T *enthalpy = patch->pool[Ind("enthalpy")];
    DECLARE_AND_EMPTY_FIELD(denthalpy_D2)
    DECLARE_AND_EMPTY_FIELD(denthalpy_D1)
    DECLARE_AND_EMPTY_FIELD(denthalpy_D0)
    denthalpy_D2->v = Partial_Derivative(enthalpy,"z");
    denthalpy_D1->v = Partial_Derivative(enthalpy,"y");
    denthalpy_D0->v = Partial_Derivative(enthalpy,"x");
    */
    
  }/* end of FOR_ALL_PATCHES(p,grid) */
}

#define ij(i,j) ((j)+Nphi*(i))
/* extrapolating the given field outside of NS using Ylm method.
// field(r,theta,phi) = Sum{C_lm * r^-(l+1) * Ylm(theta,phi)} 
//
// in this method we:
// first : find the minimum radius of NS R_min.
// second: finding the coeffs Clm using:
//         field(R_min,theta,phi) = Sum{Clm * r^-(l+1) * Ylm(theta,phi)}
// third : we interpolate using the Ylm expansion. */
static void extrapolate_outsideNS_CS_Ylm_method(Grid_T *const grid,const char *const field_name)
{
  const double FRACTION = 1.;/* if you wanna use R_min = FRACTION*R_min */
  const unsigned lmax = (unsigned)Pgeti("NS_surface_Ylm_expansion_max_l");
  unsigned Ntheta,Nphi;/* total number of theta and phi points */
  double *field_R_min = 0;/* field(R_min,theta,phi) */
  double *realClm,*imagClm;/* Clm coeffs */
  double R_min = DBL_MAX;/* min radius of NS */
  unsigned ijk,i,j,l,m,lm,p;
  
  /* initialize tables */
  init_Legendre_root_function();
  
  Ntheta  = Nphi = 2*lmax+1;
  field_R_min = alloc_double(Ntheta*Nphi);
  
  /* find the minimum radius of on the NS surface */
  for (i = 0; i < Ntheta; ++i)
  {
    double theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      double phi = j*2*M_PI/Nphi;
      Patch_T *patch = 0;
      double X[3],x[3],r;
      
      /* find patch and X,Y,Z at NS surface in which theta and phi take place */
      find_XYZ_and_patch_of_theta_phi_NS_CS(X,&patch,theta,phi,grid);
      
      /* finding x */
      assert(x_of_X(x,X,patch));
      x[0] -= patch->c[0];
      x[1] -= patch->c[1];
      x[2] -= patch->c[2];  
      r     = root_square(3,x,0);
      if (r < R_min)
      {
        R_min = r;
      }
    }
  }/* end of for (i = 0; i < Ntheta; ++i) */
  
  R_min = FRACTION*R_min;
  
  /* populate field(R_min,theta,phi) */
  for (i = 0; i < Ntheta; ++i)
  {
    double theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      double phi = j*2*M_PI/Nphi;
      Patch_T *patch = 0;
      double X[3],x[3];
      
      /* find patch for the given theta and phi */
      find_XYZ_and_patch_of_theta_phi_NS_CS(X,&patch,theta,phi,grid);
      
      /* r = R_min(sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^) */
      x[0]  = R_min*sin(theta)*cos(phi);
      x[1]  = R_min*sin(theta)*sin(phi);
      x[2]  = R_min*cos(theta);
      x[0] += patch->c[0];
      x[1] += patch->c[1];
      x[2] += patch->c[2];
      
      assert(X_of_x(X,x,patch));
        
      /* find field at the (X,Y,Z) */
      Interpolation_T *interp = init_interpolation();
      interp->field = patch->pool[Ind(field_name)];
      interp->XYZ_dir_flag = 1;
      interp->X            = X[0];
      interp->Y            = X[1];
      interp->Z            = X[2];
      plan_interpolation(interp);
      field_R_min[ij(i,j)] = execute_interpolation(interp);
      
      free_interpolation(interp);
    }
  }/* end of for (i = 0; i < Ntheta; ++i) */
  
  /* finding Clm in the expansion:
  // field(r,theta,phi) = Sum{Clm * r^-(l+1) * Ylm(theta,phi)} 
  // note: the coeffs need to be multiplied by R_min^(l+1) */
  realClm = alloc_ClmYlm(lmax);
  imagClm = alloc_ClmYlm(lmax);
  get_Ylm_coeffs(realClm,imagClm,field_R_min,Ntheta,Nphi,lmax);
  
  /* multiplying coeff by R_min^(l+1) */
  for (l = 0; l <= lmax; ++l)
  {
    double mult = pow(R_min,l+1);
    for (m = 0; m <= l; ++m)
    {
      lm = lm2n(l,m);
      realClm[lm] *= mult;
      imagClm[lm] *= mult;
    }
  }
  
  /* having found Clm's we using Ylm interpolation
  // field(r,theta,phi) = Sum{C_lm * r^-(l+1) * Ylm(theta,phi)} 
  // to extrapolate the field outside the NS. */
  FOR_ALL_PATCHES(p,grid)
  {
    /* surrounding patch */
    Patch_T *patch = grid->patch[p];
    const unsigned nn = patch->nn;
    Field_T *field;
    double *v;
    double x[3],r,theta,phi;
    
    if (!IsItNSSurroundingPatch(patch))
      continue;
      
    field = patch->pool[Ind(field_name)];
    empty_field(field);
    field->v = alloc_double(patch->nn);
    v = field->v;
    
    /* interpolate */
    for (ijk = 0; ijk < nn; ++ijk)
    {
      x[0]   = patch->node[ijk]->x[0]-patch->c[0];
      x[1]   = patch->node[ijk]->x[1]-patch->c[1];
      x[2]   = patch->node[ijk]->x[2]-patch->c[2];
      r      = root_square(3,x,0);
      theta  = acos(x[2]/r);
      phi    = arctan(x[1],x[0]);
      v[ijk] = interpolate_Clm_r_Ylm_3d(realClm,imagClm,lmax,r,theta,phi);
    }
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  /* free */
  _free(field_R_min);
  _free(realClm);
  _free(imagClm);
}
#ifdef ij
#undef ij
#endif

/* given r, theta, phi, it interpolates using:
// field(r,theta,phi) = Sum{Clm * r^-(l+1) * Ylm(theta,phi)} */
static double interpolate_Clm_r_Ylm_3d(double *const realClm,double *const imagClm,const unsigned lmax,const double r,const double theta,const double phi)
{
  double interp = 0;
  unsigned l,m,lm;
  double rpow;
  
  /* multiplying coeff by r^-(l+1) */
  for (l = 0; l <= lmax; ++l)
  {
    rpow = pow(r,l+1);
    for (m = 0; m <= l; ++m)
    {
      lm = lm2n(l,m);
      realClm[lm] /= rpow;
      imagClm[lm] /= rpow;
    }
  }
  
  interp = interpolation_Ylm(realClm,imagClm,lmax,theta,phi);
  
  /* retuning back the coeffs => multiplying coeff by r^(l+1) */
  for (l = 0; l <= lmax; ++l)
  {
    rpow = pow(r,l+1);
    for (m = 0; m <= l; ++m)
    {
      lm = lm2n(l,m);
      realClm[lm] *= rpow;
      imagClm[lm] *= rpow;
    }
  }
  
  return interp;
}

/* extrapolating phi, dphi and W in NS surrounding coords 
// in case they are needed for interpolation to the next grid
// or in calculation of enthalpy at NS surrounding patches.
// for extrapolation we demand:
// the fields spread out in the same fashion as it is changing toward 
// the NS sarface. what we have :
// f(r_out) = (f(r_in) + df)*exp(g(r)), in which df is f(r2)-f(r1), 
// r1 = FACTOR*r2 and r2 is the radius of NS surface, and g(r) is 
// a function of r_out/r2 to control the radial trend of field.
// a = (r_max-r2)/(r2-r1)
// b = (r1*r_max-Pow2(r2))/(r1-r2)
// r_out = a*r_in + b, where (r_in) r_out is r (inside)outside NS and
// r_max is the max radius in NS surrounding patch.
// in effect, it means we emulate the trend of the fields inside of NS
// from r1 to r2 in NS surrounding, from r2 to r_max. */
static void extrapolate_outsideNS_CS_slop_method(Grid_T *const grid)
{
  const double FACTOR = 0.9;/* r1 = FACTOR*r2 */
  const double att    = 0.1;/* exp(-att*(r2-r1)) */
  //const double EXP    = 1;/* (r_out/r2)^EXP */
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    /* surrounding patch */
    Patch_T *patch = grid->patch[p];
    if (!IsItNSSurroundingPatch(patch))
      continue;
     
    /* irrotational part of fluid */
    REALLOC_v_WRITE_v(phi)
    DECLARE_AND_EMPTY_FIELD(dphi_D2)
    DECLARE_AND_EMPTY_FIELD(dphi_D1)
    DECLARE_AND_EMPTY_FIELD(dphi_D0)
    
    /* spin part of fluid W^i */
    REALLOC_v_WRITE_v(W_U0)
    REALLOC_v_WRITE_v(W_U1)
    REALLOC_v_WRITE_v(W_U2)
    
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
        r_max = root_square(3,x,0);

        /* calculate r2 at NS surface */
        ijk = L(n,i,j,0);
        x[0] = patch->node[ijk]->x[0]-patch->c[0];
        x[1] = patch->node[ijk]->x[1]-patch->c[1];
        x[2] = patch->node[ijk]->x[2]-patch->c[2];
        r2   = root_square(3,x,0);

        /* calulate r1,a,b */
        r1 = FACTOR*r2;/* small scale to define r1 */
        a  = (r_max-r2)/(r2-r1);
        b  = (r1*r_max-Pow2(r2))/(r1-r2);

        /* find the value of phi and W^{i} at r2 */
        THETA = acos(x[2]/r2);
        PHI   = arctan(x[1],x[0]);
        
        /* check of r1 is inside NS */
        double x1[3];
        x1[0] = r1*sin(THETA)*cos(PHI)+NS_patch->c[0];
        x1[1] = r1*sin(THETA)*sin(PHI)+NS_patch->c[1];
        x1[2] = r1*cos(THETA)+NS_patch->c[2];
        Needle_T *needle = alloc_needle();
        needle->grid = grid;
        needle->x    = x1;
        needle_in(needle,NS_patch);
        point_finder(needle);
        if (!needle->Nans)
        {
          Error0("The FACTOR is too larg for extrapolation.\n");
        }
        free_needle(needle);
        /* end of check */
        
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
          
          r_out = root_square(3,x,0);
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
    if (1)/* decrease enthalpy exponentially from the NS surface */
    {
      /* note: patch refers to NS_surrounding patch */
      
      REALLOC_v_WRITE_v(enthalpy)
      /* find the corresponding NS patch to be used for extrapolation */
      affix = regex_find("_[[:alpha:]]{2,5}$",patch->name);/* finding the side of the patch */
      assert(affix);
      sprintf(stem,"left_NS%s",affix);
      free(affix);
      NS_patch = GetPatch(stem,grid);
      
      /* prepare interpolation arguments */
      Interpolation_T *interp_h = init_interpolation();
      interp_h->field = NS_patch->pool[LookUpField_E("enthalpy",NS_patch)];
      interp_h->XYZ_dir_flag = 1;
      plan_interpolation(interp_h);

      /* smearing the fields exponentially */
      for (i = 0; i < n[0]; ++i)
      {
        for (j = 0; j < n[1]; ++j)
        {
          /* calculate r at the NS surface */
          ijk = L(n,i,j,0);
          x[0] = patch->node[ijk]->x[0]-patch->c[0];
          x[1] = patch->node[ijk]->x[1]-patch->c[1];
          x[2] = patch->node[ijk]->x[2]-patch->c[2]; 
          r1   = root_square(3,x,0);
          
          X_of_x(X,patch->node[ijk]->x,NS_patch);

          interp_h->X  = X[0];
          interp_h->Y  = X[1];
          interp_h->Z  = X[2];

          double h_i   = execute_interpolation(interp_h);

          for (k = 0; k < n[2]; ++k)
          {
            ijk = L(n,i,j,k);
            x[0] = patch->node[ijk]->x[0]-patch->c[0];
            x[1] = patch->node[ijk]->x[1]-patch->c[1];
            x[2] = patch->node[ijk]->x[2]-patch->c[2];
            
            r2 = root_square(3,x,0);
            double e = exp(-att*(r2-r1));
            enthalpy[ijk] = h_i*e;
            
          }/* end of for (k = 0 ; k < n[2]; ++k) */
        }/* end of for (j = 0; j < n[1]; ++j) */
      }/* end of for (i = 0; i < n[0]; ++i) */
      free_interpolation(interp_h);
    }
    else if (0)/* use metric and matter fields to make enthalpy */
      Tij_IF_CTS_enthalpy(patch);
    else
      Error0(NO_OPTION);
      
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
  pr_line_custom('=');
  printf("{ Initializing fields & grid using TOV & Kerr-Schild solutions ...\n");
  
  Grid_T *grid = 0;
  struct Grid_Params_S *GridParams = init_GridParams();/* adjust some pars for construction of grid */
  GridParams->grid_prev = 0;
   
  /* solve fields for a TOV star located at left side of y axis */
  TOV_T *tov = TOV_init();
  tov->bar_m = Pgetd("NS_baryonic_mass");
  tov->description = "estimating NS";
  tov = TOV_solution(tov);
  const double ns_R = tov->rbar[tov->N-1];
  
  /* basics of Kerr Schild black hole located at right side of y axis */
  pr_line_custom('=');
  printf("{ Acquiring Black Hole properties ...\n");
  const double bh_chi_x    = Pgetd("BH_chi_U0");
  const double bh_chi_y    = Pgetd("BH_chi_U1");
  const double bh_chi_z    = Pgetd("BH_chi_U2");
  const double bh_irr_mass = Pgetd("BH_irreducible_mass");
  const double bh_R        = 2*bh_irr_mass;/* approximate initial radius */
  const double bh_chi      = sqrt(Pow2(bh_chi_x)+Pow2(bh_chi_y)+Pow2(bh_chi_z));
  const double bh_a        = bh_chi*bh_irr_mass;
  /* check size of bh_chi */
  if (GRT(bh_chi,1))
    Error0("BH spin is too large!\n");
  printf("BH properties:\n");
  printf("--> BH radius (Kerr-Schild Coords.) ~ %+e\n",bh_R);
  printf("--> BH irreducible mass             ~ %+e\n",bh_irr_mass);
  printf("--> BH dimensionless spin (x comp.) = %+e\n",bh_chi_x);
  printf("--> BH dimensionless spin (y comp.) = %+e\n",bh_chi_y);
  printf("--> BH dimensionless spin (z comp.) = %+e\n",bh_chi_z);
  printf("--> BH approximate net spin         ~ %+e\n",bh_a);
  printf("} Acquiring Black Hole properties ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
 
  /* center of rotation (approx. Center of Mass) */
  const double D = Pgetd("BH_NS_separation");
  const double C_BH = 0.5*D;/* center of BH patch, it's on +y axis */
  const double C_NS = -C_BH;/* center of NS it's on -y axis*/
  const double ns_mass = tov->ADM_m;/* NS adm mass */
  const double y_CM = (ns_mass*C_NS+bh_irr_mass*C_BH)/(ns_mass+bh_irr_mass);
  
  /* adding some parameters: */
  /* -> pars for adjusting NS surface and interpolation for the next grid */
  Pseti("did_resolution_change?",1);
  Pseti("did_NS_surface_change?",1);
  Pseti("did_AH_surface_change?",1);
  Pseti("did_NS_surface_finder_work?",1);/* if surface finder was failed 0 */
  Pseti("use_previous_data",0);
  
  /* -> center of mass */
  Psetd("x_CM",0);
  Psetd("y_CM",y_CM);
  Psetd("z_CM",0);
  //Psetd("y_CM0",y_CM);
  //Psetd("x_CM0",0);
  
  /* -> NS properties */
  Psetd("NS_center_x",0);
  Psetd("NS_center_y",C_NS);
  Psetd("NS_center_z",0);
  Psetd("NS_ADM_mass",1);
  
  /* -> BH properties */  
  Psetd("r_excision",bh_R);
  Psetd("BH_net_spin",bh_a);
  Psetd("BH_center_x",0);
  Psetd("BH_center_y",C_BH);
  Psetd("BH_center_z",0);
  Psetd("BH_Vz",0);/* BH velocity in z direction */
  
  /* -> approximate BH_Omega, the angular frequency of the horizon,
  // is a free vector that determines the spin of BH
  // and it is related to the dimensionless spin by:
  // BH_chi = 4*BH_mass*BH_Omega (PRD 86 084033). */
  Psetd("BH_Omega_U0",bh_chi_x/(4*bh_irr_mass));
  Psetd("BH_Omega_U1",bh_chi_y/(4*bh_irr_mass));
  Psetd("BH_Omega_U2",bh_chi_z/(4*bh_irr_mass));

  /* -> the Constant of the integration of Euler equation */
  Psetd("Euler_equation_constant",0);
  
  /* -> central rho0 */
  Psetd("rho_center",1E-3);
  
  /* -> errors: */
  Psetd("largest_L2norm_error",0);/* max error of calculated L2norms, 
                                  // residual or constraint violation and etc. */
  
  /* combining these two geometry to create the grid */
  GridParams->Max_R_NS_l = ns_R;
  GridParams->R_BH_r     = bh_R;
  GridParams->a_BH       = DBL_MAX;/* to catch error! *///bh_chi*bh_mass;
  GridParams->NS_R_type  = "PerfectSphere";
  GridParams->BH_R_type  = "PerfectSphere";
  grid = creat_bbn_grid_CS(GridParams);
  
  /* adding all of the fields needed for construction of Initial Data */
  bbn_add_fields(grid);
  
  /* populating the free data part of initial data that we chose ourself */
  bbn_populate_free_data(grid);
  
  /* initialize the fields using TOV and Kerr-Schild solution. */
  init_field_TOV_plus_KerrSchild(grid,tov);
  
  /* make normal vectorn on BH horizon 
  // note: this MUST be before "bbn_partial_derivatives_fields" */
  bbn_make_normal_vector_on_BH_horizon(grid);
  
  /* taking partial derivatives of the fields needed for equations */
  bbn_partial_derivatives_fields(grid);
  
  /* update u0, _J^i, _E and _S */
  Tij_IF_CTS_psi6Sources(grid);
  
  /* update u0 derivatives */
  unsigned p;
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    bbn_update_derivative_u0(patch);
  }
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  bbn_update_Aij(grid);
  
  /* find Euler equation const using enthalpy of TOV star and other fields */
  find_Euler_eq_const_TOV_KerrSchild(grid);
  
  /* add some parameters for momentum and its adjustments */
  Observable_T *obs = init_observable(grid,bbn_plan_obs_CS,bbn_free_obs_CS);
  obs->quantity     = "ADM(P,J)|BBN";
  plan_observable(obs);
  Psetd("Px_ADM",obs->Px(obs));
  Psetd("Py_ADM",obs->Py(obs));
  Psetd("Pz_ADM",obs->Pz(obs));
  Psetd("Jx_ADM",obs->Jx(obs));
  Psetd("Jy_ADM",obs->Jy(obs));
  Psetd("Jz_ADM",obs->Jz(obs));
  Psetd("v*_boost_x",0);
  Psetd("v*_boost_y",0);
  Psetd("v*_boost_z",0);
  
  /* freeing */
  free_Grid_Params_S(GridParams);
  TOV_free(tov);
  free_observable(obs);
  
  printf("} Initializing fields & grid using TOV & Kerr-Schild solutions ==> Done.\n");
  pr_line_custom('=');
  
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
  double Euler_C = 0;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if(!IsItNSPatch(patch))
      continue;
    
    ijk = 0;/* for an arbitrary point */
    READ_v(enthalpy)
    READ_v(_gammaI_U0U2)
    READ_v(_gammaI_U0U0)
    READ_v(_gammaI_U0U1)
    READ_v(_gammaI_U1U2)
    READ_v(_gammaI_U1U1)
    READ_v(_gammaI_U2U2)
    READ_v(W_U1)
    READ_v(W_U0)
    READ_v(W_U2)
    READ_v(dphi_D2)
    READ_v(dphi_D1)
    READ_v(dphi_D0)
    READ_v(Beta_U1)
    READ_v(Beta_U0)
    READ_v(Beta_U2)
    READ_v(psi)
    READ_v(u0)

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
  Psetd("Euler_equation_constant",Euler_C);
}

/* make normal vectorn on BH horizon */
void bbn_make_normal_vector_on_BH_horizon(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Making normal vector on BH horizon ...\n");
  
  const double BH_center_x = Pgetd("BH_center_x");
  const double BH_center_y = Pgetd("BH_center_y");
  const double BH_center_z = Pgetd("BH_center_z");
  unsigned p,nn,ijk;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItHorizonPatch(patch))
      continue;
      
    nn = patch->nn;
    
    READ_v(_gamma_D2D2)
    READ_v(_gamma_D0D2)
    READ_v(_gamma_D0D0)
    READ_v(_gamma_D0D1)
    READ_v(_gamma_D1D2)
    READ_v(_gamma_D1D1)
    
    /* normal vector on horizon */
    REALLOC_v_WRITE_v(_HS_U0);
    REALLOC_v_WRITE_v(_HS_U1);
    REALLOC_v_WRITE_v(_HS_U2);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x = patch->node[ijk]->x[0]-BH_center_x;
      double y = patch->node[ijk]->x[1]-BH_center_y; 
      double z = patch->node[ijk]->x[2]-BH_center_z;
      double r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
      
      /* minus sign to point outside the black hole */
      //_HS_U0[ijk] = dq2_dq1(patch,_c_,_x_,ijk);
      //_HS_U1[ijk] = dq2_dq1(patch,_c_,_y_,ijk);
      //_HS_U2[ijk] = dq2_dq1(patch,_c_,_z_,ijk);
      
      /* the jacobian method (above) has discontinuity at the plane
      // away from the horizon due to cubed spherical setup. */
      _HS_U0[ijk] = x/r;
      _HS_U1[ijk] = y/r;
      _HS_U2[ijk] = z/r;
      
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
     
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  printf("} Making normal vector on BH horizon ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* initialize the fields using TOV and Kerr-Schild solution.
// the idea is to superimpose two fields of each solutions. */
static void init_field_TOV_plus_KerrSchild(Grid_T *const grid,const TOV_T *const tov)
{
  pr_line_custom('=');
  printf("{ Superimposing TOV and Kerr-Schild solution ...\n");
  
  const double BH_center_x = Pgetd("BH_center_x");
  const double BH_center_y = Pgetd("BH_center_y");
  const double BH_center_z = Pgetd("BH_center_z");
  const double Omega_BHNS  = Pgetd("BH_NS_angular_velocity");
  const double Omega_NS_x  = Pgetd("NS_Omega_U0");
  const double Omega_NS_y  = Pgetd("NS_Omega_U1");
  const double Omega_NS_z  = Pgetd("NS_Omega_U2");
  const double R_Schwar    = tov->r[tov->N-1];/* NS's Schwarzchild radius */
  const double M_NS  = tov->ADM_m;/* NS adm mass */
  const double C_NSx = Pgetd("NS_center_x");/* center of NS on x axis */
  const double C_NSy = Pgetd("NS_center_y");/* center of NS on y axis */
  const double C_NSz = Pgetd("NS_center_z");/* center of NS on z axis */
  const double y_CM  = Pgetd("y_CM");
  const double a_BH  = Pgetd("BH_net_spin");
  const double M_BH  = Pgetd("BH_irreducible_mass");
  unsigned p;
  
  /* populate tB tR */
  Transformation_T *tB = initialize_transformation();
  Transformation_T *tR = initialize_transformation();
  bbn_transform_populate_boost_rotation(tB,tR);
  
  /* black hole parts */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    REALLOC_v_WRITE_v(Beta_U0)
    REALLOC_v_WRITE_v(Beta_U1)
    REALLOC_v_WRITE_v(Beta_U2)
    
    READ_v(_gammaI_U0U2)
    READ_v(_gammaI_U0U0)
    READ_v(_gammaI_U0U1)
    READ_v(_gammaI_U1U2)
    READ_v(_gammaI_U1U1)
    READ_v(_gammaI_U2U2)

    ADD_FIELD(KSbeta_D0)
    ADD_FIELD(KSbeta_D1)
    ADD_FIELD(KSbeta_D2)
    REALLOC_v_WRITE_v(KSbeta_D0)
    REALLOC_v_WRITE_v(KSbeta_D1)
    REALLOC_v_WRITE_v(KSbeta_D2)
    
    ADD_FIELD(KSalpha)
    REALLOC_v_WRITE_v(KSalpha)
    
    /* beta and alpha needed */
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x,y,z,H,k0,k1,k2,kt;
      x = patch->node[ijk]->x[0]-BH_center_x;
      y = patch->node[ijk]->x[1]-BH_center_y;
      z = patch->node[ijk]->x[2]-BH_center_z;
      
      bbn_transform_get_k_and_H_KerrSchild
        (x,y,z,a_BH,M_BH,tB,tR,&kt,&k0,&k1,&k2,&H);
      
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
  free_transformation(tB);
  free_transformation(tR);  
  
  /* initialization psi, eta and matter fields: */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    REALLOC_v_WRITE_v(psi)
    REALLOC_v_WRITE_v(eta)
    READ_v(KSalpha)
    
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
      
      REALLOC_v_WRITE_v(enthalpy)
      REALLOC_v_WRITE_v(rho0)
      REALLOC_v_WRITE_v(phi)
      REALLOC_v_WRITE_v(W_U0)
      REALLOC_v_WRITE_v(W_U1)
      REALLOC_v_WRITE_v(W_U2)
      
      for (ijk = 0; ijk < nn; ++ijk)
      {
        /* note that we naturally using isotropic coords. 
        // for our coordiante, so bar in rbar is dropped */
        double x = patch->node[ijk]->x[0]-C_NSx;
        double y = patch->node[ijk]->x[1]-C_NSy;
        double z = patch->node[ijk]->x[2]-C_NSz;
        double r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
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
        phi[ijk] = -Omega_BHNS*(C_NSy-y_CM)*x;
        
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
        double x    = patch->node[ijk]->x[0]-C_NSx;
        double y    = patch->node[ijk]->x[1]-C_NSy;
        double z    = patch->node[ijk]->x[2]-C_NSz;
        double r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
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
      
     REALLOC_v_WRITE_v(B0_U0)
     REALLOC_v_WRITE_v(B0_U1)
     REALLOC_v_WRITE_v(B0_U2)
     
     WRITE_v(Beta_U0)
     WRITE_v(Beta_U1)
     WRITE_v(Beta_U2)
     
     READ_v(B1_U0)
     READ_v(B1_U1)
     READ_v(B1_U2)
     READ_v(psi)
     
     /* for outermost patches the better approximation is B0 = 0 */
     if (IsItOutermostPatch(patch))
     {
       for (ijk = 0; ijk < nn; ++ijk)
       {
         double x = patch->node[ijk]->x[0];
         double y = patch->node[ijk]->x[1];
         double z = patch->node[ijk]->x[2];
         double r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
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
  
  printf("} Superimposing TOV and Kerr-Schild solution ==> Done.\n");
  pr_clock();
  pr_line_custom('=');

}

/* given the max of radius of NS, radius of BH and their separation,
// and BH spin as the parameters, create a grid with these properties.
// NOTE: WE assume we are using cubed spherical grid.
// ->return value: grid of NS and BH in which inside of the BH excised. */
static Grid_T *creat_bbn_grid_CS(struct Grid_Params_S *const GridParams)
{
  pr_line_custom('=');
  printf("{ Creating grid ...\n");
  
  Grid_T *grid_next = alloc_grid();/* adding a new grid */
  Grid_T *grid_prev = GridParams->grid_prev;
  
  /* calculate the characteristics of this grid */
  const double Max_R_NS_l = GridParams->Max_R_NS_l;/* maximum radius of NS */
  const double R_BH_r     = GridParams->R_BH_r;
  const unsigned gn = grid_next->gn;
  const double S    = Pgetd("BH_NS_separation");
  const unsigned N_Outermost_Split = (unsigned)Pgeti("Number_of_Outermost_Split"); 
  double *R_outermost = alloc_double(N_Outermost_Split);
  double box_size_l,box_size_r;
  unsigned nlb[3]/*left box*/,n;
  char var[100] = {'\0'};
  char par[100] = {'\0'};
  char val[100] = {'\0'};
  const char *kind;
  unsigned i,p;
  
  /* finding the kind of grid */
  kind = Pgets("grid_kind");
  
  /* some checks */
  if (!strcmp_i(kind,"BBN_CubedSpherical_grid"))
    Error0("This function only works with cubed spherical grid.\n");
  if(!GRT(S,0))
    Error0("The distance between the two compact objects must be positive.\n");
  if(!GRT(Max_R_NS_l,0))
    Error0("Neutron star must have positive radius.\n");
  if(!GRT(R_BH_r,0))
    Error0("Black hole must have positive radius.\n");
  if(!LSS(2*Max_R_NS_l,S))
    Error0("The neutron star radius is too big.\n");
  if(!LSS(2*R_BH_r,S))
    Error0("The black hole radius is too big.\n");
  
  grid_next->kind = dup_s(kind);  
  /* making NS and BH surfaces function */
  NS_BH_surface_CubedSpherical_grid(grid_next,GridParams);
  
  box_size_l = Pgetd("left_central_box_length_ratio")*Max_R_NS_l;
  box_size_r = Pgetd("left_central_box_length_ratio")*R_BH_r;/* use same ratio as NS */
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R_outermost[i] = Pgetd(var);
    
    if (LSS(R_outermost[i],2*S))
      Error0("The radius of outermost patches must be greater than twice of BBN distance.");
    
    if (i > 0)
      if (LSSEQL(R_outermost[i],R_outermost[i-1]))
        Error0("The radius of outermost must be increasing.");
    
  }
  
  /* filling n */
  
  /* left box */
  nlb[0] = (unsigned)PgetiEZ("n_a");
  nlb[1] = (unsigned)PgetiEZ("n_b");
  nlb[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("left_NS_n_a");
  if (n != INT_MAX)   nlb[0] = n;
  n = (unsigned)PgetiEZ("left_NS_n_b");
  if (n != INT_MAX)   nlb[1] = n;
  n = (unsigned)PgetiEZ("left_NS_n_c");
  if (n != INT_MAX)   nlb[2] = n;
    
  if(nlb[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(nlb[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(nlb[2] == INT_MAX)
    Error0("n_c could not be set.\n");
  
  /* adding the results to the parameter data base */
  
  /* n_a, n_b, n_c */
  /* left box */
  sprintf(par,"grid%u_left_central_box_n_a",gn);
  sprintf(val,"%u",nlb[0]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_left_central_box_n_b",gn);
  sprintf(val,"%u",nlb[1]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_left_central_box_n_c",gn);
  sprintf(val,"%u",nlb[2]);
  add_parameter_string(par,val);
  
  /* size a,b,c */
  sprintf(par,"grid%u_left_central_box_size_a",gn);
  Psetd(par,box_size_l);
  
  sprintf(par,"grid%u_left_central_box_size_b",gn);
  Psetd(par,box_size_l);
  
  sprintf(par,"grid%u_left_central_box_size_c",gn);
  Psetd(par,box_size_l);
  
  /* surrounding box length */
  sprintf(par,"grid%u_surrounding_box_length",gn);
  Psetd(par,S);
  
  /* right box. NOTE: this is needed when we fill the excision region */
  nlb[0] = (unsigned)PgetiEZ("n_a");
  nlb[1] = (unsigned)PgetiEZ("n_b");
  nlb[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("right_BH_n_a");
  if (n != INT_MAX)   nlb[0] = n;
  n = (unsigned)PgetiEZ("right_BH_n_b");
  if (n != INT_MAX)   nlb[1] = n;
  n = (unsigned)PgetiEZ("right_BH_n_c");
  if (n != INT_MAX)   nlb[2] = n;
    
  if(nlb[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(nlb[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(nlb[2] == INT_MAX)
    Error0("n_c could not be set.\n");
  
  /* adding the results to the parameter data base */
  
  /* n_a, n_b, n_c */
  /* right box */
  sprintf(par,"grid%u_right_central_box_n_a",gn);
  sprintf(val,"%u",nlb[0]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_central_box_n_b",gn);
  sprintf(val,"%u",nlb[1]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_central_box_n_c",gn);
  sprintf(val,"%u",nlb[2]);
  add_parameter_string(par,val);
  
  /* size a,b,c we take it the same as left box. no biggie! */
  sprintf(par,"grid%u_right_central_box_size_a",gn);
  Psetd(par,box_size_r);
  
  sprintf(par,"grid%u_right_central_box_size_b",gn);
  Psetd(par,box_size_r);
  
  sprintf(par,"grid%u_right_central_box_size_c",gn);
  Psetd(par,box_size_r);
  
  /* R1 and R2 outermost */
  sprintf(par,"grid%u_outermost%u_R2",gn,0);
  Psetd(par,R_outermost[0]);
    
  for (i = 1; i < N_Outermost_Split; i++)
  {
    /* R1: */
    sprintf(par,"grid%u_outermost%u_R1",gn,i);
    Psetd(par,R_outermost[i-1]);
    
    /* R2: */
    sprintf(par,"grid%u_outermost%u_R2",gn,i);
    Psetd(par,R_outermost[i]);
    
  }
  
  /* assuming the center of left NS at (0,-S/2,0) */
  sprintf(par,"grid%u_left_NS_center_a",gn);
  Psetd(par,0.0);
  
  sprintf(par,"grid%u_left_NS_center_b",gn);
  Psetd(par,-S/2);
  
  sprintf(par,"grid%u_left_NS_center_c",gn);
  Psetd(par,0.0);
  
  /* assuming the center of right BH at (0,S/2,0) */
  sprintf(par,"grid%u_right_BH_center_a",gn);
  Psetd(par,0.0);
  
  sprintf(par,"grid%u_right_BH_center_b",gn);
  Psetd(par,S/2);
  
  sprintf(par,"grid%u_right_BH_center_c",gn);
  Psetd(par,0.0);
  
  free(R_outermost);
  
  /* make new patches: */
  /* to optimize the routine check if we can use some of the
  // previous grid data for the next grid */
  const int change_res_flg = Pgeti("did_resolution_change?");
  const int change_NS_flg  = Pgeti("did_NS_surface_change?");
  const int change_AH_flg  = Pgeti("did_AH_surface_change?");
  
  Pseti("use_previous_data",0);
  
  /* either the resolution is changed or it is the first grid */
  if (change_res_flg || !grid_prev)/* make geometry from scratch */
  {
    printf("~> Making patches from scratch ...\n");
    
    make_patches(grid_next);/* making patch(es) to cover the grid */
    realize_geometry(grid_next);/* realizing the geometry of whole grid
                     // including the way patches have been sewed,
                     // normal to the boundary, 
                     // outer-boundary, inner boundary and etc. */
  }
  /* NS surface and AH surface are not changed 
  // so let's use the previous grid */
  else if (!change_NS_flg && !change_AH_flg)
  {
    if (!strcmp_i(kind,"BBN_CubedSpherical_grid"))
      Error0(NO_OPTION);
    
    printf("~> Using BH, NS, filling_box and outermost of previous patches ...\n");
    
    free_grid(grid_next);
    grid_next     = grid_prev;
    grid_next->gn = gn;
    grid_prev     = 0;
    GridParams->grid_prev = 0;
    
    /* update the patches name */
    FOR_ALL_PATCHES(p,grid_next)
    {
      Patch_T *patch2 = grid_next->patch[p];
      char *name      = strstr(patch2->name,"_");/* name = "_*" */
      char name2[1000];
      
      sprintf(name2,"grid%u%s",gn,name);
      free(patch2->name);
      patch2->name = dup_s(name2);
    }
    
    /* remove the patches inside the BH */
    FOR_ALL_PATCHES(p,grid_next)
    {
      Patch_T *patch2 = grid_next->patch[p];
      
      if (!IsItInsideBHPatch(patch2))
        continue;
      
      free_patch(patch2);
    }
    grid_next->np    -= 7;/* 6 cubed spherical + 1 box */
    const unsigned np = grid_next->np;
    grid_next->patch  = 
        realloc(grid_next->patch,(np+1)*sizeof(*grid_next->patch));
    grid_next->patch[np] = 0;
    
    Pseti("use_previous_data",1);
  }
  /* only NS surface is not changed */
  else if (!change_NS_flg)
  {
    printf("~> Using NS, filling_box and outermost of previous patches ...\n");
    
    make_patches(grid_next);/* making patch(es) to cover the grid */
    /* since the resolution is not changed copy the geometry */
    move_geometry(grid_next,grid_prev);
    
    /* since the geometry NS patches are pristine
    // move patch->solving_man->jacobian */
    FOR_ALL_PATCHES(p,grid_next)
    {
      Patch_T *patch2 = grid_next->patch[p];
      Patch_T *patch1 = 0;
      char *stem      = 0;
      
      if (
          !IsItNSPatch(patch2)            && 
          !IsItNSSurroundingPatch(patch2) &&
          !IsItOutermostPatch(patch2)     &&/* outermost patches are always pristine */
          !IsItFillingBoxPatch(patch2)      /* filling box patches are always pristine */
          )
        continue;
    
      stem   = strstr(patch2->name,"_");
      stem++;
      patch1 = GetPatch(stem,grid_prev);
      
      move_solve_man_jacobian(patch2,patch1);
    }
    
  }
  /* only BH surface is not changed */
  else if (!change_AH_flg)
  {
    printf("~> Using BH, filling_box and outermost of previous patches ...\n");
   
    make_patches(grid_next);/* making patch(es) to cover the grid */
    /* since the resolution is not changed copy the geometry */
    move_geometry(grid_next,grid_prev);
    
    /* since the geometry BH surrounding patches are pristine
    // move patch->solving_man->jacobian */
    FOR_ALL_PATCHES(p,grid_next)
    {
      Patch_T *patch2 = grid_next->patch[p];
      Patch_T *patch1 = 0;
      char *stem      = 0;
      
      if (
          !IsItHorizonPatch(patch2)    &&
          !IsItOutermostPatch(patch2)  &&/* outermost patches are always pristine */
          !IsItFillingBoxPatch(patch2)   /* filling box patches are always pristine */
         )
        continue;
    
      stem   = strstr(patch2->name,"_");
      stem++;
      patch1 = GetPatch(stem,grid_prev);
      
      move_solve_man_jacobian(patch2,patch1);
    }
    
  }
  else/* if both NS and AH are changed but the resolution */
  {
    printf("~> Using filling_box and outermost of previous patches ...\n");

    make_patches(grid_next);/* making patch(es) to cover the grid */
    /* since the resolution is not changed copy the geometry */
    move_geometry(grid_next,grid_prev);
    
    /* since the geometry of filling box and outermost patches are pristine
    // move patch->solving_man->jacobian */
    FOR_ALL_PATCHES(p,grid_next)
    {
      Patch_T *patch2 = grid_next->patch[p];
      Patch_T *patch1 = 0;
      char *stem      = 0;
      
      if (
          !IsItOutermostPatch(patch2)  &&/* outermost patches are always pristine */
          !IsItFillingBoxPatch(patch2)   /* filling box patches are always pristine */
         )
        continue;
    
      stem   = strstr(patch2->name,"_");
      stem++;
      patch1 = GetPatch(stem,grid_prev);
      
      move_solve_man_jacobian(patch2,patch1);
    }
  }
    
  printf("} Creating grid ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
  return grid_next;
}

/* move the patch1->solving_man->jacobian to 
// patch2->solving_man->jacobian and put patch1->solving_man->jacobian = 0 */
static void move_solve_man_jacobian(Patch_T *const patch2,Patch_T *const patch1)
{
  assert(patch2);
  
  if (!patch2->solving_man)
  {
    patch2->solving_man = calloc(1,sizeof(*patch2->solving_man));
    IsNull(patch2->solving_man);
  }
  
  patch2->solving_man->jacobian = patch1->solving_man->jacobian;
  patch2->solving_man->nj       = patch1->solving_man->nj;
  patch1->solving_man->jacobian = 0;
  patch1->solving_man->nj       = 0;
}

/* move the geometry from previous grid to the next one and empty 
// the previous geometry, i.e. put its interface pointer to 0 */
static void move_geometry(Grid_T *const grid_next,Grid_T *const grid_prev)
{
  unsigned p2,p1;
  unsigned f,sf;
  
  assert(grid_next);
  
  if (!strcmp_i(grid_prev->kind,"BBN_CubedSpherical_grid"))
    Error0("This algorithm only works for BBN cubed spherical grid!\n");
  
  FOR_ALL_PATCHES(p2,grid_next)
  {
    Patch_T *patch1 = 0;
    Patch_T *patch2 = grid_next->patch[p2];
    char *name2     = strstr(patch2->name,"_");
    
    /* find the corresponding patch */
    FOR_ALL_PATCHES(p1,grid_prev)
    {
      patch1 = grid_prev->patch[p1];
      char *name1  = strstr(patch1->name,"_");
  
      if (!strcmp(name2,name1))
        break;
    }/* FOR_ALL_PATCHES(p2,grid_prev) */
    assert(patch1);
    
    /* move geometry */
    patch2->interface   = patch1->interface;
    patch2->innerB      = patch1->innerB;
    patch2->outerB      = patch1->outerB;
    patch2->is_a_closed = patch1->is_a_closed;
    patch2->is_b_closed = patch1->is_b_closed;
    patch2->is_c_closed = patch1->is_c_closed;
    patch1->interface   = 0;
                          
    /* update the internal pointers */
    Interface_T **face = patch2->interface;
    
    /* for all interfaces */
    FOR_ALL(f,face)
    {
      face[f]->patch = patch2;
      
      assert(!face[f]->point);/* why is it not empty??? */
      
      /* for all subfaces */
      for (sf = 0; sf < face[f]->ns; ++sf)
      {
        SubFace_T *subf = face[f]->subface[sf];
        subf->patch = patch2;
      }
    }

  }/* end of FOR_ALL_PATCHES(p2,grid_next) */
}

/* making  NS and BH surfaces function */
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  printf("{ Populating surface function for NS and BH ...\n");
  
  const double Max_R_NS_l = GridParams->Max_R_NS_l;/* maximum radius of NS */
  const double R_BH_r     = GridParams->R_BH_r;
  const double a_BH       = GridParams->a_BH;
  const double y_CM       = Pgetd("y_CM");
  const double C_BH       = 0.5*Pgetd("BH_NS_separation");/* center of BH patch it's on +y axis */
  const double Omega_BHNS = Pgetd("BH_NS_angular_velocity");
  const double g2         = 1-Pow2(-Omega_BHNS*(C_BH-y_CM));/* inverse square of Lorentz factor  */
  const double BH_center[3] = {Pgetd("BH_center_x"),Pgetd("BH_center_y")-C_BH,Pgetd("BH_center_z")};
  const double W1         = Pgetd("NS_surface_update_weight");
  const double W2         = 1-W1;
  double *R,*R_update;
  const double *R_old,*R_new;
  char par[1000] = {'\0'};
  unsigned N[3],n,i,j,k,N_total,ijk;
  Patch_T patch[1] = {0};
  struct Collocation_s coll_s[2] = {0};
  double X[3],r;
  const int same_res_flag = !Pgeti("did_resolution_change?");
  double dR_sum_square = 0;
  
  /* left NS */
  
  /* filling N */
  N[0] = (unsigned)PgetiEZ("n_a");
  N[1] = (unsigned)PgetiEZ("n_b");
  N[2] = (unsigned)PgetiEZ("n_c");
  /* check for override */
  n = (unsigned)PgetiEZ("left_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)PgetiEZ("left_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)PgetiEZ("left_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    Error0("n_c could not be set.\n");
    
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
    
    X[2] = 1;/* since we are on the NS surface */
    
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
    
    if (same_res_flag)
    {
      Grid_T *grid_prev      = GridParams->grid_prev;
      Patch_T *patch_prev    = GetPatch("left_NS_up",grid_prev);
      const int R0_ind       = LookUpField_E("surface_function",patch_prev);
      const double *const R0 = patch_prev->pool[R0_ind]->v;
      
      for (i = 0; i < N[0]; ++i)
        for (j = 0; j < N[1]; ++j)
        {
          ijk = L(N,i,j,0);
          dR_sum_square += Pow2(1-R[ijk]/R0[ijk]);
        }
      
    }/* end of if (same_res_flag) */

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
    
    if (same_res_flag)
    {
      Grid_T *grid_prev      = GridParams->grid_prev;
      Patch_T *patch_prev    = GetPatch("left_NS_down",grid_prev);
      const int R0_ind       = LookUpField_E("surface_function",patch_prev);
      const double *const R0 = patch_prev->pool[R0_ind]->v;
      
      for (i = 0; i < N[0]; ++i)
        for (j = 0; j < N[1]; ++j)
        {
          ijk = L(N,i,j,0);
          dR_sum_square += Pow2(1-R[ijk]/R0[ijk]);
        }
      
    }/* end of if (same_res_flag) */

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
    
    if (same_res_flag)
    {
      Grid_T *grid_prev      = GridParams->grid_prev;
      Patch_T *patch_prev    = GetPatch("left_NS_back",grid_prev);
      const int R0_ind       = LookUpField_E("surface_function",patch_prev);
      const double *const R0 = patch_prev->pool[R0_ind]->v;
      
      for (i = 0; i < N[0]; ++i)
        for (j = 0; j < N[1]; ++j)
        {
          ijk = L(N,i,j,0);
          dR_sum_square += Pow2(1-R[ijk]/R0[ijk]);
        }
      
    }/* end of if (same_res_flag) */

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
    
    if (same_res_flag)
    {
      Grid_T *grid_prev      = GridParams->grid_prev;
      Patch_T *patch_prev    = GetPatch("left_NS_front",grid_prev);
      const int R0_ind       = LookUpField_E("surface_function",patch_prev);
      const double *const R0 = patch_prev->pool[R0_ind]->v;
      
      for (i = 0; i < N[0]; ++i)
        for (j = 0; j < N[1]; ++j)
        {
          ijk = L(N,i,j,0);
          dR_sum_square += Pow2(1-R[ijk]/R0[ijk]);
        }
      
    }/* end of if (same_res_flag) */

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
    
    if (same_res_flag)
    {
      Grid_T *grid_prev      = GridParams->grid_prev;
      Patch_T *patch_prev    = GetPatch("left_NS_left",grid_prev);
      const int R0_ind       = LookUpField_E("surface_function",patch_prev);
      const double *const R0 = patch_prev->pool[R0_ind]->v;
      
      for (i = 0; i < N[0]; ++i)
        for (j = 0; j < N[1]; ++j)
        {
          ijk = L(N,i,j,0);
          dR_sum_square += Pow2(1-R[ijk]/R0[ijk]);
        }
      
    }/* end of if (same_res_flag) */

    
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
    
    if (same_res_flag)
    {
      Grid_T *grid_prev      = GridParams->grid_prev;
      Patch_T *patch_prev    = GetPatch("left_NS_right",grid_prev);
      const int R0_ind       = LookUpField_E("surface_function",patch_prev);
      const double *const R0 = patch_prev->pool[R0_ind]->v;
      
      for (i = 0; i < N[0]; ++i)
        for (j = 0; j < N[1]; ++j)
        {
          ijk = L(N,i,j,0);
          dR_sum_square += Pow2(1-R[ijk]/R0[ijk]);
        }
      
    }/* end of if (same_res_flag) */

  }
  else
    Error0(NO_OPTION);
  
  free(R);
  
  /* should we change the NS surface */
  const double NS_surf_tolerance = Pgetd("NS_surface_tolerance");
  const double dR_rms = sqrt(dR_sum_square);
  const int    did_NS_surface_finder_work = Pgeti("did_NS_surface_finder_work?");
  /* if dR_rms is small           or 
  // if surface finder was failed 
  // do not change the NS surface function */
  if (
      (same_res_flag && GridParams->grid_prev) && 
      (dR_rms < NS_surf_tolerance || !did_NS_surface_finder_work) 
     )
  {
    if (dR_rms < NS_surf_tolerance && did_NS_surface_finder_work)/* prints only if this is the case */
      printf("~> |R_2 - R1|/|R_1| = %g < Tolerance = %g:\n",dR_rms,NS_surf_tolerance);
    if (!did_NS_surface_finder_work)
      printf("~> NS surface finder was failed:\n");
    printf("~> No changes in NS surface\n");

    Pseti("did_NS_surface_change?",0);
    
    /* update the surface function accordingly: */
    sprintf(par,"grid%u_left_NS_surface_function_up",grid->gn-1);
    R = Pgetdd(par);
    sprintf(par,"grid%u_left_NS_surface_function_up",grid->gn);
    update_parameter_array(par,R,N_total);
    
    sprintf(par,"grid%u_left_NS_surface_function_down",grid->gn-1);
    R = Pgetdd(par);
    sprintf(par,"grid%u_left_NS_surface_function_down",grid->gn);
    update_parameter_array(par,R,N_total);
    
    sprintf(par,"grid%u_left_NS_surface_function_back",grid->gn-1);
    R = Pgetdd(par);
    sprintf(par,"grid%u_left_NS_surface_function_back",grid->gn);
    update_parameter_array(par,R,N_total);
    
    sprintf(par,"grid%u_left_NS_surface_function_front",grid->gn-1);
    R = Pgetdd(par);
    sprintf(par,"grid%u_left_NS_surface_function_front",grid->gn);
    update_parameter_array(par,R,N_total);
    
    sprintf(par,"grid%u_left_NS_surface_function_left",grid->gn-1);
    R = Pgetdd(par);
    sprintf(par,"grid%u_left_NS_surface_function_left",grid->gn);
    update_parameter_array(par,R,N_total);
    
    sprintf(par,"grid%u_left_NS_surface_function_right",grid->gn-1);
    R = Pgetdd(par);
    sprintf(par,"grid%u_left_NS_surface_function_right",grid->gn);
    update_parameter_array(par,R,N_total);
    
  }
  /* else it is changed */
  else
  {
    Pseti("did_NS_surface_change?",1);
    
    if (GridParams->grid_prev)
    {
      if (dR_rms > 0)
        printf("~> |R_2 - R1|/|R_1| = %g >= Tolerance = %g:\n",dR_rms,NS_surf_tolerance);
      
      /* update the surface in relaxed fashion */
      if (same_res_flag && did_NS_surface_finder_work)
      {
        /* update the surface function accordingly: */
        sprintf(par,"grid%u_left_NS_surface_function_up",grid->gn-1);
        R_old = Pgetdd(par);
        sprintf(par,"grid%u_left_NS_surface_function_up",grid->gn);
        R_new = Pgetdd(par);
        R_update = alloc_double(N_total);
        for (i = 0; i < N[0]; ++i)
          for (j = 0; j < N[1]; ++j)
            for (k = 0; k < N[2]; ++k)
            {
              ijk = L(N,i,j,k);
              R_update[ijk] = W1*R_new[ijk]+W2*R_old[ijk];
            }
            
        update_parameter_array(par,R_update,N_total);
        free(R_update);
        
        sprintf(par,"grid%u_left_NS_surface_function_down",grid->gn-1);
        R_old = Pgetdd(par);
        sprintf(par,"grid%u_left_NS_surface_function_down",grid->gn);
        R_new = Pgetdd(par);
        R_update = alloc_double(N_total);
        for (i = 0; i < N[0]; ++i)
          for (j = 0; j < N[1]; ++j)
            for (k = 0; k < N[2]; ++k)
            {
              ijk = L(N,i,j,k);
              R_update[ijk] = W1*R_new[ijk]+W2*R_old[ijk];
            }
            
        update_parameter_array(par,R_update,N_total);
        free(R_update);
                
        sprintf(par,"grid%u_left_NS_surface_function_back",grid->gn-1);
        R_old = Pgetdd(par);
        sprintf(par,"grid%u_left_NS_surface_function_back",grid->gn);
        R_new = Pgetdd(par);
        R_update = alloc_double(N_total);
        for (i = 0; i < N[0]; ++i)
          for (j = 0; j < N[1]; ++j)
            for (k = 0; k < N[2]; ++k)
            {
              ijk = L(N,i,j,k);
              R_update[ijk] = W1*R_new[ijk]+W2*R_old[ijk];
            }
            
        update_parameter_array(par,R_update,N_total);
        free(R_update);
                
        sprintf(par,"grid%u_left_NS_surface_function_front",grid->gn-1);
        R_old = Pgetdd(par);
        sprintf(par,"grid%u_left_NS_surface_function_front",grid->gn);
        R_new = Pgetdd(par);
        R_update = alloc_double(N_total);
        for (i = 0; i < N[0]; ++i)
          for (j = 0; j < N[1]; ++j)
            for (k = 0; k < N[2]; ++k)
            {
              ijk = L(N,i,j,k);
              R_update[ijk] = W1*R_new[ijk]+W2*R_old[ijk];
            }
            
        update_parameter_array(par,R_update,N_total);
        free(R_update);
                
        sprintf(par,"grid%u_left_NS_surface_function_left",grid->gn-1);
        R_old = Pgetdd(par);
        sprintf(par,"grid%u_left_NS_surface_function_left",grid->gn);
        R_new = Pgetdd(par);
        R_update = alloc_double(N_total);
        for (i = 0; i < N[0]; ++i)
          for (j = 0; j < N[1]; ++j)
            for (k = 0; k < N[2]; ++k)
            {
              ijk = L(N,i,j,k);
              R_update[ijk] = W1*R_new[ijk]+W2*R_old[ijk];
            }
            
        update_parameter_array(par,R_update,N_total);
        free(R_update);
                
        sprintf(par,"grid%u_left_NS_surface_function_right",grid->gn-1);
        R_old = Pgetdd(par);
        sprintf(par,"grid%u_left_NS_surface_function_right",grid->gn);
        R_new = Pgetdd(par);
        R_update = alloc_double(N_total);
        for (i = 0; i < N[0]; ++i)
          for (j = 0; j < N[1]; ++j)
            for (k = 0; k < N[2]; ++k)
            {
              ijk = L(N,i,j,k);
              R_update[ijk] = W1*R_new[ijk]+W2*R_old[ijk];
            }
            
        update_parameter_array(par,R_update,N_total);
        free(R_update);
        
      }
    }
    printf("~> Update NS surface\n");
  }
  
  R = 0;
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
  N[0] = (unsigned)PgetiEZ("n_a");
  N[1] = (unsigned)PgetiEZ("n_b");
  N[2] = (unsigned)PgetiEZ("n_c");
  
  /* check for override */
  n = (unsigned)PgetiEZ("right_BH_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)PgetiEZ("right_BH_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)PgetiEZ("right_BH_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    Error0("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    Error0("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    Error0("n_c could not be set.\n");
   
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
    Error0("This surface function is incomplete;\n"
            "One needs to incorporate changes in BH center and in center of mass of the system.\n");
    /* surface up */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        r = sqrt(
                 (1+Pow2(X[0])+Pow2(X[1]))/
                 ((g2*Pow2(X[0])+Pow2(X[1]))/(Pow2(R_BH_r)+Pow2(a_BH)) + 1/Pow2(R_BH_r))
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
                 (1+Pow2(X[0])+Pow2(X[1]))/
                 ((Pow2(X[0])+g2*Pow2(X[1]))/(Pow2(R_BH_r)+Pow2(a_BH)) + 1/Pow2(R_BH_r))
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
                 (1+Pow2(X[0])+Pow2(X[1]))/
                 (((g2+Pow2(X[1])))/(Pow2(R_BH_r)+Pow2(a_BH)) + Pow2(X[0])/Pow2(R_BH_r))
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
                 (1+Pow2(X[0])+Pow2(X[1]))/
                 (((g2+Pow2(X[0])))/(Pow2(R_BH_r)+Pow2(a_BH)) + Pow2(X[1])/Pow2(R_BH_r))
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
                 (1+Pow2(X[0])+Pow2(X[1]))/
                 (((1+g2*Pow2(X[0])))/(Pow2(R_BH_r)+Pow2(a_BH)) + Pow2(X[1])/Pow2(R_BH_r))
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
                 (1+Pow2(X[0])+Pow2(X[1]))/
                 (((1+g2*Pow2(X[1])))/(Pow2(R_BH_r)+Pow2(a_BH)) + Pow2(X[0])/Pow2(R_BH_r))
                );
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_right",grid->gn);
    add_parameter_array(par,R,N_total);
  }
  else if (/* if BH center is instructed to change */
           (strstr_i(Pgets("P_ADM_control_method"),"Px_BH_center")  ||
            strstr_i(Pgets("P_ADM_control_method"),"Py_BH_center")
           )                                                        &&
            strcmp_i(GridParams->BH_R_type,"PerfectSphere")
          )
  {
    if (Pcmps("extrapolate_inside_BH_method","Ylm"))
      Error0("This part is not complete for Ylm BH extrapolation.\n"
              "You should find r_bh_max appropriately.\n");
      
    double r_bh_max = 0;/* this is for Ylm extrapolation */
    
    /* surface up */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        r = AH_surface_function_PerfectSphere_CS(X[0],X[1],R_BH_r,BH_center,UP);
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
        r = AH_surface_function_PerfectSphere_CS(X[0],X[1],R_BH_r,BH_center,DOWN);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_down",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface back */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        r = AH_surface_function_PerfectSphere_CS(X[0],X[1],R_BH_r,BH_center,BACK);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_back",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface front */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        r = AH_surface_function_PerfectSphere_CS(X[0],X[1],R_BH_r,BH_center,FRONT);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_front",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface left */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        r = AH_surface_function_PerfectSphere_CS(X[0],X[1],R_BH_r,BH_center,LEFT);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_left",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface right */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        r = AH_surface_function_PerfectSphere_CS(X[0],X[1],R_BH_r,BH_center,RIGHT);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_right_BH_surface_function_right",grid->gn);
    add_parameter_array(par,R,N_total);
    
    Psets("BH_R_type","PerfectSphere");
    Psetd("BH_R_size",r_bh_max);/* this is used later for Ylm extrapolation */
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
    
    Psets("BH_R_type","PerfectSphere");
    Psetd("BH_R_size",R_BH_r);/* this is used later for Ylm extrapolation */
  }
  else
    Error0(NO_OPTION);
  
  free(R);
  
  printf("} Populating surface function for NS and BH ==> Done.\n");
}

/* given (X,Y,Z) in the specified slice of NS in cubed spherical coords
// it finds the associated polar and azimuthal angels on the surface of NS */
static void find_theta_phi_of_XYZ_NS_CS(double *const theta,double *const phi,const double *const X,const Flag_T side)
{
  const double a = X[0];
  const double b = X[1];
  const double d = sqrt(1+Pow2(a)+Pow2(b));
  
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
      Error0(NO_OPTION);
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
  if (!patch->pool[h_ind]->v)/* if there is no enthalpy defined in the patch */
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
  GridParams->NS_R_type = Pgets("NS_surface_finder_method");
  
  /* find NS surface using cubed spherical points */
  //if (strstr_i(GridParams->NS_R_type,"CubedSpherical"))
    //find_NS_surface_CS_method_CS(grid,GridParams);
    
  /* find NS surface using spherical harmonic points */
  if (strstr_i(GridParams->NS_R_type,"SphericalHarmonic"))
    find_NS_surface_Ylm_method_CS(grid,GridParams);
  else
    Error0(NO_OPTION);
}

/* extrapolate fluid fields outside NS */
static void extrapolate_fluid_fields_outsideNS(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("{ Extrapolating fluid fields outside NS ...\n");
  unsigned p;
  
  if (strcmp_i(grid->kind,"BBN_CubedSpherical_grid"))
  {
    if (Pcmps("extrapolate_fluid_fields_method","phi:exp_continuity,enthalpy:Ylm"))
    {
      /* extrapolate phi and W => but don't make enthalpy */
      extrapolate_outsideNS_CS_exp_continuity_method(grid);
      /* extrapolate enthalpy  */
      extrapolate_outsideNS_CS_Ylm_method(grid,"enthalpy");
      
      /* calculating the derivatives of the enthalpy */
      FOR_ALL_PATCHES(p,grid)
      {
        /* surrounding patch */
        Patch_T *patch = grid->patch[p];
        if (!IsItNSSurroundingPatch(patch))
          continue;
        
        Field_T *enthalpy = patch->pool[Ind("enthalpy")];
        DECLARE_AND_EMPTY_FIELD(denthalpy_D2)
        DECLARE_AND_EMPTY_FIELD(denthalpy_D1)
        DECLARE_AND_EMPTY_FIELD(denthalpy_D0)
        denthalpy_D2->v = Partial_Derivative(enthalpy,"z");
        denthalpy_D1->v = Partial_Derivative(enthalpy,"y");
        denthalpy_D0->v = Partial_Derivative(enthalpy,"x");
      }

    }
    else
      Error0(NO_OPTION);
      
    if (0)/* this method is sucks! */
      extrapolate_outsideNS_CS_slop_method(grid);
 
  }
  else
    Error0(NO_OPTION);
  
  printf("} Extrapolating fluid fields outside NS ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* populating surface function of the apparent horizon for
// perfect sphere with radius R and arbitrary center. 
// note: it assumes the origin of the coordinate system is at (0,0,0):
// (x-c[0])^2+(y-c[1])^2+(z-c[2])^2 = R^2.
// ->return value: surface function value at the give specific (a,b) cubed spherical coords. */
static double AH_surface_function_PerfectSphere_CS(const double a,const double b,const double R,const double *const c,const Flag_T side)
{
  const double x0 = c[0];
  const double y0 = c[1];
  const double z0 = c[2];
  double S = 0;
  
  switch (side)
  {
    case UP:
      S = a*x0 + b*y0 + z0 + Sqrt(Power(a*x0 + b*y0 + z0,2) - 
     (1 + Power(a,2) + Power(b,2))*
      (-Power(R,2) + Power(x0,2) + Power(y0,2) + Power(z0,2)));
    break;
    case DOWN:
      S = -(b*x0) - a*y0 + z0 - Sqrt(Power(b*x0 + a*y0 - z0,2) - 
     (1 + Power(a,2) + Power(b,2))*
      (-Power(R,2) + Power(x0,2) + Power(y0,2) + Power(z0,2)));
    break;
    case LEFT:
      S = -(a*x0) + y0 - b*z0 - Sqrt(Power(a*x0 - y0 + b*z0,2) - 
     (1 + Power(a,2) + Power(b,2))*
      (-Power(R,2) + Power(x0,2) + Power(y0,2) + Power(z0,2)));
    break;
    case RIGHT:
      S = b*x0 + y0 + a*z0 + Sqrt(Power(b*x0 + y0 + a*z0,2) - 
     (1 + Power(a,2) + Power(b,2))*
      (-Power(R,2) + Power(x0,2) + Power(y0,2) + Power(z0,2)));
    break;
    case BACK:
      S = x0 - b*y0 - a*z0 - Sqrt(Power(-x0 + b*y0 + a*z0,2) - 
     (1 + Power(a,2) + Power(b,2))*
      (-Power(R,2) + Power(x0,2) + Power(y0,2) + Power(z0,2)));
    break;
    case FRONT:
      S = x0 + a*y0 + b*z0 + Sqrt(Power(x0 + a*y0 + b*z0,2) - 
     (1 + Power(a,2) + Power(b,2))*
      (-Power(R,2) + Power(x0,2) + Power(y0,2) + Power(z0,2)));
    break;
    default:
      Error0(NO_JOB);
  }/* end of switch */

  return fabs(S)/sqrt(Pow2(a)+Pow2(b)+1);
}

/* free given grid and parameters related to the given grid */
extern Parameter_T **parameters_global;
void bbn_free_grid_and_its_parameters(Grid_T *grid)
{
  if (!grid)/* if grid is empty do nothing */
    return;
    
  const int keep_grid = Pgeti("use_previous_data");
  char suffix[100] = {'\0'};
  unsigned i,np,gn;
  
  if (!keep_grid)
  {
    gn = grid->gn;
    free_grid(grid);
  }
  else
  {
    gn   = grid->gn-1;
    grid = 0;
  }
    
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
  
  sprintf(suffix,"grid%u_",gn);/* parameters related to this grid */
  for (i = 0; i < np;)/* no increment */
  {
    if (strstr(parameters_global[i]->lv,suffix))
    {
      /* note: the last par is put in palce of removed par
      // so don't increment i */
      free_parameter(parameters_global[i]->lv);
      np--;
    }
    else
      i++;
  }
  
}
