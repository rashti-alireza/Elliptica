/*
// Alireza Rashti
// February 2020
*/

#include "bbn_checkpoint.h"


/* if you set a parameter as below in the par file it uses this parameter
// and disregards its value obtained from the checkpoint:
// e.g.
// modify_checkpoint_par: n_a = 4(x6)
// modify_checkpoint_par: Solving_Max_Number_of_Iteration = 0  */
static unsigned n_modified_checkpoint_par;/* number of modify_checkpoint_par */
static Parameter_T **modified_checkpoint_par;/* modified pars in par file
                                           // to be used after loading of
                                           // the checkpoint file. */

/* write checkpoint for the given grid 
// NOTE: the order of writing and reading is crucial */
void bbn_write_checkpoint(Grid_T *const grid)
{
  /* print some descriptions */
  pr_line_custom('=');
  printf("{ Writing checkpoint file ...\n");
  
  if (!grid)
  {
    printf("~> The given grid is empty.\n");
    printf("} Writing checkpoint ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
    return;
  }
  
  FILE *file = 0;
  const char *const out_dir = Pgets("iteration_output");
  const unsigned sol_it_n = (unsigned)PgetiEZ("solving_iteration_number");
  char file_path[MAX_ARR];
  char msg[MAX_ARR];
  char *const p_msg = msg;/* to avoid GCC warning for FWriteP_bin */
  const double dt  = Pgetd("write_checkpoint_every");/* unit is hour */
  const double now = get_time_sec()/(3600);
  static double last_checkpoint_was = 0;/* the time where the last 
                                        // checkpoint happened in hours */
  
  if (dt+last_checkpoint_was > now && !Pgeti("STOP"))
  {
    printf("~> It's early for writing checkpoint.\n");
    printf("} Writing checkpoint ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
    return;
  }
  
  last_checkpoint_was = now;
  
  /* NOTE: the order of writing is super crucial for reading */
  
  /* write header */
  write_header(grid);
  
  /* write all parameters in the checkpoint file */
  write_parameters(grid);
  
  /* write all fields value in the checkpoint file */
  write_fields(grid);
  
  /* successful message at the end of the checkpoint file */
  sprintf(file_path,"%s/%s_temp",out_dir,CHECKPOINT_FILE_NAME);
  file = fopen(file_path,"a");
  pointerEr(file);
  sprintf(msg,"%s",END_MSG);
  FWriteP_bin(p_msg,strlen(msg)+1);
  fclose(file);
  
  /* replace checkpoint file with the previous */
  move_checkpoint_file();
  
  /* write bbn properties mainly used for id reader */
  bbn_print_properties(grid,sol_it_n,out_dir,"w",0);
  
  printf("} Writing checkpoint file ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* ->return value: if the chekpoint is completed 1, otherwise 0. */
int bbn_IsCheckpointFileCompleted(const char *const file_path)
{
  int ret = 0;
  FILE *file;
  char msg[MAX_ARR];
  int msg_len = (int)strlen(END_MSG)+1;
  
  file = fopen(file_path,"r");
  pointerEr(file);
  
  fseek(file,-msg_len,SEEK_END);
  assert(fread(msg,(unsigned)msg_len,1,file));
  
  if (strstr(msg,END_MSG))
    ret = 1;
  else
    ret = 0;
      
  fclose(file);
  
  return ret;
}

/* write header of checkpoint file */
static void write_header(const Grid_T *const grid)
{
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  unsigned np;

  sprintf(file_path,"%s/%s_temp",folder,CHECKPOINT_FILE_NAME);
  if (!access(file_path,F_OK))/* if file exists */
    Error0("File already exists.\n");
  
    
  file = fopen(file_path,"w");
  pointerEr(file);
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
    
  /* NO white spaces since I'll use fscanf */
  fprintf(file,"%s\n",ALLOC_HEADER);
  fprintf(file,"number_of_parameters=%u\n",np);
  //fprintf(file,"number_of_patches=%u\n",grid->np);
  fprintf(file,"grid_number=%u\n",grid->gn);
  fprintf(file,"grid_kind=%s\n",grid->kind);
  fprintf(file,"%s\n",ALLOC_FOOTER);
  fclose(file);
}

/* replace checkpoint file with the previous one */
static void move_checkpoint_file(void)
{
  const char *const folder = Pgets("iteration_output");
  char command[2*MAX_ARR];
  
  sprintf(command,"mv %s/%s_temp %s/%s",
          folder,CHECKPOINT_FILE_NAME,folder,CHECKPOINT_FILE_NAME);
  shell_command(command);
  
  printf("checkpoint file path:\n%s/%s\n",folder,CHECKPOINT_FILE_NAME);
  fflush(stdout);
}

/* write all of the pertinent parameters in the checkpoint file */
static void write_parameters(const Grid_T *const grid)
{
  printf ("~> Writing parameters in checkpoint file ...\n");
  fflush(stdout);
  
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  char title_line[MAX_ARR] = {'\0'};
  char *const p_title_line = title_line;/* to avoid GCC warning for FWriteP_bin */
  unsigned i,np;

  sprintf(file_path,"%s/%s_temp",folder,CHECKPOINT_FILE_NAME);
  file = fopen(file_path,"ab");
  pointerEr(file);
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
    
  /* NOTE the order is super crucial for reading part */
  sprintf(title_line,"%s",PARAM_HEADER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  
  for (i = 0; i < np; ++i)
  {
    Parameter_T *p = parameters_global[i];
    
    FWriteP_bin(p->lv,strlen(p->lv)+1);
    FWriteP_bin(p->rv,strlen(p->rv)+1);
    FWriteP_bin(p->rv_ip,strlen(p->rv_ip)+1);
    FWriteV_bin(p->rv_double,1);
    FWriteP_bin(p->rv_array,p->rv_n);
    FWriteV_bin(p->rv_n,1);
    FWriteV_bin(p->iterative,1);
    FWriteV_bin(p->double_flg,1);
  }
  
  sprintf(title_line,"%s",PARAM_FOOTER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  
  fclose(file);
  
  UNUSED(grid);
}

/* write all of the fields in the checkpoint file */
static void write_fields(const Grid_T *const grid)
{
  printf("~> Writing fields in checkpoint file ...\n");
  fflush(stdout);
  
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  char title_line[MAX_ARR] = {'\0'};
  char *const p_title_line = title_line;/* to avoid GCC warning for FWriteP_bin */
  unsigned p;
  
  sprintf(file_path,"%s/%s_temp",folder,CHECKPOINT_FILE_NAME);
  file = fopen(file_path,"ab");
  pointerEr(file);
  
  /* NOTE the order is crucial for reading part */
  
  sprintf(title_line,"%s",FIELD_HEADER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned f,count_nfld;
    
    if (IsItInsideBHPatch(patch))
      continue;
      
    /* count number fields we want to save */
    count_nfld = 0;
    for (f = 0; f < patch->nfld; ++f)
    {
      Field_T *field = patch->pool[f];
      if (DoSaveField(field))
        count_nfld++;
    }
    FWriteV_bin(count_nfld,1);
    
    for (f = 0; f < patch->nfld; ++f)
    {
      Field_T *field = patch->pool[f];
      if (DoSaveField(field))
      {
        FWriteP_bin(field->name,strlen(field->name)+1);
        FWriteP_bin(field->v,nn);
      }
    }
  }
  
  sprintf(title_line,"%s",FIELD_FOOTER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  
  fclose(file);  
}

/* ->return value: if field->name matches 1, otherwise 0 */
static int DoSaveField(const Field_T *const field)
{
  int ret = 0;
  const char *const kind = field->patch->grid->kind;
  const char *const name = field->name;
  
  if (strcmp_i(kind,"BBN_CubedSpherical_grid"))
  {
    if (
        !strcmp(name,"phi")   ||
        !strcmp(name,"psi")   ||
        !strcmp(name,"eta")   ||
        !strcmp(name,"B0_U0") ||
        !strcmp(name,"B0_U1") ||
        !strcmp(name,"B0_U2") ||
        !strcmp(name,"enthalpy")
        )
      ret = 1;
  }
  return ret;
}

/* read checkpoint file and creat grid and parameters accordingly */
Grid_T *bbn_init_from_checkpoint(FILE *const file)
{
  Grid_T *grid = 0;
  struct checkpoint_header alloc_info[1] = {0};
  
  /* reading the header for allocations */
  read_header(alloc_info,file);
  
  /* alloc grid and parameters */
  alloc_db(alloc_info);
  
  /* read parameters */
  read_parameters(alloc_info,file);
  
  grid = alloc_info->grid;
  
  /* make the patches */
  make_patches(grid);
  
  /* realizing the geometry */
  realize_geometry(grid);
  
  /* adding all of the fields needed for construction of Initial Data */
  bbn_add_fields(grid);
  
  /* populating the free data part of initial data that we chose ourself */
  bbn_populate_free_data(grid);
  
  /* read fields content from the checkpoint file */
  read_fields(alloc_info,file);
  
  /* initialzing some mediate field */
  init_mediate_field(grid);
  
  /* make normal vectorn on BH horizon 
  // note: this MUST be before "bbn_partial_derivatives_fields" */
  bbn_make_normal_vector_on_BH_horizon(grid);
  
  /* taking partial derivatives of the fields needed for equations */
  bbn_partial_derivatives_fields(grid);

  /* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S */
  bbn_update_stress_energy_tensor(grid,1);
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  bbn_update_Aij(grid);
  
  return grid;
}

/* initialzing some mediate field */
static void init_mediate_field(Grid_T *const grid)
{
  const unsigned np = grid->np;
  unsigned p;
  
  /* W_U[0-2],Beta_U[0-2],B1_U[0-2] */
  const double Omega_NS_x = Pgetd("NS_Omega_U0");
  const double Omega_NS_y = Pgetd("NS_Omega_U1");
  const double Omega_NS_z = Pgetd("NS_Omega_U2");
  const double C_NS = Pgetd("NS_Center_y");
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < np; ++p)
  {
    Patch_T *patch = grid->patch[p];
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

/* reading the header for allocations */
static void read_header(struct checkpoint_header *const alloc_info,FILE *const file)
{  
  printf("~> Reading checkpoint file header ...\n");
  fflush(stdout);
  
  char line[MAX_ARR] = {'\0'};

  fseek(file,0,SEEK_SET);
  
  /* check the header */
  fscanf(file,"%s",line);
  if (strcmp(line,ALLOC_HEADER))
      Error0("No header found. Checkpoint file got a problem.\n");
    
  /* read allocations */
  while (strcmp(line,ALLOC_FOOTER))
  {
    fscanf(file,"%s",line);
    if (!strcmp(line,ALLOC_FOOTER))
      break;
      
    char *v = strstr(line,"=");/* v -> "=..." */
    if (!v)
      Error0("No value found. Checkpoint file got a problem.\n");
    v++;
    
    if (strstr(line,"number_of_parameters"))
    {
      alloc_info->npar = (unsigned)atoi(v);
    }
    //else if (strstr(line,"number_of_patches"))
    //{
      //alloc_info->npatch = (unsigned)atoi(v);
    //}
    else if (strstr(line,"grid_number"))
    {
      alloc_info->grid_number = (unsigned)atoi(v);
    }
    else if (strstr(line,"grid_kind"))
    {
      alloc_info->grid_kind = dup_s(v);
    }
    else
      Error0(NO_OPTION);
  }
  /* read the binary parts */
  fseek(file,ftell(file)+1,SEEK_SET);/* +1 since fscanf won't read \n */
}

/* find and save modified checkpoint pars specified at the par file
// also delete these pars from the par db */
static void find_and_save_modified_checkpoint_pars(void)
{
  modified_checkpoint_par   = 0;/* global var in this file */
  n_modified_checkpoint_par = 0;/* global var in this file */
  const char *const keyword_prefix = "modify_checkpoint_par:";
  char str[MAX_ARR],*pstr;
  unsigned np,nmpar;
  
  /* find the modified pars and save them */
  np    = 0;
  nmpar = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
  {
    /* find */
    if (strstr_i(parameters_global[np]->lv,keyword_prefix))
    {
      /* save */
      modified_checkpoint_par = 
        realloc(modified_checkpoint_par,(nmpar+1)*sizeof(*modified_checkpoint_par));
      pointerEr(modified_checkpoint_par);
      modified_checkpoint_par[nmpar] = parameters_global[np];
      parameters_global[np] = 0;
      parameters_global[np] = calloc(1,sizeof(*parameters_global[np]));
      pointerEr(parameters_global[np]);
      
      /* trim the prefix */
      pstr = strstr(modified_checkpoint_par[nmpar]->lv,":");/* ->: */
      pstr++;
      assert(pstr);
      sprintf(str,"%s",pstr);
      free(modified_checkpoint_par[nmpar]->lv);
      modified_checkpoint_par[nmpar]->lv = dup_s(str);
      
      nmpar++;
    }
    
    np++;
  }
  
  n_modified_checkpoint_par = nmpar;
}

/* alloc parameters and grid */
static void alloc_db(struct checkpoint_header *const alloc_info)
{
  printf("~> Allocating parameters and patches ...\n");
  fflush(stdout);
  
  const unsigned grid_number = alloc_info->grid_number,
                 npar        = alloc_info->npar;
  Grid_T *grid = 0;
  unsigned i;
  
  /* find and save modified checkpoint pars specified at the pars file */
  find_and_save_modified_checkpoint_pars();
  
  /* allocate parameters */
  free_parameter_db();
  parameters_global = calloc(npar+1,sizeof(*parameters_global));
  pointerEr(parameters_global);
  for (i = 0; i < npar; ++i)
  {
    parameters_global[i] = calloc(1,sizeof(*parameters_global[i]));
    pointerEr(parameters_global[i]);
  }
  
  /* allocate grid */
  free_grid_db();
  grid       = alloc_grid();
  grid->gn   = grid_number;
  //grid->np   = npatch;
  grid->kind = alloc_info->grid_kind;
  
  alloc_info->grid = grid;
}

/* given the name of the checkpoint file, and the checkpoint file. 
// -> return value: a parameter read from the given checpoint file */
Parameter_T *bbn_parameter_query_from_checkpoint_file(const char *const par_name,FILE *const file)
{
  Parameter_T *par = 0;
  char line[MAX_ARR] = {'\0'};
  char *match_str = 0;
  unsigned i,npar;
  int found;
  
  fseek(file,0,SEEK_SET);
  
  /* check the header */
  fscanf(file,"%s",line);
  if (strcmp(line,ALLOC_HEADER))
      Error0("No header found. Checkpoint file got a problem.\n");
    
  /* read allocations */
  while (strcmp(line,ALLOC_FOOTER))
  {
    fscanf(file,"%s",line);
    if (!strcmp(line,ALLOC_FOOTER))
      break;
      
    char *v = strstr(line,"=");/* v -> "=..." */
    if (!v)
      Error0("No value found. Checkpoint file got a problem.\n");
    v++;
    
    if (strstr(line,"number_of_parameters"))
    {
      npar = (unsigned)atoi(v);
    }
  }
    
  /* read the binary parts */
  fseek(file,ftell(file)+1,SEEK_SET);/* +1 since fscanf won't read \n */
  
  /* is the cursor matched? */
  FReadP_bin(match_str);
  if (strcmp(match_str,PARAM_HEADER))
    Error0("It could not find the parameter header.\n");
  _free(match_str);
  
  found = 0;
  /* start reading one by one */
  for (i = 0; i < npar; ++i)
  {
    Parameter_T *p = calloc(1,sizeof(*p));
    pointerEr(p);
    
    FReadP_bin(p->lv);
    FReadP_bin(p->rv);
    FReadP_bin(p->rv_ip);
    FReadV_bin(p->rv_double);
    FReadP_bin(p->rv_array);
    FReadV_bin(p->rv_n);
    FReadV_bin(p->iterative);
    FReadV_bin(p->double_flg);
    
    if (strcmp_i(p->lv,par_name))
    {
      par   = p;
      found = 1;
      break;
    }
    free_given_parameter(p);
  }
  
  if (!found)
  {
    free_given_parameter(par);
    
    /* is the cursor matched? */
    FReadP_bin(match_str);
    if (strcmp(match_str,PARAM_FOOTER))
      Error0("It could not find the parameter footer.\n");
    _free(match_str);
  }
 
  return par;
}

/* read parameters */
static void read_parameters(struct checkpoint_header *const alloc_info,FILE *const file)
{
  /* read parameter contents */
  printf("~> Reading parameters from checkpoint file ...\n");
  fflush(stdout);
  
  const unsigned npar = alloc_info->npar;
  char *match_str;
  unsigned i;

  /* is the cursor matched? */
  FReadP_bin(match_str);
  if (strcmp(match_str,PARAM_HEADER))
    Error0("It could not find the parameter header.\n");
  _free(match_str);
  
  /* start reading one by one */
  for (i = 0; i < npar; ++i)
  {
    Parameter_T *p = parameters_global[i];
    
    FReadP_bin(p->lv);
    FReadP_bin(p->rv);
    FReadP_bin(p->rv_ip);
    FReadV_bin(p->rv_double);
    FReadP_bin(p->rv_array);
    FReadV_bin(p->rv_n);
    FReadV_bin(p->iterative);
    FReadV_bin(p->double_flg);
    
    //test
    if (0)
    {
      printf("%s         = %s\n",p->lv,p->rv);
      printf("rv_ip      = %s\n",p->rv_ip);
      printf("rv_double  = %f\n",p->rv_double);
      printf("rv_array   = %p\n",(void *)p->rv_array);
      printf("rv_n       = %u\n",p->rv_n);
      printf("iterative  = %u\n",p->iterative);
      printf("double_flg = %u\n",p->double_flg);
    }
    //end
  }
  /* is the cursor matched? */
  FReadP_bin(match_str);
  if (strcmp(match_str,PARAM_FOOTER))
    Error0("It could not find the parameter footer.\n");
  _free(match_str);
  
  /* incorporate the modified parameter in the parameter file
  // into parameter data base */
  incorporate_modified_checkpoint_par();
  
  /* free modified_checkpoint_par */
  free_modified_checkpoint_par();
  
  /* set the following parameters to default value */
  Pseti("did_resolution_change?",1);
  Pseti("did_NS_surface_change?",1);
  Pseti("did_AH_surface_change?",1);
  Pseti("use_previous_data",0);
}

/* free modified_checkpoint_par */
static void free_modified_checkpoint_par(void)
{
  unsigned i;
  
  for (i = 0; i < n_modified_checkpoint_par; i++)
  {
    _free(modified_checkpoint_par[i]->lv);
    _free(modified_checkpoint_par[i]->rv);
    _free(modified_checkpoint_par[i]->rv_ip);
    _free(modified_checkpoint_par[i]->rv_array);
    free(modified_checkpoint_par[i]);
  }
  _free(modified_checkpoint_par);
  modified_checkpoint_par   = 0;
  n_modified_checkpoint_par = 0;
}

/* incorporate the modified parameter in the parameter file
// into parameter data base */
static void incorporate_modified_checkpoint_par(void)
{
  unsigned np,i,n_found;
  
  /* free par "total_iterations_ip" since the new par file 
  // might have more iterations */
  free_parameter("total_iterations_ip");
  
  /* find the modified pars and save them */
  np      = 0;
  n_found = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
  {
    /* find */
    if (n_found < n_modified_checkpoint_par)
    {
      for (i = 0; i < n_modified_checkpoint_par; i++)
      {
        if (strcmp_i(parameters_global[np]->lv,modified_checkpoint_par[i]->lv))
        {
          printf("-> Modified parameter from checkpoint = %s\n",parameters_global[np]->lv);
          
          /* we must not have array type */
          assert(!parameters_global[np]->rv_array);
          
          _free(parameters_global[np]->rv);
          _free(parameters_global[np]->rv_ip);
          
          parameters_global[np]->rv         = modified_checkpoint_par[i]->rv;
          modified_checkpoint_par[i]->rv    = 0;
          parameters_global[np]->rv_double  = 0;
          parameters_global[np]->double_flg = 0;/* it's important to put this to 0 */
          parameters_global[np]->rv_ip      = modified_checkpoint_par[i]->rv_ip;
          modified_checkpoint_par[i]->rv_ip = 0;
          parameters_global[np]->iterative  = modified_checkpoint_par[i]->iterative;
            
          n_found++;
          break;
        }
      }
    }
    else
      break;
        
    np++;
  }
  
  /* update total iterations */
  total_iterations_ip();
}
 
/* read fields from the checkpoint file */
static void read_fields(struct checkpoint_header *const alloc_info,FILE *const file)
{  
  printf("~> Reading fields from checkpoint file ...\n");
  fflush(stdout);
  
  Grid_T *const grid = alloc_info->grid;
  char *match_str;
  unsigned p;
  
  /* is the cursor matched? */
  FReadP_bin(match_str);
  if (strcmp(match_str,FIELD_HEADER))
    Error0("It could not find the field header.\n");
  _free(match_str);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned f,count_nfld;
    
    if (IsItInsideBHPatch(patch))
      continue;
    
    FReadV_bin(count_nfld);
    assert(count_nfld);
    
    for (f = 0; f < count_nfld; ++f)
    {
      Field_T *field = 0;
      char *name = 0;
      double *v  = 0;
      
      FReadP_bin(name);
      FReadP_bin(v);
      //FReadP_bin(attr);
      
      field = patch->pool[Ind(name)];
      _free(name);
      field->v = v;
      v = 0;
    }
  }

  /* is the cursor matched? */
  FReadP_bin(match_str);
  if (strcmp(match_str,FIELD_FOOTER))
    Error0("It could not find the field footer.\n");
  _free(match_str);
}
