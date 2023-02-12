/*
// Alireza Rashti
// February 2020
*/

#include "checkpoint.h"


/* if you set a parameter as below in the par file it uses this parameter
// and disregards its value obtained from the checkpoint and if
// the parameter does not exist it is added:
// e.g.
// `CHECKPOINT_SET_PARAM_` n_a = 4(x6)
// `CHECKPOINT_SET_PARAM_` Solving_Max_Number_of_Iteration = 0 
// 
// if you wanna change the output directories after loading from
// a checkpoint file change the following in parameter files
// and note that the folders must be made already:
// `CHECKPOINT_SET_PARAM_`      top_directory = path1
// `CHECKPOINT_SET_PARAM_` `P_` my_directory  = path2
// `CHECKPOINT_SET_PARAM_` `P_` Diagnostics   = path3
//
// where `CHECKPOINT_SET_PARAM_` is the macro to instruct the 
// checkpoint param to change and `P_` is the prefix of the project. */


static Uint n_modified_checkpoint_par;/* number of modified checkpoint par */
static Parameter_T **modified_checkpoint_par;/* modified pars in par file
                                           // to be used after loading of
                                           // the checkpoint file. */

/* write checkpoint for the given phys->grid.
// NOTE: the order of writing and reading is crucial */
void write_checkpoint(Physics_T *const phys,const char *const out_dir)
{
  FUNC_TIC
  
  /* if no checkpoint at all requested */
  if (Pcmps("checkpoint_every","never"))
  {
    printf(Pretty0"No checkpoint requested.\n");
    FUNC_TOC
    return;
  }
  if (!phys)
  {
    printf(Pretty0"The given physics is empty.\n");
    FUNC_TOC
    return;
  }
  
  Grid_T *const grid = phys->grid;
  if (!grid)
  {
    printf(Pretty0"The given grid is empty.\n");
    FUNC_TOC
    return;
  }
  
  FILE *file = 0;
  char file_path[MAX_ARR];
  char msg[MAX_ARR];
  char *const p_msg = msg;/* to avoid GCC warning for FWriteP_bin */
  const double dt  = Pgetd("checkpoint_every");/* unit is hour */
  const double now = get_time_sec()/(3600.);
  static double last_checkpoint_was = 0;/* the time where the last 
                                        // checkpoint happened in hours */
  
  if (dt+last_checkpoint_was > now)
  {
    printf(Pretty0"It's early to write checkpoint.\n");
    FUNC_TOC
    return;
  }
  
  last_checkpoint_was = now;
  
  /* NOTE: the order of writing is super crucial for reading */
  
  /* write header */
  write_header(grid,out_dir);
  
  /* write all parameters in the checkpoint file */
  write_parameters(grid,out_dir);
  
  /* write all fields value in the checkpoint file */
  write_fields(grid,out_dir);
  
  /* successful message at the end of the checkpoint file */
  sprintf(file_path,"%s/%s_temp",out_dir,CHECKPOINT_FILE_NAME);
  file = Fopen(file_path,"a");
  sprintf(msg,"%s",END_MSG);
  FWriteP_bin(p_msg,strlen(msg)+1);
  Fclose(file);
  
  /* replace checkpoint file with the previous */
  move_checkpoint_file(out_dir);
  
  FUNC_TOC
}

/* ->return value: if the chekpoint is completed 1, otherwise 0. */
int is_checkpoint_sound(const char *const file_path)
{
  int ret = 0;
  FILE *file;
  char msg[MAX_ARR];
  int msg_len = (int)strlen(END_MSG)+1;
  
  file = Fopen(file_path,"r");
  
  fseek(file,-msg_len,SEEK_END);
  assert(fread(msg,(Uint)msg_len,1,file));
  
  if (strstr(msg,END_MSG))
    ret = 1;
  else
    ret = 0;
      
  Fclose(file);
  
  return ret;
}

/* write header of checkpoint file */
static void write_header(const Grid_T *const grid,const char *const folder)
{
  FILE *file = 0;
  char file_path[MAX_ARR];
  Uint np;

  sprintf(file_path,"%s/%s_temp",folder,CHECKPOINT_FILE_NAME);
  if (!access(file_path,F_OK))/* if file exists */
    Error0("File already exists.\n");
  
    
  file = Fopen(file_path,"w");
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
    
  /* NO white spaces since I'll use fscanf */
  fprintf(file,"%s\n",ALLOC_HEADER);
  fprintf(file,"number_of_parameters=%u\n",np);
  //fprintf(file,"number_of_patches=%u\n",grid->np);
  fprintf(file,"grid_number=%u\n",grid->gn);
  fprintf(file,"grid_kind=%s\n",Pgets("grid_kind"));
  fprintf(file,"%s\n",ALLOC_FOOTER);
  Fclose(file);
}

/* replace checkpoint file with the previous one */
static void move_checkpoint_file(const char *const folder)
{
  char command[2*MAX_ARR];
  
  sprintf(command,"mv %s/%s_temp %s/%s",
          folder,CHECKPOINT_FILE_NAME,folder,CHECKPOINT_FILE_NAME);
  shell_command(command);
  
  printf(Pretty0"checkpoint file path:\n%s/%s\n",
                 folder,CHECKPOINT_FILE_NAME);
  fflush(stdout);
}

/* write all of the pertinent parameters in the checkpoint file */
static void write_parameters(const Grid_T *const grid,const char *const folder)
{
  printf (Pretty0"Writing parameters in checkpoint file ...\n");
  fflush(stdout);
  
  FILE *file = 0;
  char file_path[MAX_ARR];
  char title_line[MAX_ARR] = {'\0'};
  char *const p_title_line = title_line;/* to avoid GCC warning for FWriteP_bin */
  Uint i,np;

  sprintf(file_path,"%s/%s_temp",folder,CHECKPOINT_FILE_NAME);
  file = Fopen(file_path,"ab");
  
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
  
  Fclose(file);
  
  UNUSED(grid);
}

/* write all of the fields in the checkpoint file */
static void write_fields(const Grid_T *const grid,const char *const folder)
{
  printf(Pretty0"Writing fields in checkpoint file ...\n");
  fflush(stdout);
  
  FILE *file = 0;
  const char *const list   = Pgets("checkpoint_save");
  char file_path[MAX_ARR];
  char title_line[MAX_ARR] = {'\0'};
  char regex[MAX_ARR]      = {'\0'};
  char *const p_title_line = title_line;/* to avoid GCC warning for FWriteP_bin */
  Uint p;
  
  sprintf(file_path,"%s/%s_temp",folder,CHECKPOINT_FILE_NAME);
  file = Fopen(file_path,"ab");
  
  /* NOTE the order is crucial for reading part */
  
  sprintf(title_line,"%s",FIELD_HEADER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Uint nn = patch->nn;
    Uint f,count_nfld;
    
    /* count number fields we want to save */
    count_nfld = 0;
    for (f = 0; f < patch->nfld; ++f)
    {
      Field_T *field = patch->fields[f];
      
      sprintf(regex,"\\b%s\\b",field->name);
      if (patch->fields[f]->v && regex_search(regex,list))
        count_nfld++;
    }
    FWriteV_bin(count_nfld,1);
    
    for (f = 0; f < patch->nfld; ++f)
    {
      Field_T *field = patch->fields[f];
      
      sprintf(regex,"\\b%s\\b",field->name);
      if (patch->fields[f]->v && regex_search(regex,list))
      {
        FWriteP_bin(field->name,strlen(field->name)+1);
        FWriteP_bin(field->v,nn);
      }
    }
    
  }
  
  sprintf(title_line,"%s",FIELD_FOOTER);
  FWriteP_bin(p_title_line,strlen(title_line)+1);
  
  Fclose(file);  
}

/* read checkpoint file and creat grid and parameters accordingly
// no patch andfields at this stage, patch and fields must be 
// added first by the project. */
static Grid_T *read_grid_and_params(FILE *const file)
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
  
  /* free */
  Free(alloc_info->grid_kind);
  
  return grid;
}
  
/* reading the header for allocations */
static void read_header(struct checkpoint_header *const alloc_info,FILE *const file)
{  
  printf(Pretty0"Reading checkpoint file header ...\n");
  fflush(stdout);
  
  char line[MAX_ARR] = {'\0'};
  int sys_ret;
  
  fseek(file,0,SEEK_SET);
  
  /* check the header */
  sys_ret = fscanf(file,"%s",line);
  if (strcmp(line,ALLOC_HEADER))
      Error0("No header found. Checkpoint file got a problem.\n");
    
  /* read allocations */
  while (strcmp(line,ALLOC_FOOTER))
  {
    sys_ret = fscanf(file,"%s",line);
    if (!strcmp(line,ALLOC_FOOTER))
      break;
      
    char *v = strstr(line,"=");/* v -> "=..." */
    if (!v)
      Error0("No value found. Checkpoint file got a problem.\n");
    v++;
    
    if (strstr(line,"number_of_parameters"))
    {
      alloc_info->npar = (Uint)atoi(v);
    }
    //else if (strstr(line,"number_of_patches"))
    //{
      //alloc_info->npatch = (Uint)atoi(v);
    //}
    else if (strstr(line,"grid_number"))
    {
      alloc_info->grid_number = (Uint)atoi(v);
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
  UNUSED(sys_ret)
}

/* find and save modified checkpoint pars specified at the par file
// also delete these pars from the par db */
static void find_and_save_modified_checkpoint_pars(void)
{
  modified_checkpoint_par   = 0;/* global var in this file */
  n_modified_checkpoint_par = 0;/* global var in this file */
  const char *const keyword_prefix = CHECKPOINT_SET_PARAM_ ;
  char str[MAX_ARR],*pstr;
  Uint np,nmpar;
  
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
      IsNull(modified_checkpoint_par);
      modified_checkpoint_par[nmpar] = parameters_global[np];
      parameters_global[np] = 0;
      parameters_global[np] = calloc(1,sizeof(*parameters_global[np]));
      IsNull(parameters_global[np]);
      
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
  printf(Pretty0"Allocating parameters and patches ...\n");
  fflush(stdout);
  
  const Uint grid_number = alloc_info->grid_number,
                 npar    = alloc_info->npar;
  Grid_T *grid = 0;
  Uint i;
  
  /* find and save modified checkpoint pars specified at the pars file */
  find_and_save_modified_checkpoint_pars();
  
  /* allocate parameters */
  free_parameter_db();
  parameters_global = calloc(npar+1,sizeof(*parameters_global));
  IsNull(parameters_global);
  for (i = 0; i < npar; ++i)
  {
    parameters_global[i] = calloc(1,sizeof(*parameters_global[i]));
    IsNull(parameters_global[i]);
  }
  
  /* allocate grid */
  free_grid_db();
  grid       = alloc_grid();
  grid->gn   = grid_number;
  //grid->np   = npatch;
  grid->kind = set_grid_kind(alloc_info->grid_kind);
  
  alloc_info->grid = grid;
}

/* given the name of the checkpoint file, and the checkpoint file. 
// -> return value: a parameter read from the given checpoint file, null if not found. */
Parameter_T *parameter_query_from_checkpoint(const char *const par_name,FILE *const file)
{
  Parameter_T *par = 0;
  char line[MAX_ARR] = {'\0'};
  char *match_str = 0;
  Uint i,npar = 0;
  int found;
  int sys_ret;
  
  fseek(file,0,SEEK_SET);
  
  /* check the header */
  sys_ret = fscanf(file,"%s",line);

  if (strcmp(line,ALLOC_HEADER))
      Error0("No header found. Checkpoint file got a problem.\n");
    
  /* read allocations */
  while (strcmp(line,ALLOC_FOOTER))
  {
    sys_ret = fscanf(file,"%s",line);
    if (!strcmp(line,ALLOC_FOOTER))
      break;
      
    char *v = strstr(line,"=");/* v -> "=..." */
    if (!v)
      Error0("No value found. Checkpoint file got a problem.\n");
    v++;
    
    if (strstr(line,"number_of_parameters"))
    {
      npar = (Uint)atoi(v);
    }
  }
    
  /* read the binary parts */
  fseek(file,ftell(file)+1,SEEK_SET);/* +1 since fscanf won't read \n */
  
  /* is the cursor matched? */
  FReadP_bin(match_str);
  if (strcmp(match_str,PARAM_HEADER))
    Error0("It could not find the parameter header.\n");
  Free(match_str);
  
  found = 0;
  /* start reading one by one */
  for (i = 0; i < npar; ++i)
  {
    Parameter_T *p = calloc(1,sizeof(*p));
    IsNull(p);
    
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
    Free(match_str);
  }
  
  UNUSED(sys_ret)
  return par;
}

/* read parameters */
static void read_parameters(struct checkpoint_header *const alloc_info,FILE *const file)
{
  /* read parameter contents */
  printf(Pretty0"Reading parameters from checkpoint file ...\n");
  fflush(stdout);
  
  const Uint npar = alloc_info->npar;
  char *match_str;
  Uint i;

  /* is the cursor matched? */
  FReadP_bin(match_str);
  if (strcmp(match_str,PARAM_HEADER))
    Error0("It could not find the parameter header.\n");
  Free(match_str);
  
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
  Free(match_str);
  
  /* incorporate the modified parameter in the parameter file
  // into parameter data base */
  incorporate_modified_checkpoint_par();
  
  /* free modified_checkpoint_par */
  free_modified_checkpoint_par();
}

/* free modified_checkpoint_par */
static void free_modified_checkpoint_par(void)
{
  Uint i;
  
  for (i = 0; i < n_modified_checkpoint_par; i++)
  {
    Free(modified_checkpoint_par[i]->lv);
    Free(modified_checkpoint_par[i]->rv);
    Free(modified_checkpoint_par[i]->rv_ip);
    Free(modified_checkpoint_par[i]->rv_array);
    free(modified_checkpoint_par[i]);
  }
  Free(modified_checkpoint_par);
  modified_checkpoint_par   = 0;
  n_modified_checkpoint_par = 0;
}

/* incorporate the modified parameter in the parameter file
// into parameter data base */
static void incorporate_modified_checkpoint_par(void)
{
  Uint np,i,n_found;
  
  /* free par "total_iterations_ip" since the new par file 
  // might have more iterations */
  free_parameter("total_iterations_ip");
  
  /* if the parameter already exists: */
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
          printf(Pretty0"Modified parameter from checkpoint = %s\n",parameters_global[np]->lv);
          
          /* we must not have array type */
          assert(!parameters_global[np]->rv_array);
          
          Free(parameters_global[np]->rv);
          Free(parameters_global[np]->rv_ip);
          
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
  /* if the parameter is new: */
  /* add the new parameter */
  for (i = 0; i < n_modified_checkpoint_par; i++)
  {
    /* if does not exist */
    if (!get_parameter(modified_checkpoint_par[i]->lv))
    {
      printf(Pretty0"Adding new parameter               = %s\n",modified_checkpoint_par[i]->lv);
          
      add_parameter(modified_checkpoint_par[i]->lv,modified_checkpoint_par[i]->rv);
    }
  }
  
  /* update total iterations */
  total_iterations_ip();
}
 
/* read fields from the checkpoint file */
void read_fields_from_checkpoint_file(Physics_T *const phys,FILE *const file)
{
  FUNC_TIC
  
  Grid_T *const grid = phys->grid;
  char *match_str;
  Uint p;
  
  /* is the cursor matched? */
  FReadP_bin(match_str);
  if (strcmp(match_str,FIELD_HEADER))
    Error0("It could not find the field header.\n");
  Free(match_str);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    Uint f,count_nfld;
    
    FReadV_bin(count_nfld);
    
    if(!count_nfld)
    {
      printf(Pretty0"No field has been found to read!\n");
    }
    
    for (f = 0; f < count_nfld; ++f)
    {
      Field_T *field = 0;
      char *name = 0;
      double *v  = 0;
      
      FReadP_bin(name);
      FReadP_bin(v);
      //FReadP_bin(attr);
      
      field = patch->fields[Ind(name)];
      Free(name);
      field->v = v;
      v = 0;
    }
  }

  /* is the cursor matched? */
  FReadP_bin(match_str);
  if (strcmp(match_str,FIELD_FOOTER))
    Error0("It could not find the field footer.\n");
  Free(match_str);

  FUNC_TOC
}


/* ->: checkpoint file whose position ready for fields to be read.
// loading the grid and parameters from checkpoint file.
// No check at this stage, we assume every thing is consistence.
// NOTE: this only load the grid and parameter without patches and fields.
// patches and fields must be added separately and then be read by function:
// read_fields_from_checkpoint_file using the returned file.
// the path of checkpoint file is saved in "checkpoint_file_path". */
void *open_checkpoint_file_then_read_grid_and_params(Physics_T *const phys)
{
  FUNC_TIC
  
  FILE *file   = 0;
  const char *const checkpoint_file_path = Pgets("checkpoint_file_path");
  
  if (access(checkpoint_file_path,F_OK))/* if file does not exist */
    Errors("Checkpoint file does not exist at\n%s\n",checkpoint_file_path);
    
  file = Fopen(checkpoint_file_path,"r");
  IsNull(file);
  
  phys->grid = read_grid_and_params(file);
  
  FUNC_TOC
  return file;
}

/* check if there is a consistence checkpoint file to be used 
// for initialization. it uses leading number of cur_out_dir
// to find out if in previus numbers there is a sound checkpoint file.
// -> return value : 1 if exists, 0 otherwise.
// also it sets par "checkpoint_file_path" to be used for loading. */
int can_we_use_checkpoint(const char *const cur_out_dir)
{
  int ret = 0;
  const int FOLDER_TYPE   = 4;
  const Uint MAX_LIST_NUM = 10;
  DIR *prev_dir;
  struct dirent *ent;
  struct stat st = {0};/* status of files */
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
  Uint count,i;
  
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
        if (is_checkpoint_sound(str))
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
  
  /* set the path of checkpoint file and some preparation */
  if (ret)
  {
    Psets("checkpoint_file_path",prev_data_file_path);
    
    printf(Pretty0"checkpoint file found at:\n%s\n",prev_data_file_path);
    
    /* remove the current directory */
    sprintf(str,"rm -rf %s",cur_out_dir);
    shell_command(str);
  }
  
  return ret;
}
