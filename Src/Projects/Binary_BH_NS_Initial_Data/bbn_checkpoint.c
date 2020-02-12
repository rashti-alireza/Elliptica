/*
// Alireza Rashti
// February 2020
*/

#include "bbn_checkpoint.h"

/* write checkpoint for the given grid */
void bbn_write_checkpoint(const Grid_T *const grid)
{
  /* print some descriptions */
  pr_line_custom('=');
  printf("{ Writing checkpoint file ...\n");
  
  const double dt  = Pgetd("checkpoint_dt_hours");
  const double now = get_time_sec()/(3600);
  static double last_checkpoint_was = 0;/* the time where the last 
                                        // checkpoint happened in hours */
  
  /* some checks */
  /*if (dt+last_checkpoint_was < now)
  {
    printf("~> It's early for writing checkpoint.\n");
    printf("} Writing checkpoint ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
    return;
  }*/
  
  /* some checks */
  if (!grid)
  {
    printf("~> The given grid is empty.\n");
    printf("} Writing checkpoint ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
    return;
  }
  
  last_checkpoint_was = now;
  
  /* NOTE: the order of writing is crucial */
  
  /* write header */
  write_header(grid);
  
  /* write all parameters in the checkpoint file */
  write_parameters(grid);
  
  /* write all fields value in the checkpoint file */
  write_fields(grid);
  
  /* replace checkpoint file with the previous */
  move_checkpoint_file();
  
  printf("} Writing checkpoint file ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* write header of checkpoint file */
static void write_header(const Grid_T *const grid)
{
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  unsigned np;

  sprintf(file_path,"%s/%s_temp",folder,checkpoint_file_name);
  if (access(file_path,F_OK) != -1)/* if file exists */
    abortEr("File already exists.\n");
    
  file = fopen(file_path,"w");
  pointerEr(file);
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
    
  /* NO white spaces since I'll use fscanf */
  fprintf(file,"%s\n",ALLOC_HEADER);
  fprintf(file,"number_of_parameters=%u\n",np);
  fprintf(file,"number_of_patches=%u\n",grid->np);
  fprintf(file,"grid_number=%u\n",grid->gn);
  fprintf(file,"grid_kind=%s\n",grid->kind);
  fprintf(file,"%s\n",ALLOC_FOOTER);
  fclose(file);
}

/* replace checkpoint file with the previous one */
static void move_checkpoint_file(void)
{
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  char command[2*MAX_ARR];
  
  sprintf(file_path,"%s/%s",folder,checkpoint_file_name);
  sprintf(command,"mv %s/%s_temp %s/%s",
          folder,checkpoint_file_name,folder,checkpoint_file_name);
  printf("shell command:\n$ mv %s_temp %s\n",checkpoint_file_name,checkpoint_file_name);
  fflush(stdout);
  system(command);
}

/* write all of the pertinent parameters in the checkpoint file */
static void write_parameters(const Grid_T *const grid)
{
  printf ("~> Writing all parameters in the checkpoint file ...\n");
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  char title_line[MAX_ARR];
  char *const p_title_line = title_line;/* defined to avoid warning */
  unsigned i,np;

  sprintf(file_path,"%s/%s_temp",folder,checkpoint_file_name);
  file = fopen(file_path,"ab");
  pointerEr(file);
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
    
  /* NOTE the order is crucial for reading part */
  sprintf(title_line,"%s",PARAM_HEADER);
  Write(p_title_line,strlen(title_line)+1);
  for (i = 0; i < np; ++i)
  {
    Parameter_T *p = parameters_global[i];
    
    Write(p->lv,strlen(p->lv)+1);
    Write(p->rv,strlen(p->rv)+1);
    Write(p->rv_ip,strlen(p->rv_ip)+1);
    Write(&p->rv_double,1);
    Write(p->rv_array,p->rv_n);
    Write(&p->rv_n,1);
    Write(&p->iterative,1);
    Write(&p->double_flg,1);
  }
  sprintf(title_line,"%s",PARAM_FOOTER);
  Write(p_title_line,strlen(title_line)+1);
  fclose(file);
  
  UNUSED(grid);
}

/* write all of the fields in the checkpoint file */
static void write_fields(const Grid_T *const grid)
{
  printf("~> Writing all fields in the checkpoint file ...\n");
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  char title_line[MAX_ARR];
  char *const p_title_line = title_line;/* defined to avoid warning */
  unsigned p;
  
  sprintf(file_path,"%s/%s_temp",folder,checkpoint_file_name);
  file = fopen(file_path,"ab");
  pointerEr(file);
  
  /* NOTE the order is crucial for reading part */
  
  sprintf(title_line,"%s",FIELD_HEADER);
  Write(p_title_line,strlen(title_line)+1);
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned f;
    Write(patch->name,strlen(patch->name)+1);
    Write(&patch->nn,sizeof(patch->nn));
    Write(&patch->nfld,sizeof(patch->nfld));
    
    for (f = 0; f < patch->nfld; ++f)
    {
      Field_T *field = patch->pool[f];
      Write(field->name,strlen(field->name)+1);
      Write(field->v,nn);
      Write(field->attr,strlen(field->attr)+1);
    }
  }
  sprintf(title_line,"%s",FIELD_FOOTER);
  Write(p_title_line,strlen(title_line)+1);
  fclose(file);  
}

/* read checkpoint file and creat grid and parameters accordingly */
Grid_T *bbn_read_checkpoint(void)
{
  /* print some descriptions */
  pr_line_custom('=');
  printf("{ Reading checkpoint file ...\n");
  
  Grid_T *grid = 0;
  FILE *file = 0;
  struct checkpoint_header alloc_info[1] = {0};
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  
  sprintf(file_path,"%s/%s",folder,checkpoint_file_name);
  if (!(access(file_path,F_OK) != -1))/* if file dosn't exist */
    abortEr("File does not exist.\n");
  
  file = fopen(file_path,"r");
  pointerEr(file);
  
  /* reading the header for allocations */
  read_header(alloc_info,file);
  
  /* alloc grid, patches and parameters */
  alloc_db(alloc_info);
 
  /* read parameters */
  read_parameters(alloc_info,file);
  
  /* read fields content */
  read_fields(alloc_info,file);
  
  fclose(file);

  printf("} Reading checkpoint file ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
  return grid;
}

/* reading the header for allocations */
static void read_header(struct checkpoint_header *const alloc_info,FILE *const file)
{  
  printf("~> Reading checkpoint file header ...\n");
  
  char line[MAX_ARR] = {'\0'};

  /* check the header */
  fscanf(file,"%s",line);
  if (strcmp(line,ALLOC_HEADER))
      abortEr("No header found. Checkpoint file got a problem.\n");
    
  /* read allocations */
  while (strcmp(line,ALLOC_FOOTER))
  {
    fscanf(file,"%s",line);
    if (!strcmp(line,ALLOC_FOOTER))
      break;
      
    char *v = strstr(line,"=");/* v -> "=..." */
    if (!v)
      abortEr("No value found. Checkpoint file got a problem.\n");
    v++;
    
    if (strstr(line,"number_of_parameters"))
    {
      alloc_info->npar = (unsigned)atoi(v);
    }
    else if (strstr(line,"number_of_patches"))
    {
      alloc_info->npatch = (unsigned)atoi(v);
    }
    else if (strstr(line,"grid_number"))
    {
      alloc_info->grid_number = (unsigned)atoi(v);
    }
    else if (strstr(line,"grid_kind"))
    {
      alloc_info->grid_kind = dup_s(v);
    }
    else
      abortEr(NO_OPTION);
  }
  /* read the binary parts */
  fseek(file,ftell(file)+1,SEEK_SET);/* +1 since fscanf won't read \n */
}

/* alloc grid, patches and parameters */
static void alloc_db(struct checkpoint_header *const alloc_info)
{
  printf("~> Allocating parameters and patches ...\n");
  
  const unsigned npatch = alloc_info->npatch,
            grid_number = alloc_info->grid_number,
            npar        = alloc_info->npar;
  Grid_T *grid = 0;
  unsigned i;
  
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
  grids_global  = 0;
  grid       = alloc_grid();
  grid->gn   = grid_number;
  grid->np   = npatch;
  grid->kind = alloc_info->grid_kind;
  
  /* allocate patches */
  grid->patch = calloc(npatch+1,sizeof(*grid->patch));
  pointerEr(grid->patch);
  for (i = 0; i < npatch; ++i)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    pointerEr(grid->patch[i]);
  }
 
  alloc_info->grid = grid;
}

/* read parameters */
static void read_parameters(struct checkpoint_header *const alloc_info,FILE *const file)
{
  /* read parameter contents */
  printf("~> Reading parameters content ...\n");
  
  const unsigned npar = alloc_info->npar;
  char *match_str;
  unsigned match_str_s,i;

  /* is the cursor matched? */
  ReadP(match_str,match_str_s);
  if (strcmp(match_str,PARAM_HEADER))
    abortEr("It could not find the parameter header.\n");
  _free(match_str);
  
  /* start reading one by one */
  for (i = 0; i < npar; ++i)
  {
    Parameter_T *p = parameters_global[i];
    unsigned mem_size;
    
    ReadP(p->lv,          mem_size);
    ReadP(p->rv,          mem_size);
    ReadP(p->rv_ip,       mem_size);
    ReadV(&p->rv_double,  mem_size);
    ReadP(p->rv_array,    mem_size);
    ReadV(&p->rv_n,       mem_size);
    ReadV(&p->iterative,  mem_size);
    ReadV(&p->double_flg, mem_size);
    
    //test
    printf("%s         = %s\n",p->lv,p->rv);
    printf("rv_ip      = %s\n",p->rv_ip);
    printf("rv_double  = %f\n",p->rv_double);
    printf("rv_array   = %p\n",(void *)p->rv_array);
    printf("rv_n       = %u\n",p->rv_n);
    printf("iterative  = %u\n",p->iterative);
    printf("double_flg = %u\n\n\n",p->double_flg);
    //end
  }
  /* is the cursor matched? */
  ReadP(match_str,match_str_s);
  if (strcmp(match_str,PARAM_FOOTER))
    abortEr("It could not find the parameter footer.\n");
  _free(match_str);
}
 
/* read fields content */
static void read_fields(struct checkpoint_header *const alloc_info,FILE *const file)
{  
  printf("~> Reading fields content ...\n");
  Grid_T *const grid = alloc_info->grid;
  char *match_str;
  unsigned match_str_s;
  unsigned p;
  
  /* is the cursor matched? */
  ReadP(match_str,match_str_s);
  if (strcmp(match_str,FIELD_HEADER))
    abortEr("It could not find the field header.\n");
  _free(match_str);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned f;
    unsigned mem_size;
    ReadP(patch->name,mem_size);
    ReadV(&patch->nn,mem_size);
    ReadV(&patch->nfld,mem_size);
    
    /* allocate fields */
    assert(patch->nfld);
    patch->pool = calloc(patch->nfld,sizeof(*patch->pool));
    pointerEr(patch->pool);
    
    for (f = 0; f < patch->nfld; ++f)
    {
      patch->pool[f] = calloc(1,sizeof(*patch->pool[f]));
      pointerEr(patch->pool[f]);
      
      Field_T *field = patch->pool[f];
      ReadP(field->name,mem_size);
      ReadP(field->v,mem_size);
      ReadP(field->attr,mem_size);
      //test
      printf("field = %s, attr = %s\n",field->name,field->attr);
      //end
    }
  }

  /* is the cursor matched? */
  ReadP(match_str,match_str_s);
  if (strcmp(match_str,FIELD_FOOTER))
    abortEr("It could not find the field footer.\n");
  _free(match_str);
}
 