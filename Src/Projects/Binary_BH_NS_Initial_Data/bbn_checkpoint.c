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
  unsigned np,p;

  sprintf(file_path,"%s/%s_temp",folder,checkpoint_file_name);
  if (access(file_path,F_OK) != -1)/* if file exists */
    abortEr("File already exists.\n");
    
  file = fopen(file_path,"w");
  pointerEr(file);
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
    
  /* NO white spaces since I'll use fscanf */
  fprintf(file,"number_of_parameters=%u\n",np);
  fprintf(file,"number_of_patches=%u\n",grid->np);
  fprintf(file,"grid_number=%u\n",grid->gn);
  fprintf(file,"grid_kind=%s\n",grid->kind);
  fprintf(file,HEADER_DONE);
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
  unsigned i,np;

  sprintf(file_path,"%s/%s_temp",folder,checkpoint_file_name);
  file = fopen(file_path,"ab");
  pointerEr(file);
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
    
  /* NOTE the order is crucial for reading part */
  Write(&"#parameters#",strlen("#parameters#")+1);
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
  unsigned p;
  
  sprintf(file_path,"%s/%s_temp",folder,checkpoint_file_name);
  file = fopen(file_path,"ab");
  pointerEr(file);
  
  /* NOTE the order is crucial for reading part */
  Write(&"#fields#",strlen("#fields#")+1);
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned f;
    Write(patch->name,strlen(patch->name)+1);
    Write(&patch->nfld,sizeof(patch->nfld));
    
    for (f = 0; f < patch->nfld; ++f)
    {
      Field_T *field = patch->pool[f];
      Write(field->name,strlen(field->name)+1);
      Write(field->v,nn);
      Write(field->attr,strlen(field->attr)+1);
    }
  }
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
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  char line[MAX_ARR] = {'\0'};
  char *grid_kind = 0;
  unsigned npatch,grid_number,npar,i;
  long cursor;
  
  sprintf(file_path,"%s/%s",folder,checkpoint_file_name);
  if (!(access(file_path,F_OK) != -1))/* if file dosn't exist */
    abortEr("File does not exist.\n");
  
  file = fopen(file_path,"r");
  pointerEr(file);
  
  /* read headers */
  while (strcmp(line,HEADER_DONE))
  {
    fscanf(file,"%s",line);
    if (!strcmp(line,HEADER_DONE))
      break;
      
    char *v = strstr(line,"=");/* v -> "=..." */
    if (!v)
      abortEr("No value found. Checkpoint file got a problem.\n");
    v++;
    
    if (strstr(line,"number_of_parameters"))
    {
      npar = (unsigned)atoi(v);
    }
    else if (strstr(line,"number_of_patches"))
    {
      npatch = (unsigned)atoi(v);
    }
    else if (strstr(line,"grid_number"))
    {
      grid_number = (unsigned)atoi(v);
    }
    else if (strstr(line,"grid_kind"))
    {
      grid_kind = dup_s(v);
    }
    else
      abortEr(NO_OPTION);
  }
  cursor = ftell(file);
  fclose(file);
  
  printf("~> Allocating parameters and patches ...\n");
  
  /* allocate parameters */
  free_parameter_db();
  parameters_global = calloc(npar+1,sizeof(*parameters_global));
  pointerEr(parameters_global);
  for (i = 0; i < npar; ++i)
    parameters_global = alloc_parameter(&parameters_global);
  
  /* allocate grid */
  grid       = alloc_grid();
  grid->gn   = grid_number;
  grid->np   = npatch;
  grid->kind = grid_kind;
  
  /* allocate patches */
  grid->patch = calloc(npatch+1,sizeof(*grid->patch));
  pointerEr(grid->patch);
  for (i = 0; i < npatch; ++i)
  {
    grid->patch[i] = calloc(1,sizeof(*grid->patch[i]));
    pointerEr(grid->patch[i]);
  }
  
  printf("} Reading checkpoint file ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
  return grid;
}
