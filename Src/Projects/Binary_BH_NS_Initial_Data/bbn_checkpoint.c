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
  /* no white spaces */
  fprintf(file,"number_of_parameters=%u\n",np);
  fprintf(file,"number_of_patches=%u\n",grid->np);
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    fprintf(file,"patch=%s\n",patch->name);
    fprintf(file,"number_of_field=%u\n",patch->nfld);
  }
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
  unsigned p;
  
  sprintf(file_path,"%s/%s",folder,checkpoint_file_name);
  if (!(access(file_path,F_OK) != -1))/* if file dosn't exist */
    abortEr("File does not exist.\n");
  
  file = fopen(file_path,"r");
  pointerEr(file);
  
  /* read headers */
  while (!strcmp(line,HEADER_DONE))
  {
    fscanf(file,"%[^\n]",line);
    printf("%s",line);
  }
  
  printf("} Reading checkpoint file ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
  return grid;
}
