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
  printf("{ Writing checkpoint ...\n");
  
  const char *const checkpoint_file_name = "checkpoint.dat";
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
  write_header(grid,checkpoint_file_name);
  
  /* write all parameters in the checkpoint file */
  write_parameters(grid,checkpoint_file_name);
  
  /* write all fields value in the checkpoint file */
  write_fields(grid,checkpoint_file_name);
  
  /* replace checkpoint file with the previous */
  move_checkpoint_file(checkpoint_file_name);
  
  printf("} Writing checkpoint ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* write header of checkpoint file */
static void write_header(const Grid_T *const grid,const char *const file_name)
{
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  unsigned np,p;

  sprintf(file_path,"%s/%s_temp",folder,file_name);
  if (access(file_path,F_OK) != -1)/* if file exists */
    abortEr("File already exists.\n");
    
  file = fopen(file_path,"w");
  pointerEr(file);
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
  fprintf(file,"number of parameters = %u\n",np);
  fprintf(file,"number of patches = %u\n",grid->np);
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    fprintf(file,"patch = %s\n",patch->name);
    fprintf(file,"number of field = %u\n",patch->nfld);
  }
  
  fclose(file);
}

/* replace checkpoint file with the previous one */
static void move_checkpoint_file(const char *const file_name)
{
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  char command[2*MAX_ARR];
  
  sprintf(file_path,"%s/%s",folder,file_name);
  sprintf(command,"mv %s/%s_temp %s/%s",
          folder,file_name,folder,file_name);
  printf("shell command:\n$ mv %s_temp %s\n",file_name,file_name);
  fflush(stdout);
  system(command);
}

/* write all of the pertinent parameters in the checkpoint file */
static void write_parameters(const Grid_T *const grid,const char *const file_name)
{
  printf ("~> Writing all parameters in the checkpoint file ...\n");
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  unsigned i,np;

  sprintf(file_path,"%s/%s_temp",folder,file_name);
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
static void write_fields(const Grid_T *const grid,const char *const file_name)
{
  printf("~> Writing all fields in the checkpoint file ...\n");
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  unsigned p;
  
  sprintf(file_path,"%s/%s_temp",folder,file_name);
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

