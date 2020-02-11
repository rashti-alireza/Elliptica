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
  if (dt+last_checkpoint_was < now)
  {
    printf("~> It's early for writing checkpoint.\n");
    printf("} Writing checkpoint ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
    return;
  }
  
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
  
  /* open checkpoint file */
  //checkpoint_file = open_binary_file("checkpoint.BIN")
  
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

/* replace checkpoint file with the previous one */
static void move_checkpoint_file(const char *const file_name)
{
  return;
  
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  char command[2*MAX_ARR];
  
  sprintf(file_path,"%s/%s",folder,file_name);
  sprintf(command,"mv %s/%s_temp %s/%s",
          folder,file_name,folder,file_name);
  printf("shell command:\n$ mv %s_temp %s\n\n",file_name,file_name);
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
  
  sprintf(file_path,"%s/%s_temp",folder,file_name);
  if (access(file_path,F_OK) != -1)/* if file exists */
    abortEr("File already exists.\n");
    
  file = fopen(file_path,"wb");
  pointerEr(file);
  
  /* write parameters */
  unsigned np;

  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
  {
    Parameter_T *par = parameters_global[np];
    fwrite(par,sizeof(*par),1,file);
  }
  
  fclose(file);
  
  UNUSED(grid);
}

/* write all of the fields in the checkpoint file */
static void write_fields(const Grid_T *const grid,const char *const file_name)
{
  return;
  UNUSED(grid);
  
  printf ("~> Writing all fields in the checkpoint file ...\n");
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[MAX_ARR];
  
  sprintf(file_path,"%s/%s_temp",folder,file_name);
  if (access(file_path,F_OK) != -1)/* if file exists */
  {
    file = fopen(file_path,"w");
    pointerEr(file);
  }
  else
    abortEr("File does not already exist.\n");
    
  fclose(file);  
}

