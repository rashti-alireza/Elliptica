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
  
  FILE *checkpoint_file = 0;
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
  checkpoint_file = open_binary_file("checkpoint.BIN")
  
  /* write all parameters in the checkpoint file */
  write_parameters(grid,checkpoint_file);
  
  /* write all fields value in the checkpoint file */
  write_fields(grid,checkpoint_file);
  
  /* close and replace checkpoint file */
  close_binary_file(checkpoint_file,"checkpoint.BIN");
  
  printf("} Writing checkpoint ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* open checkpoint file with the given name fname */
static void *open_binary_file(const char *const fname)
{
  FILE *file = 0;
  const char *const folder = Pgets("iteration_output");
  char file_path[100];
  
  sprintf(file_path,"%s/%s_temp",folder,fname);
  if (access(file_path,F_OK) != -1)/* if file exists */
    abortEr("File already exists.\n");
    
  file = fopen(file_path,"wb");
  pointerEr(file);
  
  return file;
}

/* write all of the pertinent parameters in the checkpoint file */
static void write_parameters(const Grid_T *const grid,FILE *const file)
{
  printf ("~> Writing all parameters in the checkpoint file ...\n");
}

/* write all of the fields in the checkpoint file */
static void write_fields(const Grid_T *const grid,FILE *const file)
{
  printf ("~> Writing all fields in the checkpoint file ...\n");
}

