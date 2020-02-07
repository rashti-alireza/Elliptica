/*
// Alireza Rashti
// May 2018
*/

#include "global_variables.h"

/* initiate global variables */
int init_global_variables(const char *const path)
{
  initial_time_global = time(0);
  grids_global        = 0;
  parameters_global   = 0;
  projects_global     = 0;
  make_path_global();
  find_inputfile_name(path);
  
  return EXIT_SUCCESS;
}

/* finding inputfile name */
static void find_inputfile_name(const char *const path)
{
  char *last;
  const char *p;
  char name[MAX_ARR] = {'\0'};
  int i;
  
  /* does the file exist */
  FILE *f_par = fopen(path,"r");
  if (!f_par)
    abortEr("The parameter file does not exist.\n");
  fclose(f_par);
  
  last = strrchr(path,'.');
  if (!last)
    abortEr("The name of the input file must have some extension.\n");
  
  i = 0;
  for (p = path; p != last; p++)
  {
    if (*p == '.' || *p == '/')
      continue;
      
    name[i] = *p;
    i++;
  }
  
  inputfile_name_global = dup_s(name);
  
  /*TEST_START
    //printf("globale_inputfile_name = %s\n",global_inputfile_name);
  end */
  
}

/* making global_path which shows the location of input file */
static void make_path_global(void)
{
  char dir[MAX_ARR] = {'\0'};
  char *p;
  
  p = getcwd(dir,sizeof(dir));
  pointerEr(p);
  
  path_global = dup_s(dir);
  
  /*TEST_START
    //printf("globale_path = %s\n",global_path);
  end */
 
}

