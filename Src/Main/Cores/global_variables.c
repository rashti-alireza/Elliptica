/*
// Alireza Rashti
// May 2018
*/

#include "global_variables.h"

/* initiate global variables */
int global_variables_init(char *const path)
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
static void find_inputfile_name(char *const path)
{
  char *last,*p;
  char name[MAX_ARR] = {'\0'};
  int i;
  
  last = strrchr(path,'.');
  
  if (strstr(path,".in") == 0)
  {
    abortEr("The name of input file must have extension\".in\".\n");
  }
  
  i = 0;
  for (p = path; p != last; p++)
  {
    if (*p == '.' || *p == '/')
      continue;
      
    name[i] = *p;
    i++;
  }
  
  inputfile_name_global = strdup(name);
  
  //TEST_START
    //printf("globale_inputfile_name = %s\n",global_inputfile_name);
  //end
  
}

/* making global_path */
static void make_path_global(void)
{
  char dir[MAX_ARR] = {'\0'};
  char *p;
  
  p = getcwd(dir,sizeof(dir));
  pointerEr(p);
  
  path_global = strdup(dir);
  
  //TEST_START
    //printf("globale_path = %s\n",global_path);
  //end
 
}

