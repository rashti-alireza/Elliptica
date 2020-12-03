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
  find_relative_root_path(path);
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
  FILE *f_par = Fopen(path,"r");
  if (!f_par)
    Error0("The parameter file does not exist.\n");
  Fclose(f_par);
  
  last = strrchr(path,'.');
  if (!last)
    Error0("The name of the input file must have some extension.\n");
  
  p = strrchr(path,'/');
  if (p) p++;
  else   p = path;
  
  Psets("parameter_file_name",p);
  
  i = 0;
  while(p && p != last)
  {
    name[i] = *p;
    i++;
    p++;
  }
  
  Psets("parameter_file_name_stem",name);
}

/* finding the relative root path i.e. 
// where the location of input file is. */
static void find_relative_root_path(const char *const path)
{
  char dir[MAX_ARR] = {'\0'};
  char *p;
  
  sprintf(dir,"%s",path);
  p = strrchr(dir,'/');
  
  if (p) 
    p[0] = '\0';
  else
    sprintf(dir,".");
  
  Psets("relative_root_path",dir);
}

