/*
// Alireza Rashti
// May 2018
*/

#include "global_variables.h"

/* initiate global variables */
int global_variables_init(char *const path)
{
  global_parameter = 0;
  make_global_path();
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
  
  if (last == 0)
  {
    abortEr("The name of input file must have extension.\n");
  }
  
  i = 0;
  for (p = path; p != last; p++)
  {
    if (*p == '.' || *p == '/')
      continue;
      
    name[i] = *p;
    i++;
  }
  
  global_inputfile_name = strdup(name);
  
  //TEST_START
    //printf("globale_inputfile_name = %s\n",global_inputfile_name);
  //end
  
}

/* making global_path */
static void make_global_path(void)
{
  char dir[MAX_ARR] = {'\0'};
  char *p;
  
  p = getcwd(dir,sizeof(dir));
  pointerEr(p);
  
  global_path = strdup(dir);
  
  //TEST_START
    //printf("globale_path = %s\n",global_path);
  //end
 
}

