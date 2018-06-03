/*
// Alireza Rashti
// May 2018
*/

#include "global_variables.h"

/* initiate global variables */
int global_variables_init(char *path)
{
  global_parameter = 0;
  make_global_path(path);
  find_inputfile_name(path);
  
  return EXIT_SUCCESS;
}

/* finding inputfile name */
static void find_inputfile_name(char *path)
{
  char *first,*last,*p;
  char name[1000] = {'\0'};
  
  first = strrchr(path,'\\');
  pointerEr(first);
  
  last = strrchr(path,'.');
  pointerEr(last);
  
  if (last <= first)
  {
    fprintf(stderr,ERROR_MASSAGE"The name of input file must be like:"
    " \"input.in\"\n");
    printf(ERROR_MASSAGE_EXIT);
    abort();
    
  }
  
  for (p = first; *p != '.'; p++)
    name[p-first] = first[p-first];
  
  global_inputfile_name = strdup(name);
  
}

/* making global_path */
static void make_global_path(char *path)
{
  char *first,*last,*p;
  char dir[1000] = {'\0'};
  
  first = strchr(path,'\\');
  pointerEr(first);
  last = strrchr(path,'\\');
  pointerEr(last);
  
  for (p = first; p < last; p++)
    dir[p-first] = first[p-first];
  
  global_path = strdup(dir);
}

