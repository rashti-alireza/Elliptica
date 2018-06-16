/*
// Alireza Rashti
// May 2018
*/

#include "projects.h"

/* adding the desired project for run */
void add_project(ProjFunc *projfunc, char *const name, char *const des)
{
  pointerEr(projfunc);
  
  Project_T *proj;

  proj = get_project(name);
  if (proj)
    abortEr_s("This project \"%s\" has already been added!",name);
  
  proj = alloc_project(&projects_global);
  
  proj->name = strdup(name);
  proj->des  = strdup(des);
  proj->func = projfunc;

}

/* having project name, it returns a pointer to 
// the corresponding project func
*/
void *get_project(char *const proj_name)
{
  int i;
  
  i = 0;
  while (projects_global != 0 && projects_global[i] != 0)
  {
    if (strcmp(projects_global[i]->name,proj_name) == 0)
      return projects_global[i]->func;
    i++;
  }
  
  return 0;
}

