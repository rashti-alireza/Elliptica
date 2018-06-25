/*
// Alireza Rashti
// May 2018
*/

#include "projects.h"

/* executing the project */
int project_execute(Project_T *proj)
{
  int ret;
  
  /* printing project name */
  pr_comment(proj->name);
  
  /* calling projectc func */
  ret = proj->func();
  
  return ret;
}


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
Project_T *get_project(char *const proj_name)
{
  int i;
  
  i = 0;
  while (projects_global != 0 && projects_global[i] != 0)
  {
    if (strcmp_i(projects_global[i]->name,proj_name))
      return projects_global[i];
    i++;
  }
  
  return 0;
}

