/*
// Alireza Rashti
// May 2018
*/

#include "projects.h"

/* executing the project */
int project_execute(Project_T *const proj)
{
  int ret;
  
  /* printing project name */
  pr_comment(proj->name);
  
  /* calling projectc func */
  ret = proj->func();
  
  /* printing project name */
  pr_comment(proj->name);
  
  return ret;
}


/* adding the desired project for run */
void add_project(ProjFunc *const projfunc, const char *const name, const char *const des)
{
  if (projfunc == 0)
   abortEr("Empty project.\n");
  
  Project_T *proj;

  proj = get_project(name);
  if (proj)
    abortEr_s("This project \"%s\" has already been added!",name);
  
  proj = alloc_project(&projects_global);
  
  proj->name = dup_s(name);
  proj->des  = dup_s(des);
  proj->func = projfunc;

}

/* having project name, it returns a pointer to 
// the corresponding project func
*/
Project_T *get_project(const char *const proj_name)
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

