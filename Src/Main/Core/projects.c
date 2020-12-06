/*
// Alireza Rashti
// May 2018
*/

#include "projects.h"

/* executing the project */
int execute_project(Project_T *const proj)
{
  int ret;
  
  /* printing project name */
  pr_comment(proj->name);
  
  /* calling projectc func */
  ret = proj->func(0);
  
  /* printing project name */
  pr_comment(proj->name);
  
  return ret;
}


/* adding the desired project for run */
void add_project(ProjFunc *const projfunc, const char *const name, const char *const des)
{
  if (projfunc == 0)
   Error0("Empty project.\n");
  
  Project_T *proj;

  proj = get_project(name);
  if (proj)
    Errors("This project \"%s\" has already been added!",name);
  
  proj = alloc_project(&projects_global);
  
  proj->name = dup_s(name);
  proj->des  = dup_s(des);
  proj->func = projfunc;

}

/* free the whole data base of project */
void free_project_db(void)
{
  Uint np;
  
  np = 0;
  while (projects_global != 0 && projects_global[np] != 0)
  {
    Free(projects_global[np]->name);
    Free(projects_global[np]->des);
    free(projects_global[np]);
    np++;
  }
  Free(projects_global);
  
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

/* adding 2 block of memory for project data base 
// and putting the last block to null and 
// returning pointer to one before the last block
*/
void *alloc_project(Project_T ***const mem)
{
  Uint i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  IsNull((*mem));
  
  (*mem)[i] = malloc(sizeof(*(*mem)[i]));
  IsNull((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}



