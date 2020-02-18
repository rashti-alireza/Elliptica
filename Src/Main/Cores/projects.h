#include "core_lib.h"
#include "error_handling_lib.h"

extern Project_T **projects_global;

Project_T *get_project(const char *const proj_name);
int execute_project(Project_T *const proj);
void add_project(ProjFunc *const projfunc, const char *const name, const char *const des);
void free_project_db(void);
void *alloc_project(Project_T ***const mem);
