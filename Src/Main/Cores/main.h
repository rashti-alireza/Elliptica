#include "core_lib.h"
#include "error_handling_lib.h"

int make_parameters(char *const path);
int global_variables_init(char *const path);
int projects_data_base(void);
void *get_project(char *const proj_name);
int project_execute(Project_T *proj);