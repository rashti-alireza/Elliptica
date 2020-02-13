#include "core_lib.h"
#include "error_handling_lib.h"

int make_parameters(const char *const path);
int init_global_variables(const char *const path);
int create_db_projects(void);
void *get_project(const char *const proj_name);
int execute_project(Project_T *const proj);
void free_project_db(void);
void free_parameter_db(void);
void free_grid_db(void);
static void pr_logo(void);
static void pr_info(void);
