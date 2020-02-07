#include "core_lib.h"
#include "error_handling_lib.h"

int make_parameters(const char *const path);
int global_variables_init(const char *const path);
int projects_data_base(void);
void *get_project(const char *const proj_name);
int project_execute(Project_T *const proj);
void free_db_project(void);
static void pr_logo(void);
static void pr_info(void);
