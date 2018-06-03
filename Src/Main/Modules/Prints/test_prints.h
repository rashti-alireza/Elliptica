#include "core_lib.h"
#include "error_handling_lib.h"

extern char *global_path;
extern Parameter_T **global_parameter;

void *get_parameter_value(char *const par_name,Flag_T kind, double *value);