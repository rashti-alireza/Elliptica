#include "core_lib.h"
#include "error_handling_lib.h"

extern char *inputfile_name_global;

Parameter_T *get_parameter(char *const par_name);
static void set_default(char *lhs,char *rhs);
