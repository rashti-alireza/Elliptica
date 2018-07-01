#include "core_lib.h"
#include "error_handling_lib.h"

extern char *inputfile_name_global;

Parameter_T *get_parameter(const char *const par_name);
static void set_default(const char *const lhs,const char *const rhs);
void add_parameter(const char *const lv, const char *const rv);
void set_default_parameter(void);
