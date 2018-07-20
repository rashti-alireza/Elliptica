#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"

Variable_T *add_variable(const char *const name,const char *const attribute,Patch_T *const patch,const Flag_T alloc_flg);
void remove_variable(const char *const name,Patch_T *const patch);
void add_attribute(Variable_T *const var,const char *const attribute);
int LookUpVar(const char *const name,Patch_T *const patch);


