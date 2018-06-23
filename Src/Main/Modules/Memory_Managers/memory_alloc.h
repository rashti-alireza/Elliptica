#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"

static void alloc_patches_cartesian(Grid_T *grid);
static void alloc_nodes_cartesian(Grid_T *grid);
Parameter_T *get_parameter(char *const par_name);
