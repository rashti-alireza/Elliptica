#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "text_tools_lib.h"
#include "macros_lib.h"
#include "maths_general_lib.h"
#include "coordinates_lib.h"
#include <silo.h>


void pr_fields(const Grid_T *const grid);
static void get_field_vs_coords_3d(const char *const par,char ***const field_name,char ***const c0,char ***const c1,char ***const c2,unsigned *const nf);
