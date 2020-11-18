#include "core_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "fields_lib.h"

#define STR_SIZE1 (99)
#define STR_SIZE2 (999)

void 
update_parameters_and_directories
  (
   const unsigned main_loop_iter,
   const char *const dir_name_format/* eg: "BBN_%s_%ux%ux%u" */
  );

void free_grid_and_its_parameters(Grid_T *grid);

