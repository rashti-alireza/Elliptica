#include "bbn_headers.h"
#include "physics_observables_lib.h"
#include <unistd.h>


void bbn_study_initial_data(Grid_T *const grid);
void bbn_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder);
void bbn_print_residual_norms(Grid_T *const grid,const unsigned iteration, const char *const folder);
