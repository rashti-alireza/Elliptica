#include "sbh_headers.h"
#include "maths_field_analysis.h"

void sbh_study_initial_data(Grid_T *const grid);
void sbh_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder);
void sbh_convergence_test(const Grid_T *const grid);
