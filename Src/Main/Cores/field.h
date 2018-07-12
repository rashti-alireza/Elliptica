#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "approximation_lib.h"

Field_T *init_field(const char *const name,Grid_T *const grid);
void add_field(Field_T *const f,Grid_T *const grid);
Field_T *get_field_S(const char *const name,Grid_T *const grid);
double *make_coeffs(Field_T *const f);