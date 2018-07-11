#include "core_lib.h"
#include "error_handling_lib.h"

void free_2d(void *mem0);
void free_matrix(void *mem0, const unsigned long c);
void free_needle(Needle_T *needle);
void free_points(Grid_T *const grid);
void free_func_PtoV(sFunc_PtoV_T **func);
void free_field(Field_T *f);
