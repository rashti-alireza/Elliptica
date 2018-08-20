#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"

#define FACE_NUM 6

extern Grid_T **grids_global;

void *alloc_parameter(Parameter_T ***const mem);
void *alloc_project(Project_T ***const mem);
void *alloc_grid(void);
void alloc_patches(Grid_T *const grid);
void alloc_nodes(Grid_T *const grid);
void alloc_interface(Patch_T *const patch);
void *alloc_point(const unsigned s);
void *alloc_sFunc_PtoV(sFunc_PtoV_T ***const mem);
void *alloc_needle(void);
static void alloc_patches_Cartesian_grid(Grid_T *const grid);
Parameter_T *get_parameter(const char *const par_name);
double *alloc_double(const unsigned N);
void *alloc_sFunc_Patch2Pdouble(sFunc_Patch2Pdouble_T ***const mem);
double **alloc_matrix(const long unsigned R,const long unsigned C);
Solution_T *alloc_solution(Patch_T *const patch, const unsigned n);
