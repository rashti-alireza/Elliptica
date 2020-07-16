#include "core_lib.h"
#include "error_handling_lib.h"
#include <umfpack.h>
#include "maths_linear_algebra_lib.h"
#include "maths_matrix_solvers_lib.h"

int direct_solver_umfpack_di(void *vp);
int direct_solver_umfpack_dl(void *vp);
void umfpack_error_di(const double *const Control,const int status,const char *const file,const int line);
void umfpack_error_dl(const double *const Control,const long status,const char *const file,const int line);
static void umfpack_failed(const int status,const char *const file,const int line);
int direct_solver_series_umfpack_di(void *vp);
int direct_solver_series_umfpack_dl(void *vp);
Umfpack_T *init_umfpack(void);
void free_umfpack(Umfpack_T *umf);






