#include "core_lib.h"
#include "error_handling_lib.h"
#include <suitesparse/umfpack.h>

int direct_solver_umfpack_di(void *vp);
int direct_solver_umfpack_dl(void *vp);
void umfpack_error_di(const double *const Control,const int status,const char *const file,const int line);
void umfpack_error_dl(const double *const Control,const long status,const char *const file,const int line);
static void umfpack_failed(const int status,const char *const file,const int line);
int direct_solver_series_umfpack_di(void *vp);



