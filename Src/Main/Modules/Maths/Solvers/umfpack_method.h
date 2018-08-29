#include "core_lib.h"
#include "error_handling_lib.h"
#include <suitesparse/umfpack.h>

int direct_solver_umfpack_di(void *vp);
int direct_solver_umfpack_dl(void *vp);
static void umfpack_error(const int status,const char *const file,const int line);

