#ifndef physics_equation_LIB_H
#define physics_equation_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PHYSICS_T;

/* analyzing solution after solve */
typedef void fFunc_analyze_solution_T(struct PHYSICS_T *const phys,const int iteration);

int eq_main(struct PHYSICS_T *const phys);


#endif


