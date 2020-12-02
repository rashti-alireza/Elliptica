#ifndef maths_complex_LIB_H
#define maths_complex_LIB_H
#include "elliptica_system_lib.h"


#include <complex.h>
/* #include <tgmath.h> for the type generate macros. */
#undef I
#define imagI _Complex_I

void *alloc_double_complex(const Uint N);

#endif




