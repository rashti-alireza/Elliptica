#ifndef maths_complex_LIB_H
#define maths_complex_LIB_H


#include <complex.h>
/* #include <tgmath.h> for the type generate macros. */
#undef I
#define imagI _Complex_I

void *alloc_double_complex(const unsigned N);

#endif



