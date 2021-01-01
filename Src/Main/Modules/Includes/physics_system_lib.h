#ifndef physics_system_LIB_H
#define physics_system_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PHYSICS_T;

int sys_main(struct PHYSICS_T *const phys);

void 
sys_print_properties
  (Physics_T *const phys,
  const char *const params,
  FILE *const file,
  const int pr_screen);

#endif


