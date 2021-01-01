#ifndef physics_blackhole_LIB_H
#define physics_blackhole_LIB_H
#include "elliptica_system_lib.h"

/* forward declaration */
struct PHYSICS_T;

int bh_main(struct PHYSICS_T *const phys);

void 
bh_interpolating_fields_on_a_line
  (
  Physics_T *const phys/* physics of interest */,
  const char *const sfields_name/* comma separated fields */,
  const char *const dir/* output directory */,
  const char *const stem_g/* if stem of a metric given => test det(g) > 0 */
  );
  
void 
bh_print_properties
  (Physics_T *const phys,
  const char *const params,
  FILE *const file,
  const int pr_screen);
  
double bh_calculate_expansion_on_AH(Physics_T *const phys);

#endif


