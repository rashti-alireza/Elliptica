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
  const int test_det_adm_g/* if 1, it tests det(adm_g) > 0 */
  );

#endif


