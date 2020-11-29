#include "bh_header.h"

int bh_main(Physics_T *const phys);
static int start_off_black_hole(Physics_T *const phys);
static int add_black_hole_params(Physics_T *const phys);
static int add_black_hole_fields(Physics_T *const phys);
static int tune_black_hole_radius(Physics_T *const phys);
static int find_black_hole_surface(Physics_T *const phys);

