#include "bh_header.h"

int bh_main(Physics_T *const phys);
static int start_off_black_hole(Physics_T *const phys);
static int set_black_hole_params(Physics_T *const phys);
static int add_black_hole_fields(Physics_T *const phys);
static int tune_black_hole_radius(Physics_T *const phys);
static int tune_black_hole_spin(Physics_T *const phys);
static int find_black_hole_surface(Physics_T *const phys);
static int update_conformal_normal_vector_on_AH(Physics_T *const phys);
static int update_inner_BC_values(Physics_T *const phys);


