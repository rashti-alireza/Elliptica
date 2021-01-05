#include "star_header.h"

int star_main(Physics_T *const phys);
static int tune_star_Euler_constant(Physics_T *const phys);
static int extrapolate_star_matter(Physics_T *const phys);
static int find_star_surface(Physics_T *const phys);
static int tune_star_force_balance_equation(Physics_T *const phys);
static int tune_star_center(Physics_T *const phys);
static int add_star_fields(Physics_T *const phys);
static int set_star_params(Physics_T *const phys);
static int start_off_star(Physics_T *const phys);

