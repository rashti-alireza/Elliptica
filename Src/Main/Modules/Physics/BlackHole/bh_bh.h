#include "bh_header.h"

int bh_tune_black_hole_radius(Physics_T *const phys);
int bh_find_black_hole_surface(Physics_T *const phys);
int bh_fill_inside_black_hole(Physics_T *const phys);
static void tune_BH_radius_irreducible_mass_perfect_s2(Physics_T *const phys);
static void find_bh_surface_perfect_s2(Physics_T *const phys);
static void start_off_KerrSchild_perfect_s2(Physics_T *const phys);
int bh_start_off(Physics_T *const phys);
int bh_add_params(Physics_T *const phys);
int bh_add_fields(Physics_T *const phys);




