#include "core_lib.h"

Grid_T *add_grid(void);
static int make_patches(Grid_T *grid);
static int make_coords(Grid_T *grid);
static int realize_geometry(Grid *grid);
