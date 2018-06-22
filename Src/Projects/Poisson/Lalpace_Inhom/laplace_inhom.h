#include "core_lib.h"

Grid_T *laplace_inhom_make_grid(void);
int laplace_inhom_solve_eq(Grid_T *gird);
int laplace_inhom_pr_answer(Grid_T *grid);
int laplace_inhom_clean_up(Grid_T *grid);
