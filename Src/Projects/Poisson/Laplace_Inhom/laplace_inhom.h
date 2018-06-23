#include "core_lib.h"

Grid_T *Laplace_Inhom_make_grid(void);
int Laplace_Inhom_solve_eq(Grid_T *gird);
int Laplace_Inhom_pr_answer(Grid_T *grid);
int Laplace_Inhom_clean_up(Grid_T *grid);
