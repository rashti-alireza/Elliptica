#include "core_lib.h"
#include "maths_approximation_lib.h"
#include "maths_analytic_lib.h"

int Laplace_Inhom(void);
Grid_T *Laplace_Inhom_make_grid(void);
int Laplace_Inhom_solve_eq(Grid_T *const gird);
int Laplace_Inhom_pr_answer(Grid_T *const grid);
int Laplace_Inhom_clean_up(Grid_T *const grid);
