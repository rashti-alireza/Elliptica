#include "core_lib.h"
#include "maths_approximation_lib.h"
#include "maths_special_functions_lib.h"

#define STR_LEN_MAX 400

int Poisson0(void);
Grid_T *poisson0_make_grid(void);
int poisson0_solve_eq(Grid_T *const grid);
int poisson0_pr_answer(Grid_T *const grid);
int poisson0_clean_up(Grid_T *const grid);
int poisson0_analyze_answer(const Grid_T *const grid);
