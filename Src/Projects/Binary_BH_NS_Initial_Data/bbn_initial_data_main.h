#include "bbn_headers.h"

#define STR_LEN_MAX 400

Grid_T *bbn_initialize_next_grid(Grid_T *const grid_prev);
void bbn_solve_elliptic_eqs(Grid_T *const grid);
void bbn_free_grid_and_its_parameters(Grid_T *grid);
void bbn_set_default_parameters(void);
void bbn_write_checkpoint(const Grid_T *const grid);
void bbn_construct_id(void);
void bbn_elliptic_eqs_convergence_test(void);
static void update_parameters_and_directories(const unsigned main_loop_iter);



