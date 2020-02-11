#include "bbn_headers.h"

int Binary_BH_NS_Initial_Data(void);
Grid_T *bbn_initialize_next_grid(Grid_T *const grid_prev);
void bbn_solve_elliptic_eqs(Grid_T *const grid);
void bbn_free_grid_and_its_parameters(Grid_T *grid);
void bbn_set_default_parameters(void);
void bbn_write_checkpoint(const Grid_T *const grid);
static void Elliptic_Eqs_Convergence_Test_BBN(void);
static void update_parameters_and_directories(const unsigned iter);



