void make_jacobian_eq(Grid_T *const grid, const char * const* types);
void test_make_jacobian_eq(Grid_T *const grid, const char * const* types);
void *init_eq(void);
void add_eq(sEquation_T ***const data_base, fEquation_T *const eq,const char *const name);
void populate_solution_man(Grid_T *const grid,sEquation_T **const field_eq,sEquation_T **const bc_eq);

